# Install packages 
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("IlluminaHumanMethylation450kanno.ilmn12.hg19", dependencies=TRUE)
BiocManager::install("IlluminaHumanMethylationEPICanno.ilm10b4.hg19",dependencies=TRUE)
BiocManager::install(c("RnBeads","limma"))
install.packages("readxl")
install.packages("tidyverse")
install.packages("preprocessCore")
install.packages("minfi")  

# Load the libraries
library(readxl)
library(tidyverse)
library(RnBeads)
library(preprocessCore)
library(minfi)
library(limma)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

#Load the Dataset
file_path <- "C:/Users/tulasi/Downloads/colorectal cancer - project/GSE_220160.xlsx" 
methylation_data <- read_excel(file_path)


#Data preprocessing
# remove missing values
methylation_data <- na.omit(as.data.frame(methylation_data))

# Select only numeric columns (assuming non-numeric columns are first 4)
numeric_data <- methylation_data[, 5:ncol(methylation_data)]  # Exclude first 4 columns
str(methylation_data)

# Convert to a numeric matrix
beta_values <- as.matrix(numeric_data)
mode(beta_values) <- "numeric"

if (!is.numeric(beta_values)) stop("Beta values are not numeric")

# Apply quantile normalization
normalized_data <- normalize.quantiles(beta_values)

# Convert back to a data frame
normalized_df <- as.data.frame(normalized_data)

# Restore column names (samples) and row names (CpG sites)
colnames(normalized_df) <- colnames(numeric_data)
rownames(normalized_df) <- methylation_data[[1]]  # Restore probe IDs
normalized_df <- normalized_df[!duplicated(rownames(normalized_df)), ]


#Filter Out CpGs on Allosomes (X & Y Chromosomes)
# Check chromosomes in the annotation
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
annotation <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

table(annotation$chr)

# Convert chromosome column to character
annotation$chr <- as.character(annotation$chr)  

# Extract autosomal CpGs
autosomal_cpgs <- rownames(annotation)[annotation$chr %in% paste0("chr", 1:22)]
length(autosomal_cpgs)  # Ensure this is > 0

# Ensure beta_values has correct row names (CpG IDs)
rownames(beta_values) <- methylation_data[[1]]  

# Check the row names to verify they are now correct CpG IDs
head(rownames(beta_values))

# Check if any CpGs in beta_values match autosomal CpGs
matching_cpgs <- intersect(rownames(beta_values), autosomal_cpgs)
length(matching_cpgs)  # Should now be > 0 if the filtering works

# Filter the beta_values to include only autosomal CpGs
beta_values_filtered <- beta_values[matching_cpgs, ]

# Check dimensions after filtering
dim(beta_values_filtered)  # Should show a reduced number of rows (not 0)


#Differentially Methylation Analysis 
# Manually define sample groups
sample_groups <- c("cancer", "cancer", "cancer", "cancer", "cancer", "cancer", "cancer", "cancer",
                   "polyp", "polyp", "polyp", "polyp", "polyp", "polyp", "polyp", "polyp")

sample_groups <- factor(sample_groups, levels = c("polyp", "cancer"))  
table(sample_groups)

# Convert Beta Values to M-values (logit transformation)
m_values <- log2((beta_values + 0.001) / (1 - beta_values + 0.001))

# Create Design Matrix
design <- model.matrix(~ sample_groups)

# Fit Linear Model
fit <- lmFit(m_values, design)
fit <- eBayes(fit)

# Extract Differentially Methylated Probes (DMPs)
dmps <- topTable(fit, coef = 2, adjust.method = "fdr", number = Inf)

# Apply Delta Beta & FDR Thresholds
dmps_filtered <- dmps[abs(dmps$logFC) > 0.1 & dmps$adj.P.Val < 0.6, ]

hist(dmps$adj.P.Val, breaks = 50, main = "Distribution of Adjusted P-Values", xlab = "Adjusted P-Value")

# Load gene annotation data
probe_gene_map <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)[, c("Name", "UCSC_RefGene_Name")]

# Merge with methylation data
df_annotated <- merge(methylation_data, probe_gene_map, by.x = "ID_REF", by.y = "Name", all.x = TRUE)

# Split multiple gene names into separate rows
df_annotated$UCSC_RefGene_Name <- strsplit(df_annotated$UCSC_RefGene_Name, ";")
df_annotated <- data.frame(
  ID_REF = rep(df_annotated$ID_REF, sapply(df_annotated$UCSC_RefGene_Name, length)),
  Gene = unlist(df_annotated$UCSC_RefGene_Name)
)

df_annotated$Gene <- trimws(toupper(df_annotated$Gene))
df_annotated <- distinct(df_annotated)

# Check if significant genes are available
if (!exists("df_annotated") || nrow(df_annotated) == 0) {
  stop("Error: df_annotated is empty. Ensure gene annotation is correct.")
}

if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
  install.packages("BiocManager")
  BiocManager::install("clusterProfiler")
}
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  BiocManager::install("org.Hs.eg.db")
}

# Load packages
library(clusterProfiler)
library(org.Hs.eg.db)

significant_genes <- unique(df_annotated$Gene[df_annotated$ID_REF %in% rownames(dmps_filtered)])
entrez_ids <- bitr(significant_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

if (nrow(entrez_ids) == 0) {
  stop("??? Warning: No valid gene mappings found. Check input gene list.")
}

unmapped_genes <- significant_genes[!significant_genes %in% entrez_ids$SYMBOL]
print(unmapped_genes)
significant_genes <- toupper(trimws(significant_genes))  # Convert to uppercase and remove spaces

filtered_entrez_gene_list <- na.omit(entrez_ids$ENTREZID)  # Remove NA values
filtered_entrez_gene_list <- unique(filtered_entrez_gene_list)  # Remove duplicates

#Pathway analysis
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(biomaRt)

kegg_results <- enrichKEGG(
  gene = filtered_entrez_gene_list,  # Your list of Entrez IDs
  organism = "hsa",  # hsa for Homo sapiens
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)

# GO enrichment
BiocManager::install(c("DOSE"))
go_results <- enrichGO(
  gene = filtered_entrez_gene_list,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "ALL",
  pAdjustMethod = "fdr",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)

# Print the structure
str(go_results)

# Check the number of significant terms found
if (is.null(go_results) || nrow(go_results@result) == 0) {
  print("??? Warning: No significant GO terms found.")
} else {
  head(go_results@result)  # Display top GO terms
}

# Visualization
barplot(go_results, showCategory = 10, title = "GO Enrichment Analysis")
go_results <- pairwise_termsim(go_results) 
treeplot(go_results)


#Marker analysis
install.packages("devtools", dependencies = TRUE)
install.packages("caret", dependencies = TRUE)
install.packages("randomForest", dependencies = TRUE)

library(devtools)
library(caret)
library(randomForest)

# Prepare data for classification
# Convert methylation beta values to a data frame with sample labels
df_classification <- as.data.frame(t(beta_values_filtered))  # Transpose so samples are rows
df_classification$Group <- sample_groups  # Add group labels

# Ensure the class labels are factors
df_classification$Group <- as.factor(df_classification$Group)

# Create training and test sets
set.seed(123)  # For reproducibility
train_index <- createDataPartition(df_classification$Group, p = 0.8, list = FALSE)
train_data <- df_classification[train_index, ]
test_data <- df_classification[-train_index, ]

# Define RFE control
ctrl <- rfeControl(functions = rfFuncs, method = "cv", number = 5)

# Run Recursive Feature Elimination
rfe_result <- rfe(train_data[, -ncol(train_data)], train_data$Group,
                  sizes = c(10, 20, 50, 100),  # Number of features to try
                  rfeControl = ctrl)

# Get top selected features (biomarkers)
selected_features <- predictors(rfe_result)
print(selected_features)

# Train Random Forest classifier on selected markers
rf_model <- randomForest(Group ~ ., data = train_data[, c(selected_features, "Group")], importance = TRUE)

# Predict on test data
rf_predictions <- predict(rf_model, newdata = test_data[, selected_features])

# Evaluate performance
conf_matrix <- confusionMatrix(rf_predictions, test_data$Group)
print(conf_matrix)

# Plot variable importance
varImpPlot(rf_model, main = "Top Differentially Methylated Markers")




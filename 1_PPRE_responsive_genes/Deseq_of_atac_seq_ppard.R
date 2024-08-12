# Install packages if not already installed
# install.packages("DESeq2")
# install.packages("ggplot2")
# install.packages("dplyr")

# Load required libraries
library(DESeq2)
library(ggplot2)
library(dplyr)

# Check dimensions of count data matrix
ncol_countData <- ncol(atac_ppard_pca)
length_sampleTypes <- length(sample_names)

# Compare dimensions
print(ncol_countData)
print(length_sampleTypes)

# Round count data matrix to integers
atac_ppard_pca_int <- round(atac_ppard_pca)

# Create a DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = as.matrix(atac_ppard_pca_int),
                              colData = DataFrame(sample_type = sample_names),
                              design = ~ sample_type)

# Perform variance stabilization transformation
dds <- DESeq(dds)

# Get normalized counts (regularized log transformation)
normalized_counts <- assay(rlog(dds))


# Perform PCA on normalized counts
pca_result <- prcomp(normalized_counts, scale. = TRUE)

# Extract PCA results (PC1 and PC2)
pca_df <- as.data.frame(pca_result$x[, 1:2])
pca_df_cols <- tibble::rownames_to_column(pca_df, "peak_id")

atac_peak_id_ensembl_id <- atac_genes %>% select(peak_id, ensembl_gene_id, sample_type)

pca_df_samples <- merge(pca_df_cols, atac_peak_id_ensembl_id, by = "peak_id")

pca_df$Sample_Type <- sample_names # Add sample type information to PCA results
rownames(pca_df) <- colnames(normalized_counts)  # Use sample names as row names

pca_df_names <- pca_df %>%
  mutate(
    modified_sample_names = gsub("_\\d+$", "", Sample_Type)
  )

# Create PCA plot with color mapping based on Sample_Type
ggplot(pca_df_cols, aes(x = PC1, y = PC2)) +
  geom_point(size = 3) +
  scale_color_viridis_d() +
  labs(title = "PCA Plot",
       x = "PC1",
       y = "PC2",
       color = "Sample Type") +
  theme_minimal()


# Convert sample types and counts to data frame
data_df <- data.frame(sample_name = sample_names, t(normalized_counts)) %>%
  mutate(
    sample_name = gsub("_\\d+$", "", sample_name)
  )

# Calculate variance and identify top variable genes per sample type
top_var_genes_per_type <- data_df %>%
  group_by(sample_name) %>%
  summarize(across(starts_with("ID_"), var)) %>%
  pivot_longer(cols = starts_with("ID_"), names_to = "peak_id", values_to = "Variance") %>%
  group_by(sample_name) %>%
  top_n(10, Variance) %>%
  arrange(sample_name, desc(Variance))



ppard_atac_peak_ensembl <- atac_ppard %>% select(ensembl_gene_id, peak_id, hgnc_symbol)

top_var_genes_per_type_gene_name <- merge(top_var_genes_per_type, ppard_atac_peak_ensembl, by = "peak_id", all.x = TRUE)

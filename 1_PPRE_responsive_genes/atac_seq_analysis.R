library(tidyverse)

atac_ppard <- read_delim("~/Documents/COURSES/UU/MRP/GSE178984_202010625_dubois_atac_normcounts.txt/GSE178984_202010625_dubois_atac_normcounts.txt", 
                         delim = "\t", col_names = TRUE) %>%
  pivot_longer(cols = 7:24, names_to = "sample_type", values_to = "norm_count")
  

atac_genes <- atac_ppard %>%
  mutate(
    base_sample_type = str_extract(sample_type, "^[A-Za-z]+"),
    replicate = str_extract(sample_type, "\\d+$")
  ) %>%
  group_by(ensembl_gene_id, base_sample_type) %>%
  summarize(
    mean_norm_count = mean(norm_count)
  )

## PCA plot

pca_atac <- atac_genes %>%
  select(ensembl_gene_id, mean_norm_count, base_sample_type) %>%
  spread(base_sample_type, mean_norm_count) %>%
  column_to_rownames(var = "ensembl_gene_id")

pca_atac_scaled <- scale(pca_atac)

pca_atac_result <- prcomp(pca_atac_scaled, scale. = TRUE)

pca_atac_df <- as.data.frame(pca_atac_result$x[, 1:2])  # Keep only the first two principal components

# Add row names as a column (gene_id)
pca_atac_df$ensembl_gene_id <- rownames(pca_atac)

# Rename columns
colnames(pca_atac_df) <- c("PC1", "PC2", "ensembl_gene_id")

pca_atac_df <- merge(pca_atac_df, atac_genes, by = "ensembl_gene_id", all.x = TRUE)

ggplot(pca_atac_df, aes(x = PC1, y = PC2, color = ensembl_gene_id)) +
  geom_point() +
  labs(title = "PCA Plot",
       x = paste("PC1", pca_atac_result$eig[1, 2], "% variance", sep = ""),
       y = paste("PC2", pca_atac_result$eig[2, 2], "% variance", sep = "")) +
  theme_minimal()

ggplot(atac_genes, aes(x = hgnc_symbol, y = mean_norm_count)) +
  geom_point()

atac_summary <- atac_ppard %>%  
  group_by(sample_type) %>%
  summarize(
    n = n(),
    norm_count = mean(norm_count)
    
  )

control <- c("Control_1", "Control_2", "Control_3")
LCFA <- c("LCFA_1", "LCFA_2", "LCFA_3")

atac_control <- atac_ppard %>%
  filter(sample_type %in% control)

atac_control

atac_LCFA <- atac_ppard %>%
  filter(sample_type %in% LCFA)


ggplot(atac_ppard, aes(ensembl_gene_id, norm_count)) +
  geom_boxplot() +
  facet_grid(~sample_type) 

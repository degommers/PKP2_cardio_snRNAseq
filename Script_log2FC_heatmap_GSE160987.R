# Packages
library(tidyverse)
library(pheatmap)
library(tsne)
library(reshape2)
library(scales)
library(DESeq2)

# Read and preprocess gene expression data
GSE160987_20201015_Dubois_counts <- read.csv("~/Documents/COURSES/UU/MRP/GSE178984_GSE160987_dubois_atac_normcounts.txt/GSE160987_20201015_Dubois_counts.csv")

# Filter out invalid gene symbols
GSE160987_20201015_Dubois_counts_rm <- GSE160987_20201015_Dubois_counts %>%
  filter(!grepl("\\d{1,2}-[A-Za-z]{3}$", gene_symbol) & str_detect(gene_symbol, "[A-Z]"))

# Set row names and convert to matrix
GSE160987_matrix <- as.matrix(GSE160987_20201015_Dubois_counts_rm[, -1])
rownames(GSE160987_matrix) <- GSE160987_20201015_Dubois_counts_rm[, 1]

# Select samples for 28 days
GSE160987_matrix_28days <- GSE160987_matrix[, 1:30]

# Replace NA values with zero (adjust as needed)
GSE160987_matrix_28days[is.na(GSE160987_matrix_28days)] <- 0

# Read and preprocess sample metadata
col_data_GSE160987 <- read.delim("~/Documents/COURSES/UU/MRP/GSE178984_GSE160987_dubois_atac_normcounts.txt/col_data_RNAseq_PPARd_paper.txt") %>%
  mutate(
    sample_name = str_replace_all(sample_name, "[^[:alnum:]]", "."),
    run = str_replace_all(run, "[^[:alnum:]]", "."),
    treatment = str_replace_all(treatment, "[^[:alnum:]]", "_")
  )

# Set factors for treatment and run columns
col_data_GSE160987$treatment <- factor(col_data_GSE160987$treatment)
col_data_GSE160987$run <- factor(col_data_GSE160987$run)

# Prepare sample metadata for DESeq2
col_data_GSE160987_rownames <- col_data_GSE160987 %>%
  select(-1) %>%
  as.data.frame()

# Set row names using sample names
rownames(col_data_GSE160987_rownames) <- col_data_GSE160987$sample_name



# Create DESeqDataSet object
GSE160987_dds <- DESeqDataSetFromMatrix(
  countData = GSE160987_matrix_28days,
  colData = col_data_GSE160987_rownames,
  design = ~ treatment
)

# Set reference level for treatment
GSE160987_dds$treatment <- relevel(GSE160987_dds$treatment, ref = "Control")

# Perform DESeq2 analysis
dds <- DESeq(GSE160987_dds)
contrasts_to_examine <- resultsNames(dds) 


# Initialize an empty list to store results
results_list <- list()


log2foldchange_list <- list()
# Loop through each contrast in contrasts_to_examine
for (contrast_name in contrasts_to_examine[2:9]) {
  print(contrast_name)
  
  factorname <- str_extract(contrast_name, "^(treatment)")
  numerator <- str_extract(contrast_name, "(LCFA___PPAR[agd]_Ag|LCFA|PPAR[agd]_A[ng])")
  denominator <- str_extract(contrast_name, "(Control)$")
  
  print(factorname)
  print(numerator)
  print(denominator)
  # Extract results using the specified contrast
  contrast_results <- results(dds, contrast = c(factorname, numerator, denominator))
  
  # Extract Log2FoldChange values and store in a matrix
  log2foldchange_matrix <- matrix(contrast_results$log2FoldChange,
                                  nrow = length(contrast_results$log2FoldChange),
                                  ncol = 1,
                                  dimnames = list(rownames(contrast_results), NULL))
  
  # Assign column name to the Log2FoldChange matrix
  colnames(log2foldchange_matrix) <- contrast_name
  
  # Store the Log2FoldChange matrix in the list
  log2foldchange_list[[contrast_name]] <- log2foldchange_matrix
}

# Combine Log2FoldChange matrices into a single data frame
log2foldchange_df <- do.call(cbind, log2foldchange_list)
log2foldchange_df_df <- as.data.frame(log2foldchange_df)

log2foldchange_df_1 <- log2foldchange_df_df %>% select(treatment_PPARa_Ag_vs_Control) %>% filter(treatment_PPARa_Ag_vs_Control > 0.5)
log2foldchange_df_1$gene_name <- rownames(log2foldchange_df_1)
log2foldchange_genes <- log2foldchange_df_1 %>%select(gene_name) %>% write_csv("logfoldchange_upregulated_PPARa_0_5.csv", col_names = FALSE)

matrix_LFC <- as.matrix(log2foldchange_df, rownames.force = TRUE, colnames.force = TRUE)

genes_of_interest <- c("ABCB11", "ACACB",
                       "SLC25A17",
                       "DECR1",
                       "CPT1A",
                       "EHHADH",
                       "ACADS",
                       "SESN2",
                       "HADH",
                       "SLC27A2",
                       "FABP1",
                       "ETFDH",
                       "MTLN",
                       "ADIPOQ",
                       "PPARA",
                       "ACAD11",
                       "DECR2",
                       "FABP3",
                       "PPARG",
                       "PPARD",
                       "LEP",
                       "RXRA",
                       "RXRB",
                       "RXRG"
                       )

officialgeneset <- c("ABCB11", "ACACB",
                      "SLC25A17",
                      "DECR1",
                      "CPT1A",
                      "EHHADH",
                      "ACADS",
                      "SESN2",
                      "HADH",
                      "SLC27A2",
                      "FABP1",
                      "ETFDH",
                      "MTLN",
                      "ADIPOQ",
                      "PPARA",
                      "ACAD11",
                      "DECR2",
                      "APOA5",
                      "ABCB11",
                      "FGF21",
                      "LEP",
                      "HMGCS2",
                      "ACSL1",
                      "ANGPTL4",
                      "CYP7A1")

genes_motif_ppre_paper_2010 <- c("PLIN2", "CYP7A1", "CPT2", "FABP1", "CBX1", "FGF21", 
                                 "SLC25A20", "ANGPTL4", "CPT1A", "UGT1A9", "APOA2", 
                                 "HMGCS2", "APOA5", "NR1H2", "ACSL1", "ACOX1", "ACADM")

genes_motif_ppre_motif_not_paper_2010 <- c("ACACB", "ECH1", "ALDH3A1", "APOA1", "ACADS", "ACSL3", 
                                                           "HACL1", "G6PC1", "HSD17B4", "SLC27A2", "PCTP", "ALDH9A1", 
                                                           "VLDLR", "FADS1", "AGPAT2", "ETFDH", "PNPLA2", "GK", "TXNIP")

genes_of_interest2 <- c(genes = c(genes_motif_ppre_paper_2010, genes_motif_ppre_motif_not_paper_2010))
genes_of_interest2_unique <- unique(genes_of_interest2)
matrix_LFC_genes <- subset(matrix_LFC, rownames(matrix_LFC) %in% genes_of_interest)

matrix_LFC_genes_extra <- subset(matrix_LFC, rownames(matrix_LFC) %in% genes_of_interest2_unique)
matrix_LFC_genes_official <- subset(matrix_LFC, rownames(matrix_LFC) %in% officialgeneset)
# Changing the names of the matrix to proper names:
colnames(matrix_LFC_genes) <- c("LCFA", "LCFA+PPARa (Ag)", "LCFA+PPARd (Ag)", 
                                "LCFA+PPARg (Ag)", "PPARa (Ag)", "PPARd (Ag)",
                                "PPARd (An)", "PPARg (Ag)")

colnames(matrix_LFC_genes_extra) <- c("LCFA", "LCFA+PPARa (Ag)", "LCFA+PPARd (Ag)", 
                                "LCFA+PPARg (Ag)", "PPARa (Ag)", "PPARd (Ag)",
                                "PPARd (An)", "PPARg (Ag)")

colnames(matrix_LFC_genes_official) <- c("LCFA", "LCFA+PPARa (Ag)", "LCFA+PPARd (Ag)", 
                                      "LCFA+PPARg (Ag)", "PPARa (Ag)", "PPARd (Ag)",
                                      "PPARd (An)", "PPARg (Ag)")

# Define the column names and their corresponding types
compound_types <- c(
  'LCFA' = 'LCFA',
  'LCFA+PPARa (Ag)' = 'Combined',
  'LCFA+PPARd (Ag)' = 'Combined',
  'LCFA+PPARg (Ag)' = 'Combined',
  'PPARa (Ag)' = 'Agonist',
  'PPARd (Ag)' = 'Agonist',
  'PPARd (An)' = 'Antagonist',
  'PPARg (Ag)' = 'Agonist'
)

receptor_types <- c(
  LCFA = "NA",
  `LCFA+PPARa (Ag)` = "alpha",
  `LCFA+PPARd (Ag)` = "delta",
  `LCFA+PPARg (Ag)` = "gamma",
  `PPARa (Ag)` = "alpha",
  `PPARd (Ag)` = "delta",
  `PPARd (An)` = "delta",
  `PPARg (Ag)` = "gamma"
)

# Define color palettes for compound types and receptor types
compound_colors <- c(
  LCFA = "#1f77b4",         # Blue
  Combined = "#ff7f0e",      # Orange
  Agonist = "#2ca02c",       # Green
  Antagonist = "#d62728"     # Red
)

receptor_colors <- c(
  alpha = "#9467bd",        # Purple
  delta = "#8c564b",        # Brown
  gamma = "#e377c2",        # Pink
  `NA` = "grey"               # Grey for NA (no receptor)
)

# Define the column names from your original matrix
column_names <- c(
  'LCFA',
  'LCFA+PPARa (Ag)',
  'LCFA+PPARd (Ag)',
  'LCFA+PPARg (Ag)',
  'PPARa (Ag)',
  'PPARd (Ag)',
  'PPARd (An)',
  'PPARg (Ag)'
)


# Create a vector to map each column to its corresponding type
compound_annotation <- sapply(column_names, function(col) compound_types[[col]])
receptor_annotation <- sapply(column_names, function(col) receptor_types[[col]])

# Create the annotation DataFrame for columns
annotation_col <- data.frame(
  Treatment = compound_annotation,
  Receptor = receptor_annotation,
  row.names = column_names
)

# Colors for heatmap annotation
# Create lists of colors for compound type and receptor type annotations
compound_color_list <- sapply(compound_annotation, function(comp) compound_colors[[comp]])

receptor_color_list <- sapply(receptor_annotation, function(rec) receptor_colors[[rec]])

# Organize color lists into a nested list structure
annotation_colors <- list(
  Compound_Type_Colors = compound_color_list,
  Receptor_Type_Colors = receptor_color_list
)

# Specify colors
ann_colors = list(
  Receptor = c(  alpha = "#3CAEA3",   # Teal
                 delta = "#FF5733",   # Coral
                 gamma = "#D7E8BA",
                 `NA` = "grey" ),
  Treatment = c(  LCFA = "darkgoldenrod1",         # Blue
                  Combined = "azure4",      # Orange
                  Agonist = "darkolivegreen3",       # Green
                  Antagonist = "darkred" )
)

pheatmap(matrix_LFC_genes_official, scale = "none", show_rownames = TRUE, 
         border_color = NA, fontsize = 12,
         clustering_distance_cols = "euclidean",
         clustering_distance_rows = "euclidean", 
         cutree_cols = 6, 
         cutree_rows = 4, 
         breaks = seq(-4, 4, length.out = 10),
         color = colorRampPalette(c("blue", "white", "red"))(10), 
         annotation_col = annotation_col,
         annotation_colors = ann_colors
         )


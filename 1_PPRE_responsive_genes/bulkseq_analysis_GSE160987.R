library(DESeq2)
library(tidyverse)

GSE160987_20201015_Dubois_counts <- read.csv("GSE160987_20201015_Dubois_counts.csv")

# Remove 'to date' converted gene names and empty gene symbol
GSE160987_20201015_Dubois_counts_rm <- GSE160987_20201015_Dubois_counts %>%
  filter(!grepl("\\d{1,2}-[A-Za-z]{3}$", gene_symbol), str_detect(gene_symbol, "[A-Z]") )

# First column to rownames:
GSE160987_prep_rm <- GSE160987_20201015_Dubois_counts_rm[, -1]
rownames(GSE160987_prep_rm) <- GSE160987_20201015_Dubois_counts_rm[,1]

GSE160987_matrix <- as.matrix(GSE160987_prep_rm)

# Select only 28 days samples:
GSE160987_matrix_28days <- GSE160987_matrix[,1:30]
  
cols_GSSE160987 <- colnames(GSE160987_prep_rm)


col_data_GSE160987 <- read.delim("col_data_RNAseq_PPARd_paper.txt")

col_data_GSE160987 <- col_data_GSE160987 %>%
  mutate(
    sample_name = str_replace_all(sample_name, "[^[:alnum:]]", "."),
    run = str_replace_all(run, "[^[:alnum:]]", "."),
    treatment = str_replace_all(treatment, "[^[:alnum:]]", "_")
  )

col_data_GSE160987$treatment <- factor(col_data_GSE160987$treatment) 
col_data_GSE160987$run <- factor(col_data_GSE160987$run)

# Remove the first column and use its values as row names
col_data_GSE160987_rownames <- col_data_GSE160987 %>%
  select(-1)  # Remove the first column

# Set row names using values from the first column
rownames(col_data_GSE160987_rownames) <- col_data_GSE160987$sample_name

# Convert tibble to data frame to set row names
col_data_GSE160987_rownames <- as.data.frame(col_data_GSE160987_rownames)

rownames(col_data_GSE160987)
head(GSE160987_matrix_28days)

# Check if all rownames of the meta data are in the matrix
all(rownames(col_data_GSE160987_rownames) %in% colnames(GSE160987_matrix_28days))

all(rownames(col_data_GSE160987_rownames) == colnames(GSE160987_matrix_28days))

# Replace NA values with zero (adjust as needed)
GSE160987_matrix_28days[is.na(GSE160987_matrix_28days)] <- 0

GSE160987_dds <- DESeqDataSetFromMatrix(countData = GSE160987_matrix_28days,
                                        colData = col_data_GSE160987_rownames,
                                        design = ~ treatment)
GSE160987_dds$treatment <- relevel(GSE160987_dds$treatment, ref = "Control")

dds <- DESeq(GSE160987_dds)


res2 <- as.data.frame(results(DESeq(GSE160987_dds)))
res3 <- res2[order(res2$pvalue, res2$padj),]



res <- results(dds)
res
test <- res@listData

res05 <- results(dds, alpha=0.05)

pcaPlot(count = vst)

comparisons <- resultsNames(dds)

resOrdered <- res[order(res$pvalue),]
summary(resOrdered)
summary(res05)
sum(res$padj < 0.1, na.rm=TRUE)

# Plotting
plotMA(res05, ylim = c(-2, 2))

dds_data <- plotCounts(dds, gene=which.min(res$padj), intgroup = "treatment", returnData = TRUE)


library("pheatmap")
ntd <- normTransform(dds)
sampleDists <- dist(t(assay(ntd)))

select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[,c("treatment")])
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)

plotPCA(ntd, intgroup = "treatment")

sampleDist_matrix <- as.matrix(sampleDists)
pheatmap(sampleDist_matrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists)

vsd <- vst(dds, blind = FALSE)
head(assay(vsd), 3)

expr_norm_df <- as.data.frame(assay(vsd))
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
              "DECR2")

expr_norm_df_genes <- expr_norm_df[row.names(expr_norm_df) %in% genes_oi,]


pheatmap(
  expr_norm_df_genes
)



fold_change <- res$log2FoldChange[genes_of_interest]
print(fold_change)
# Create data frame for plotting
plot_data <- data.frame(
  treatment = colData(dds)[, treatment],
  gene = rep(valid_genes, each = ncol(assay(dds))),
  fold_change = exp(fold_change)  # Transform log2FoldChange to fold change
)

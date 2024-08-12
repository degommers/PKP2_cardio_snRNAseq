# Packages
library(tidyverse)
library(pheatmap)
library(tsne)
library(reshape2)
library(scales)
library(DESeq2)
library(gridExtra)

# Read and preprocess gene expression data
GSE160987_20201015_Dubois_counts <- read.csv("GSE160987_20201015_Dubois_counts.csv")

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
col_data_GSE160987 <- read.delim("col_data_RNAseq_PPARd_paper.txt") %>%
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
                       "DECR2")

matrix_LFC_genes <- subset(matrix_LFC, rownames(matrix_LFC) %in% genes_of_interest)

# Changing the names of the matrix to proper names:
colnames(matrix_LFC_genes) <- c("LCFA", "LCFA+PPARa (Ag)", "LCFA+PPARd (Ag)", 
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

pheatmap(matrix_LFC_genes, scale = "none", show_rownames = TRUE, 
         border_color = NA, fontsize = 12,
         clustering_distance_cols = "euclidean", 
         treeheight_row = 0, 
         cutree_cols = 5, 
         cutree_rows = 5, 
         color = colorRampPalette(c("blue", "white", "red"))(10), 
         annotation_col = annotation_col,
         annotation_colors = ann_colors)






# Get DESeq2 results and sort by p-value and adjusted p-value
res_GSE160987 <- results(dds, contrast = c("treatment", "Control"))
res_GSE160987 <- lfcShrink(dds, coef=9, res=res_GSE160987)

res_sorted <- res_GSE160987[order(res_GSE160987$pvalue, res_GSE160987$padj), ]



# important functions


library(pheatmap)  # Ensure pheatmap library is loaded if not already

plotHeatmap_dg <- function(vst_data = vst, res_data = res_GSE160987, condition = "treatment", 
                           fileName = NULL, genes_of_interest = NULL, 
                           top_genes = 500, mycol = 4, myrow = 4, 
                           method = "correlation") {
  if (is.null(genes_of_interest)) {
    # Use top 'top_genes' significant genes from 'res_data'
    if (top_genes > nrow(res_data)) {
      stop("top_genes exceeds the number of genes in res_data.")
    }
    
    topGenes <- unique(rownames(res_data)[order(abs(res_data$pvalue))][1:top_genes])
  } else {
    # Use specified 'genes_of_interest' if provided
    genes_of_interest <- intersect(genes_of_interest, rownames(res_data))  # Filter to genes present in res_data
    if (length(genes_of_interest) == 0) {
      stop("No valid genes_of_interest found in res_data.")
    }
    topGenes <- genes_of_interest
  }
  
  mat <- assay(vst_data)[topGenes, ]
  mat <- (mat - rowMin(mat)) / (rowMax(mat) - rowMin(mat))
  
  df <- as.data.frame(colData(vst_data)[, condition])
  rownames(df) <- colnames(assay(vst_data))
  colnames(df) <- "treatment" 
  
  
  print(df)
  
  newColors <- colorRampPalette(grDevices::rainbow(length(unique(df$treatment))))
  mycolors <- newColors(length(unique(df$treatment)))
  names(mycolors) <- unique(df$treatment)
  mycolors <- list(treatment = mycolors)
  
  if (!is.null(fileName)) {
    tiff(paste0(fileName, "_heatmap_", ifelse(is.null(genes_of_interest), paste("top", top_genes), paste(genes_of_interest, collapse = "_")), ".tif"), width = 18, height = 10, units = 'in', res = 300, compression = 'lzw')
  }
  
  pheatmap(mat, annotation_col = df, scale = "none", show_rownames = FALSE, border_color = NA, fontsize = 12,
           clustering_distance_cols = method, treeheight_row = 0, 
           cutree_cols = mycol, cutree_rows = myrow, 
           color = colorRampPalette(c("blue", "white", "red"))(10), annotation_colors = mycolors)
  
  if (!is.null(fileName)) {
    dev.off()
  }
}

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

plotHeatmap_dg(mycol = 4, vst_data = vsd, res_data = res_GSE160987, 
               condition = "treatment", genes_of_interest = genes_of_interest)


### PCA plot function)
pcaPlot_dg <- function(count = vst, res_select = NULL, condition = "treatment", fileName = NULL){
  if (class(count) == "DESeqTransform"){
    if (!is.null(res_select))
      count = count[rownames(count) %in% rownames(res_select)]
  } else if (class(count) == "DESeqDataSet"){
    if (is.null(res_select))
      stop("res_select should be provided if the count is 'DESeqDataSet'")
    count = count[rownames(count) %in% rownames(res_select)]
  } else stop("Unknown class of count")
  mydata = assay(count)
  
  mds <- plotPCA(count, intgroup=condition)
  ### Not entirely clear what intgroup does... !!!
  # intgroup - interesting groups: a character vector of names in colData(x) to use for grouping
  
  # colors
  mdsCond <- unique(mds$data[,condition])
  newColors <- colorRampPalette(grDevices::rainbow(length(mdsCond)))
  mycolors <- newColors(length(mdsCond))
  names(mycolors) <- mdsCond
  
  # PCA plot
  g = mds + geom_point(size=4, alpha=0.8) + 
    geom_vline(xintercept = 0, linetype=2, colour="black") + 
    geom_hline(yintercept = 0, linetype=2, colour="black") +
    # xlab("PC1") + ylab("PC2") + 
    theme_bw() + 
    theme(text = element_text(size=25), axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank()) + scale_colour_manual(name="Condition", values = mycolors)
  
  # save the plot
  if (!is.null(fileName)){
    tiff(paste0(fileName,"_PCA_vsd.tif"), width=9, height=6, units='in', dpi=300, compression='lzw')
  }
  
  print(g)
  
  if (!is.null(fileName)) dev.off()
}

pcaPlot_dg(count = vsd)


# TSNE

tsnePlot_dg <- function(count = vst, res_select = NULL, condition = "treatment", fileName = NULL, ...){
  if (class(count) == "DESeqTransform"){
    if (!is.null(res_select))
      count = count[rownames(count) %in% rownames(res_select)]
  } else if (class(count) == "DESeqDataSet"){
    if (is.null(res_select))
      stop("res_select should be provided if the count is 'DESeqDataSet'")
    count = count[rownames(count) %in% rownames(res_select)]
  } else stop("Unknown class of count")
  mydata = assay(count)
  
  tsne_data <- tsne(dist(t(mydata)), ...)
  tsne_data <- as.data.frame(tsne_data)
  tsne_data$group <- colData(count)[,condition]
  
  mdsCond <- unique(tsne_data$group)
  newColors <- colorRampPalette(grDevices::rainbow(length(mdsCond)))
  mycolors <- newColors(length(mdsCond))
  names(mycolors) <- mdsCond
  
  g = ggplot(tsne_data, aes(V1, V2, color=group)) + 
    geom_point(size=4, alpha=0.8) + 
    geom_vline(xintercept = 0, linetype=2, colour="black") + 
    geom_hline(yintercept = 0, linetype=2, colour="black") +
    theme_bw() + 
    theme(text = element_text(size=25), axis.line = element_line(colour = "black"), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          panel.border = element_blank(), panel.background = element_blank()) + 
    scale_colour_manual(name="Treatment", values = mycolors)
  print(g)
  
  # save the plot
  if (!is.null(fileName))
    ggsave(paste0(fileName,"_tsne_plot.tif"), width=9, height=6, units='in', dpi=300, compression='lzw')
}

tsnePlot_dg(count = vsd)

### Expression plot

library(ggplot2)
library(reshape2)  # For melt function
library(gridExtra)  # For arranging plots

expressionPlot_multi <- function(res, dds, genes = NULL, condition = "treatment", fileName = NULL) {
  if (is.null(genes)) {
    stop("Please provide a non-empty list of genes or specify the number of top genes.")
  }
  
  if (is.numeric(genes)) {
    topGenes <- unique(rownames(res)[order(abs(res$pvalue))][1:genes])
  } else {
    topGenes <- genes
  }
  
  # Filter to genes present in dds
  topGenes <- intersect(topGenes, rownames(assay(dds)))
  if (length(topGenes) == 0) {
    stop("No valid genes found in the dataset.")
  }
  
  # Subset expression data for selected genes and treatment conditions
  temp <- assay(dds)[topGenes, ]
  temp <- cbind.data.frame(t(temp), as.character(colData(dds)[, condition]))
  colnames(temp) <- c(head(colnames(temp), -1), "treatment")
  
  # Reshape data for plotting with ggplot2
  x <- melt(temp, id.vars = "treatment")
  
  # Create individual boxplot for each gene
  plots <- lapply(unique(x$variable), function(gene) {
    ggplot(data = subset(x, variable == gene), aes(x = treatment, y = value)) +
      geom_boxplot() +
      labs(title = gene) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  })
  
  # Arrange and display plots
  grid.arrange(grobs = plots, ncol = 2)  # Adjust ncol based on the number of plots you want per row
  
  # Save the plot
  if (!is.null(fileName)) {
    ggsave(paste0(fileName, "_expression.tif"), plot = grid.arrange(grobs = plots, ncol = 2),
           width = 9, height = 6, units = 'in', res = 300, compression = 'lzw')
  }
}

# Example usage:
selected_genes <- c("Gene1", "Gene2", "Gene3")  # Replace with your list of genes of interest
expressionPlot(res = res05, dds = GSE160987_dds, genes = genes_of_interest, 
               condition = "treatment", fileName = "expression_plots_17_FAO_genes")


expressionPlot_multi(res, GSE160987_dds, genes = c("ABCB11", "ACACB",
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
                                                   "DECR2"))

#########################
#    Plot Expression    #
#########################

expressionPlot <- function(res, dds, genes = 2, condition = "treatment", fileName = NULL) 
{
  if (is.numeric(genes)) {
    topGenes <- unique(rownames(res)[order(abs(res$pvalue))][1:genes])
  } else { 
    topGenes <- genes}
  
  temp <- assay(dds)[topGenes,]
  temp <- cbind.data.frame(t(temp),as.character(colData(dds)[,condition]))
  print(temp)
  colnames(temp) <- c(head(colnames(temp),-1),"treatment")
  x <- melt(temp, id.vars = "treatment")
  g <- ggplot(data=x, aes(x=treatment,y=value)) + geom_boxplot() + 
    facet_grid(vars(variable), scales = "free")
  
  print(g)
  
  # save the plot
  if (!is.null(fileName))
    ggsave(paste0(fileName,"_expression.tif"), width=9, height=6, units='in', res=300, compression='lzw')
}

expressionPlot(res, GSE160987_dds, genes = c("ABCB11", "ACACB",
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
                                                   "DECR2"))


# Fold change expression plot

expressionFoldChangeMatrix <- function(res, dds, genes = NULL, condition = "treatment") {
  if (is.null(genes) || length(genes) == 0) {
    stop("Please provide a non-empty list of genes.")
  }
  
  # Filter genes that are present in the DESeqDataSet
  valid_genes <- intersect(genes, rownames(assay(dds)))
  if (length(valid_genes) == 0) {
    stop("No valid genes found in the DESeqDataSet.")
  }
  
  # Get fold change values (log2FoldChange) from DESeq2 results for valid genes
  fold_change <- res$log2FoldChange[valid_genes]
  
  # Create data frame for plotting
  plot_data <- data.frame(
    treatment = colData(dds)[, condition],
    gene = rep(valid_genes, each = ncol(assay(dds))),
    fold_change = fold_change
  )
  
  print(plot_data)
  
  plot_data_sum <- plot_data %>%
       group_by(treatment) %>%
       summarize(
           meanlfc = mean(fold_change)
       ) %>%
       ungroup()
  
  print(plot_data_sum)
  
  # Create ggplot object
  g <- ggplot(data = plot_data_sum, aes(x = treatment, y = fold_change)) +
    geom_point() +
    facet_wrap(~ gene, scales = "free_y", ncol = 4) +  # Facet by genes with free y-axis scales
    theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels
  
  # Print the plot
  print(g)
  
  # Create matrix of fold change values per treatment type (cols) and genes of interest (rows)
  #fold_change_matrix <- reshape2::dcast(plot_data, gene ~ treatment, value.var = "fold_change")
  
  # Set row names of the matrix to genes of interest
  #rownames(fold_change_matrix) <- valid_genes
  
  #return(fold_change_matrix)
}

# Example usage of expressionPlot_FC function with specific gene list
selected_genes <- c("ABCB11", "ACACB", "SLC25A17", "DECR1", "CPT1A",
                    "EHHADH", "ACADS", "SESN2", "HADH", "SLC27A2",
                    "FABP1", "ETFDH", "MTLN", "ADIPOQ", "PPARA",
                    "ACAD11", "DECR2")

# Call expressionPlot_FC function with DESeq2 results and specific gene list
FC_matrix <- expressionPlot_FC(res = res_GSE160987, dds = dds, genes = selected_genes, 
                  condition = "treatment")


expressionFoldChangeMatrix <- function(res, dds, genes = NULL, condition = "treatment") {
  if (is.null(genes) || length(genes) == 0) {
    stop("Please provide a non-empty list of genes.")
  }
  
  # Filter genes that are present in the DESeqDataSet
  valid_genes <- intersect(genes, rownames(assay(dds)))
  if (length(valid_genes) == 0) {
    stop("No valid genes found in the DESeqDataSet.")
  }
  
  # Get fold change values (log2FoldChange) from DESeq2 results for valid genes
  fold_change <- res$log2FoldChange[valid_genes]
  
  
  # Create data frame for plotting
  plot_data <- data.frame(
    treatment = colData(dds)[, condition],
    gene = rep(valid_genes, each = ncol(assay(dds))),
    fold_change = fold_change
  )
  
  print(plot_data)
  
  # Group by treatment and calculate mean fold change
  plot_data_sum <- plot_data %>%
    group_by(treatment, gene) %>%
    summarize(mean_fold_change = mean(fold_change)) %>%
    ungroup()
  
  
  # Create ggplot object
  g <- ggplot(data = plot_data, aes(x = gene, y = fold_change)) +
    geom_boxplot() +
    facet_wrap(~ treatment, scales = "free_y", ncol = 4) +  # Facet by genes with free y-axis scales
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_y_continuous(labels = scientific)

  
  # Print the plot
  print(g)
  
  # Optionally, you can return the summarized data frame
  return(plot_data)
}

# Example usage of expressionFoldChangeMatrix function with specific gene list
selected_genes <- c("ABCB11", "ACACB", "SLC25A17", "DECR1", "CPT1A",
                    "EHHADH", "ACADS", "SESN2", "HADH", "SLC27A2",
                    "FABP1", "ETFDH", "MTLN", "ADIPOQ", "PPARA",
                    "ACAD11", "DECR2")

# Call expressionFoldChangeMatrix function with DESeq2 results and specific gene list
FC_data <- expressionFoldChangeMatrix(res = res_GSE160987, dds = dds, genes = selected_genes, 
                                      condition = "treatment")




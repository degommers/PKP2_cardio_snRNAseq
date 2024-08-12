#### Analysis of PKP2 and Control snRNAseq ####
# Author: Demi Gommers
# Date: May, 2024
# R Version: R4.4.0
# Rtools: Rtools44
# Rstudio: 

#### STEP 0: installing and loading library ####
if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
  BiocManager::install(version = "3.19") }

setRepositories(ind = 1:3, addURLs = c('https://satijalab.r-universe.dev', 'https://bnprks.r-universe.dev/'))

install.packages(c("BPCells", "presto", "glmGamPoi", "hdf5r"))
install.packages(c("shinydashboard", "DT", "shiny", "shinyjs"))
remotes::install_github("jokergoo/ComplexHeatmap")
nBiocManager::install("DESeq2")
BiocManager::install("Nebulosa")
BiocManager::install('multtest')

BiocManager::install(c("rhdf5", "DirichletMultinomial", "CNEr", "BSgenome", 
                       "beachmat", "SingleCellExperiment", "DelayedMatrixStats", 
                       "ensembldb", "TFBSTools", "glmGamPoi", "EnsDb.Hsapiens.v86", 
                       "BSgenome.Hsapiens.UCSC.hg38"), version = "3.19")
install.packages("XML")
BiocManager::install("KEGGgraph")
BiocManager::install("clusterProfiler")
BiocManager::install("pathview")
BiocManager::install("enrichplot")

# Install the remotes package
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
install.packages('Signac')
install.packages('metap')
install.packages("scales")
install.packages("scCustomize")
install.packages("fftw3")

library(remotes)
remotes::install_github("satijalab/azimuth", ref = "master") 
remotes::install_github("satijalab/seurat-wrappers", quiet = TRUE)
remotes::install_github("10XGenomics/loupeR")

##### Load Packages #####
library(dplyr)
library(ggh4x)
library(metap)
library(purrr)
library(Seurat)
library(SeuratData)
library(rhdf5)
library(Nebulosa)
library(stringr)
library(Matrix)
library(loupeR)
library(patchwork)
library(ggplot2)
library(tools)
library(hdf5r)
library(DESeq2)
library(edgeR)
library(ape)
library(multtest)
library(tidyr)
library(emmeans)
library(ggsignif)
library(ggpubr)
library(ComplexHeatmap)
library(ggsci)
library(ggrepel)
library(tibble)
library(clusterProfiler)
library(enrichplot)
library(pathview)
library(gridExtra)
library(scales)
library(scCustomize)
library(ggpmisc)
library(ggbreak)
library(powerjoin)

#### STEP 1: Load files and other variables ####
# Set standard colors
colors_thesis <- c(
  "#e6194b", "#3cb44b", "#ffe119", "#0082c8", "#f58231", "#808080", "#911eb4", "#46f0f0", "#f032e6",
  "#800000", "#d2f53c", "#fabebe", "#008080", "#e6beff", "#aa6e28", "#fffac8", "#aaffc3",
  "#808000", "#ffd8b1", "#000080", "#ffffff", "#000000")

colors_control_PKP2 <- c("#008080", "#DC143C")

my_pal <- pal_d3("category20")(12)
print(my_pal)

pal_cm <- colorRampPalette(c("lightgrey", "#FF7F0EFF"))(8)
pal_fb <- colorRampPalette(c("lightgrey", "#1F77B4FF"))(10)
pal_ad <- colorRampPalette(c("lightgrey", "#7F7F7FFF"))(6)


# PPRE gene_list
PPRE_geneset <- c("ABCB11", "ACACB", "SLC25A17", "DECR1",  "CPT1A",  "EHHADH",  "ACADS",  "SESN2",
                  "HADH", "HADHB", "HADHA", "SLC27A2", "FABP1", "ETFDH", 
                  "SMIM37", "ADIPOQ", "PPARA", "ACAD11", "DECR2", "FABP3", "PPARG", 
                  "PPARD", "LEP", "RXRA", "RXRB", "RXRG")

# Define directory where h5 files are located
dir_path <- "~/Documents/R_seurat_analysis (copy 1)/RV0_h5" 

load_files_h5 <- function(dir_path) {
  # List all .h5 files in the directory
  h5_files <- list.files(path = dir_path, pattern = "*.h5", full.names = TRUE)
  
  # Initialize a list to store Seurat objects
  seurat_objects <- list()
  
  # Loop through each file and create a Seurat object
  for (file in h5_files) {
    print(file)
    sample_name <- file_path_sans_ext(basename(file))
    sample <- str_extract(sample_name, pattern = "H[0-9]{2}")
    print(sample)
    counts <- Read10X_h5(file, use.names = TRUE, unique.features = TRUE)
    seurat_obj <- CreateSeuratObject(counts = counts, project = sample, assay = "RNA", min.cells = 3, min.features = 200)
    seurat_objects[[sample]] <- seurat_obj
  }
  
  return(seurat_objects)
}

# Execute load_files_h5 function
seurat_objects <- load_files_h5(dir_path)

#### STEP 2: Include mitochondrial metadata and check quality ####
# Create function to annotate mitochondrial genes and create violin plots
mito_anno_vln_plots <- function(seurat_objects) {
  # Create empty list to store objects and plots
  vln_plots <- list()
  
  # Loop through individual samples in the list of seurat objects
  for (i in names(seurat_objects)) {
    
    # Extract mitochondrial genes from seurat objects
    mito.gene <- grep(pattern = "^MT-", x = rownames(seurat_objects[[i]]), value = TRUE)
    name <- names(seurat_objects[i])
    print(name)
    
    # Calculate the percentage of mitochondrial genes and add to metadata of seurat object
    percent.mito <- colSums(GetAssayData(seurat_objects[[i]], layer = "counts")[mito.gene, ]) / colSums(GetAssayData(seurat_objects[[i]], layer = "counts"))
    seurat_objects[[i]] <- AddMetaData(seurat_objects[[i]], percent.mito, col.name = "percent.mito")
    
    # Create Violin plots for features nFeature_RNA, nCount_RNA, percent.mito
    vln_plots[[i]] <- VlnPlot(object = seurat_objects[[i]], features = c("nFeature_RNA", "nCount_RNA", "percent.mito"))
    
    pdf(paste0("VLNplot_", name, ".pdf"))
    print(vln_plots[[i]])
    dev.off()
  }
  
  return(seurat_objects)
}

# Execute mito_anno_vln_plots function on seurat objects
seurat_objects <- mito_anno_vln_plots(seurat_objects) 

#### STEP 3: Merging Seurat objects and include conditions and sample names ####
# Create function to subset and merge seurat objects
subset_and_merge_seurat_objects <- function(seurat_objects, list_samples, sample_conditions) {
  
  # Create empty lists to store seurat objects and violin plots
  subset_seurat_objects <- list()
  subset_vln_plots <- list()
  
  for (i in names(seurat_objects)) {
    print(i)
    print(paste("Original cell count:", ncol(seurat_objects[[i]])))
    
    subset_seurat_objects[[i]] <- subset(seurat_objects[[i]], 
                                         subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mito < 0.01 & nCount_RNA > 0 & nCount_RNA < 25000)
    
    subset_seurat_objects[[i]] <- NormalizeData(subset_seurat_objects[[i]], layer = "counts")
    subset_vln_plots[[i]] <- VlnPlot(object = subset_seurat_objects[[i]], features = c("nFeature_RNA", "nCount_RNA", "percent.mito"))
    
    # Create separate violin plot for subsetted seurat object
    pdf(paste0("plots_excl_h55/vlnplot/VLNplot_subset", i, ".pdf"))
    print(subset_vln_plots[[i]])
    dev.off()
  }
  
  # Merge subset_seurat_objects 
  subset_seurat_merge <- merge(subset_seurat_objects[[1]], y = subset_seurat_objects[2:6], 
                               add.cell.ids = list_samples, project = "PKP2")
  
  # Convert original idents and Conditions into factors
  subset_seurat_merge@meta.data$orig.ident <- as.factor(subset_seurat_merge@meta.data$orig.ident)
  subset_seurat_merge@meta.data$Condition <- as.factor(sample_conditions[subset_seurat_merge@meta.data$orig.ident])
  
  # Normalize and Join layers 
  subset_seurat_merge <- NormalizeData(subset_seurat_merge)
  subset_seurat_merge <- JoinLayers(subset_seurat_merge)
  
  return(subset_seurat_merge)
} 


# Define variables required for subsetting and merging:
# sample conditions
sample_conditions <- c(
  "H07" = "PKP2",
  "H08" = "PKP2",
  "H09" = "PKP2",
  "H11" = "PKP2",
  "H49" = "Control",
  "H53" = "Control"
)
# Samples
list_samples <- c("H07", "H08", "H09", "H11", "H49", "H53")

# Execute subset_and_merge_seurat_objects function
subset_seurat_merge_unintegrated <- subset_and_merge_seurat_objects(seurat_objects, list_samples, sample_conditions)

# Plot quality
Idents(subset_seurat_merge_unintegrated) <- "orig.ident"
pdf("plots_excl_h55/vlnplot/VLNPLot_merged_seurat.pdf", width = 10)
VlnPlot(subset_seurat_merge, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3,
        alpha = 0.2, pt.size = 0.2)
dev.off()

define_resolution_clusters <- function(subset_seurat_merge) {
  # This function is not required to run everytime, only one to check how much clusters you want
  all.genes <- rownames(subset_seurat_merge)
  subset_seurat_merge_res <- ScaleData(subset_seurat_merge, features = all.genes)
  subset_seurat_merge_res <- FindVariableFeatures(subset_seurat_merge_res, selection.method = "vst", nfeatures = 2000)
  subset_seurat_merge_res <- RunPCA(subset_seurat_merge_res, features = VariableFeatures(object = subset_seurat_merge_res))
  
  subset_seurat_merge_res <- FindNeighbors(subset_seurat_merge_res, dims = 1:16)
  subset_seurat_merge_res <- FindClusters(subset_seurat_merge_res, resolution = seq(0.2, 1.0, by = 0.2))
  
  subset_seurat_merge_res <- RunUMAP(subset_seurat_merge_res, dims = 1:16)
  
  ## Which resolution? --> 0.2
  
  umap0.4 <- DimPlot(subset_seurat_merge_res, reduction = "umap", group.by = "RNA_snn_res.0.4")
  umap0.2 <- DimPlot(subset_seurat_merge_res, reduction = "umap", group.by = "RNA_snn_res.0.2")
  umap0.6 <- DimPlot(subset_seurat_merge_res, reduction = "umap", group.by = "RNA_snn_res.0.6")
  umap0.8 <- DimPlot(subset_seurat_merge_res, reduction = "umap", group.by = "RNA_snn_res.0.8")
  umap1.0 <- DimPlot(subset_seurat_merge_res, reduction = "umap", group.by = "RNA_snn_res.1")
  
  pdf("plots_excl_h55/UMAP/UMAP_res_0.2-1.pdf")
  print(umap0.2)
  print(umap0.4)
  print(umap0.6)
  print(umap0.8)
  print(umap1.0)
  dev.off()
  
  return(subset_seurat_merge_res)
}

# Execute define_resolution_clusters --> interprete which resolution you want to use (I chose 0.2)
subset_seurat_merge_unintegrated <- define_resolution_clusters(subset_seurat_merge_unintegrated)

##### Analysis by dimension reduction #####
### Rescale, cluster and dimension reduction

transform_seurat <- function(subset_seurat_merge, res, reduction) {
  # Extract all genes from seurat object
  all.genes <- rownames(subset_seurat_merge)
  
  # Scale seurat object
  subset_seurat_merge <- ScaleData(subset_seurat_merge, features = all.genes)
  
  # Find variable 2000 features with vst method
  subset_seurat_merge <- FindVariableFeatures(subset_seurat_merge, selection.method = "vst", nfeatures = 2000)
  
  # Run PCA
  subset_seurat_merge <- RunPCA(subset_seurat_merge, features = VariableFeatures(object = subset_seurat_merge))
  
  # Cluster the cells based on earlier defined resoltion
  subset_seurat_merge <- FindNeighbors(subset_seurat_merge, dims = 1:20)
  subset_seurat_merge <- FindClusters(subset_seurat_merge, resolution = res) # Previous check on other resultions was done (see code at line XXX)
  subset_seurat_merge <- RunUMAP(subset_seurat_merge, dims = 1:20)
  
  return(subset_seurat_merge)
}

# Execute transform_seurat
subset_seurat_merge_unintegrated <- transform_seurat(subset_seurat_merge_unintegrated, res = 0.2)


### Check right clustering --> make sure that two conditions do not cluster together ###
# Plot PCA
pdf("plots_excl_h55/PCA_condition_samples.pdf")
DimPlot(subset_seurat_merge, reduction = "pca", group.by = c("Condition", "orig.ident"))
dev.off()

# Plot dimensions
pdf("plots_excl_h55/elbowplot.pdf")
ElbowPlot(subset_seurat_merge, ndims = 50)
dev.off()

pdf("plots_excl_h55/UMAP/UMAP_plot_final_res.pdf", width = 20)
DimPlot(subset_seurat_merge, reduction = "umap", split.by = c("orig.ident"), group.by = "Condition")
DimPlot(subset_seurat_merge, reduction = "umap", split.by = c("Condition"), group.by = "orig.ident")
dev.off()

pdf("plots_excl_h55/UMAP/UMAP_unintegrated_res_0.2.pdf", width = 15, height = 12)
DimPlot(subset_seurat_merge, reduction = "umap", alpha = 0.5, group.by = "Condition")
dev.off()

# Comment: Clear differences between the two conditions is seen, so I am trying to integrate with the cca method

#### STEP 4: Integration of layers ####
# Split the layers of conditions
split_seurat_object <- function(subset_seurat_merge) {
  # Split object into separate condition layers
  print("Splitting object in to layers per condition")
  subset_seurat_merge[["RNA"]] <- split(subset_seurat_merge[["RNA"]], f = subset_seurat_merge$Condition, layers = c("data", "counts"))
  
  all.genes <- rownames(subset_seurat_merge)
  print(head(all.genes))
  subset_seurat_merge <- NormalizeData(subset_seurat_merge)
  subset_seurat_merge <- FindVariableFeatures(subset_seurat_merge)
  subset_seurat_merge <- ScaleData(subset_seurat_merge, features = all.genes)
  subset_seurat_merge <- RunPCA(subset_seurat_merge)
  
  return(subset_seurat_merge)
}

subset_seurat_merge_split <- split_seurat_object(subset_seurat_merge_unintegrated)  


# Integrate layers with CCAIntegration method
subset_seurat_merge_split <- IntegrateLayers(object = subset_seurat_merge_split, 
                                       method = CCAIntegration, 
                                       assay = "RNA", 
                                       orig.reduction = "pca", 
                                       new.reduction = "integrated.cca")

# Rejoin layers
subset_seurat_merge <- JoinLayers(subset_seurat_merge_split)
FeaturePlot(subset_seurat_merge, features =  'HADH', slot =  "scale.data")

transform_integrated_seurat <- function(subset_seurat_merge) {
  # Recluster the cells
  subset_seurat_merge <- FindNeighbors(subset_seurat_merge, reduction = "integrated.cca", dims = 1:16)
  subset_seurat_merge <- FindClusters(subset_seurat_merge, resolution = 0.2, cluster.name = "integr.clusters")
  subset_seurat_merge <- RunUMAP(subset_seurat_merge, reduction = "integrated.cca", dims = 1:16)
  
  # Plotting UMAPs with different groups/splits
  umap_split_id <- DimPlot(subset_seurat_merge, reduction = "umap", split.by = c("orig.ident"), group.by = "Condition")
  umap_split_condition <- DimPlot(subset_seurat_merge, reduction = "umap", split.by = c("Condition"), group.by = "orig.ident")
  umap_condition <- DimPlot(subset_seurat_merge, reduction = "umap", alpha = 0.5, group.by = "Condition", label = TRUE)
  umap <- DimPlot(subset_seurat_merge, reduction = "umap", alpha = 0.5, label = TRUE)
  
  # Printing UMAPs
  pdf("plots_excl_h55/UMAP/UMAP_plot_integrated_split_orig_condition.pdf", width = 20)
  print(umap_split_id)
  print(umap_split_condition)
  dev.off()
  
  pdf("plots_excl_h55/UMAP/UMAP_integrated_condition.pdf", width = 15, height = 12)
  print(umap_condition)
  dev.off()
  
  pdf("plots_excl_h55/UMAP/UMAP_integrated.pdf", width = 15, height = 12)
  print(umap)
  dev.off()
  
  return(subset_seurat_merge)
}

subset_seurat_merge <- transform_integrated_seurat(subset_seurat_merge)

##### STEP 4b: Cluster annotation #####
genes_cluster_annotation_paper <- c("PPARG", "PLIN1", "ADIPOQ", "MGST1","ADIPOQ", "GPAM", "LEP", "ACACB", 
                                    "FER", "IL18R1", 
                                             "NTM", "SLC24A3", "ANK3", "XKR4", "CDH19", "NRXN3", "NRXN1", "ARHGAP15", 
                                             "PARP8", "ANKRD44", "SKAP1", "CD247", "PTPRC", "CD163", "FRMD4B", "RBM47", 
                                             "MRC1", "FMN1", "LDB2", "PTPRM", "EGFL7", "CDH5", "PKHD1L1", "VWF", "ST6GALNAC3", 
                                             "ABCA8", "ABCA6", "NEGR1", "GSN", "DCN", "VIM", "DCN", "FRMD3", "EPS8", "PLA2G5", "NR2F2-AS1", 
                                             "DLC1", "ACTA2", "CALD1", "RBM20", "MLIP",  "TTN", "TNNT2", "TNNI3", "RYR2", "MLIP", "RBM20",
                                    "TNNI3", "RYR2", "TTN", "PLN", "SLC8A1", "VWF"
                                    )

genes_cluster_annotation_paper_combined <- unique(genes_cluster_annotation_paper)
DimPlot(subset_seurat_merge, reduction = 'umap', cols = colors_thesis)

pdf("plots_excl_h55/dotplot_marker_genes_thesis.pdf", width = 15, height = 7)
Idents(subset_seurat_merge) <- "seurat_clusters"
DotPlot(subset_seurat_merge, features = genes_cluster_annotation_paper_combined,
        col.min = 0, col.max = 4, cols = c("lightgrey", "red")) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.2),
                                                                           legend.position = "top",
                                                                           legend.title = element_text(size = 10),
                                                                           legend.title.position = "top",
                                                                           legend.text = element_text(size = 8),
                                                                           legend.box = "vertical")
  
dev.off()

# DotPlot(subset(subset_seurat_merge, subset = seurat_annotations == "Unknown"), 
        # features = c("CD38", "CDH4", "DPP10", "GPAM", "LEP", "LSAMP", "CNTNAP2", "CAMTA1", "RBFOX3"))

####
# Function to annotate clusters of given seurat object 
cluster_annotation <- c(
  "FB", # 0
  "CM", # 1
  "Mural", # 2
  "EC", # 3
  "Monocytes", # 4
  "Unknown", # 5
  "CM", # 6
  "T-cells", # 7
  "AD", # 8
  "Mural", # 9
  "Myeloid", # 10
  "FB", # 11
  "SMC", # 12
  "NC", # 13
  "CM", # 14
  "Lymphoid" # 15
)
unique(cluster_annotation)
  

annotate_clusters <- function(seurat_object, cluster_annotation) {
  Idents(seurat_object) <- "seurat_clusters"
  
  names(cluster_annotation) <- levels(seurat_object)
  cluster_annotation_names <- levels(cluster_annotation)
  seurat_object <- RenameIdents(seurat_object, cluster_annotation)
  seurat_object <- AddMetaData(seurat_object, Idents(subset_seurat_merge), col.name = "seurat_annotations")
  
  return(seurat_object)
}

plot_annotated <- function(seurat_object, path = "plots_excl_h55/", marker_genes) {
  # Making the figures
  umap_labeled <- DimPlot(subset_seurat_merge, reduction = "umap", label = TRUE) + scale_color_manual(values = colors_thesis)
  umap_split_sample <- DimPlot(subset_seurat_merge, reduction = "umap", split.by = c("orig.ident"))
  
  Idents(seurat_object) <- "seurat_annotations"
  
  order_dotplot <- c("AD", "NC", "T-cells", "Monocytes", "EC", "FB", "Mural", "CM", "SMC", "Myeloid", "Lymphoid", "Unknown")
  dotplot_marker_genes <- DotPlot(subset_seurat_merge, features = marker_genes, group.by = "seurat_annotations", cols = "Spectral") + 
    scale_y_discrete(limits = order_dotplot) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.2),
          # legend.position = "top",
          legend.title = element_text(size = 10),
          legend.title.position = "top",
          legend.text = element_text(size = 8),
          legend.box = "vertical") 
  
  # Saving plots as PDF
  pdf(paste0(path, "UMAP/UMAP_cluster_annotation.pdf"), width = 10, height = 11)
  print(umap_labeled)
  dev.off()
  
  pdf(paste0(path, "UMAP/UMAP_cluster_annotation_split_sample.pdf"), width = 25, height = 6)
  print(umap_split_sample)
  dev.off()
  
  pdf(paste0(path, "dotplot/dotplot_marker_genes_paper_cluster_annotation.pdf"), width = 12, height = 6)
  print(dotplot_marker_genes)
  dev.off()
  
  print(paste("Plots are made and saved in", path))
}

subset_seurat_merge <- annotate_clusters(subset_seurat_merge, cluster_annotation)
plot_annotated(subset_seurat_merge, path = "plots_excl_h55/", genes_cluster_annotation_paper_combined)



# All markers
pdf("plots_excl_h55/feature_plot_cell_allmarkers.pdf", width = 26, height = 26)
FeaturePlot(subset_seurat_merge, features = genes_cluster_annotation_paper_combined)
dev.off()

# Find cluster markers
cluster.markers <- FindAllMarkers(subset_seurat_merge, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
head(cluster.markers)

cluster.markers.cl <- cluster.markers %>% 
  group_by(cluster) %>%
  filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup()




##### STEP 4b: Metrics on clusters #####
# Count cells per cluster
cell_counts <- as.data.frame(table(Idents(subset_seurat_merge), subset_seurat_merge$orig.ident, subset_seurat_merge$Condition)) %>%
  rename(cluster = "Var1", sample = "Var2", Condition = "Var3") %>%
  pivot_wider(
    names_from = Condition,
    values_from = Freq
  ) %>%
  pivot_longer(cols = c(PKP2, Control), values_to = "cell_count", names_to = "Condition") %>%
  filter(cell_count > 0) %>%
  group_by(sample) %>%
  mutate(
    total_cells = sum(cell_count),
    ratio_celltype = cell_count/total_cells
  )

cell_counts_raw <- as.data.frame(table(Idents(subset_seurat_merge), subset_seurat_merge$orig.ident, subset_seurat_merge$Condition))

# Statistics
### Function to plot the cell counts per condition faceted by cluster type
boxplot_cell_counts <- function(cell_counts) {
  my_pal <- pal_d3("category20")(12)
  strip_color <- strip_themed(background_x = elem_list_rect(fill = my_pal))
  
  ggplot(cell_counts, aes(Condition, ratio_celltype, fill = Condition))+
    geom_boxplot(width = 0.4, outlier.size = 1, outlier.color = "black") + 
    geom_point(position = "jitter", alpha = 0.3, size = 1) +
    facet_grid(~cluster, scales = "free") +
    scale_fill_manual(values = c("#008080", "#DC143C")) +
    stat_compare_means(label = "p.format", size = 3, method = "wilcox.test", vjust = 3) +
    theme_pubr() +
    labs(
      y = "fraction of celltype per total cells",
      x = ""
    )  +
    ylim(0, 1) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid.major = element_line(color = "#EBEBEBEB"),
      panel.grid.minor = element_line(color = "#EBEBEBEB"),
      panel.border = element_rect(color = "grey", fill = FALSE)
      # strip.background = element_rect(fill = "white", color = "darkgrey")
    )
}

boxplot_cell_counts(cell_counts)

bar_plot_count <- function(cell_counts, fill) {
  base_plot <- ggplot(cell_counts, aes(x = sample, y = ratio_celltype)) +
    geom_col(width = 0.5, color = 'black', position = "fill", show.legend = FALSE) +
    scale_fill_d3(palette = "category20") + 
    scale_y_continuous(labels = scales::percent) +
    theme_pubr() + 
    coord_flip() +
    labs(y = "Abundance", x = "", fill = "Cell type") +
    facet_grid(Condition~., scales = "free_y") +
    theme(
      strip.text = element_blank(),
      strip.background = element_rect(fill = "white", color = "white")
    )
  
  if (fill == "subcluster") {
    base_plot <- base_plot + aes(fill = subcluster)
  } else {
    base_plot <- base_plot + aes(fill = cluster)
  }
  
  return(base_plot)
}

bar_plot_count(cell_counts, fill = "seurat_annotations")



cellcounts_total <- ggplot(data = cell_counts %>% filter(cluster == "FB"), aes(y = sample, x = cluster, label = total_cells)) +
  geom_text(size = 3) +
  theme_pubr() +
  facet_grid(Condition~., scales = "free") + 
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid.major = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank()
  )

cellcounts_total
##### STEP 4c: PPRE responsive genes #####
# Feature plot per PPRE responsive gene:
pdf("plots_excl_h55/feature_plot_PPRE_responsive_genes.pdf", width = 15, height = 100)
FeaturePlot(subset_seurat_merge, features = PPRE_geneset, split.by = "Condition")
dev.off()

pdf("plots_excl_h55/feature_plot_PPRE_responsive_genes_density.pdf", width = 30, height = 30)
plot_density(subset_seurat_merge, features = PPRE_geneset, reduction = "umap", )
dev.off()

##### STEP 5a: Aggregate Expression Clusters #####
## Aggregate expression for all to create pseudobulk
pseudo_seurat_all <- AggregateExpression(subset_seurat_merge, assays = "RNA",
                                         return.seurat = TRUE, group.by = c("Condition", "orig.ident", "seurat_annotations"))
# Create right metadata (celltype.condition)
Idents(pseudo_seurat_all) <- "Condition"
# create a new column to annotate sample-condition-celltype in the single-cell dataset
pseudo_seurat_all$celltype.condition <- paste(pseudo_seurat_all$seurat_annotations, pseudo_seurat_all$Condition, sep = "_")


pseudo_counts <- table(Idents(pseudo_seurat_all))


### DE analysis ###
# Execute DE analysis of between PKP2 and control per cluster in sc and pseudobulk
clusters <- Idents(subset_seurat_merge)
subset_seurat_merge$seurat_annotations <- clusters
subset_seurat_merge$celltype.condition <- paste(subset_seurat_merge$seurat_annotations, subset_seurat_merge$Condition, sep = "_")

# create a new column to annotate sample-condition-celltype in the single-cell dataset
subset_seurat_merge$orig.ident.condition <- paste0(subset_seurat_merge$Condition, "-", subset_seurat_merge$orig.ident)

DE_pseudo <- FindMarkers(pseudo_seurat_all,
                         max.cells.per.ident = 1000,
                         ident.1 = "PKP2", ident.2 = "Control", 
                         test.use = "DESeq2",
                         verbose = FALSE)

top30_pseudo <- DE_pseudo %>%
  arrange(p_val_adj) %>%
  head(50) %>%
  rownames()

DE_pseudo <- DE_pseudo %>%
  rownames_to_column(var = "gene_symbol") %>%
  mutate(
    diff_exp = case_when(avg_log2FC > 0.6 & p_val < 0.05 ~ "UP",
                         avg_log2FC < -0.6 & p_val < 0.05 ~ "DOWN",
                         avg_log2FC < -0.5 & avg_log2FC >= -0.6 & p_val < 0.05 ~ "semi DOWN",
                         avg_log2FC > 0.5 & avg_log2FC <= 0.6 & p_val < 0.05 ~ "semi UP",
                         TRUE ~ "NO"),
    sig_label = if_else(gene_symbol %in% top30_pseudo, gene_symbol, NA_character_),
    PPRE_genes = if_else(gene_symbol %in% PPRE_geneset, gene_symbol, NA_character_)
  )
  
pdf("plots_excl_h55/volcanoplot/pseudobulk_DE_PPRE_labels.pdf", width = 12, height = 12)
ggplot(DE_pseudo, aes(x = avg_log2FC, y = -log10(p_val), colour = diff_exp, label = PPRE_genes)) +
  geom_point(size = 6, alpha = 0.5) + theme_minimal() +
  labs(
    y = expression("-log"[10]*"p-value"),
    x = expression("log"[2]*"FC"),
    color = "Differentially expressed"
  ) +
  theme(
    axis.line = element_line(color = "black"),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    plot.title = element_text(size = 18),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14)
  ) +
  scale_color_manual(values = c("#DC0000FF", "darkgrey","orange",  "orange", "#00A087FF"), 
                     labels = c("Downregulated", "Not significant", "Semi-downregulated", "Semi-upregulated", "Upregulated")) +
  geom_label_repel(size = 6, max.overlaps = 500, force = 10, label.padding = 0.15, max.time = 1.5, show.legend = FALSE, max.iter = 15000)

dev.off()

##### STEP 5b: PKP2 expression #####
pkp2_expression_plot <- VlnPlot(subset_seurat_merge, features = "PKP2", group.by = c("Condition", "orig.ident"), layer = "data", cols = colors_control_PKP2) + 
  stat_compare_means(method = "wilcox.test", label =  "p.format")

ggsave("plots_excl_h55/pkp2_expression_vln.png", pkp2_expression_plot, dpi = 600, width = 6, height = 8)

##### STEP 5c: Finding subcluster per cluster #####
subcluster_clusters <- function(seurat_object) {
  Idents(seurat_object) <- "seurat_annotations"
  
  for (i in unique(seurat_object$seurat_annotations)) {
    print(i)
    seurat_object <- FindSubCluster(seurat_object, cluster = i, graph.name = "RNA_snn", resolution = 0.5, subcluster.name = paste0("sub.", i))
  }
  
  ## Add subcluster column to metadata of seurat object
  # Combine the subcluster columns into a single column
  metadata <- seurat_object@meta.data 
  
  head(metadata)
  metadata2 <- metadata %>%
    rownames_to_column("cell_id") %>%
    # Gather all columns starting with "sub." into long format
    pivot_longer(cols = starts_with("sub."),
                 names_to = "cluster",
                 values_to = "subcluster",
                 values_drop_na = TRUE) %>%
    # Extract the subcluster identifier from the column names
    mutate(cluster = sub("sub\\.", "", cluster)) %>%
    filter(grepl("[0-9]", subcluster)) %>%
    select(cell_id, subcluster)
  
  # Join the new column with the original metadata
  metadata <- metadata %>%
    rownames_to_column(var = "cell_id") %>%
    left_join(metadata2, by = "cell_id") %>%
    column_to_rownames(var = "cell_id")
  
  # Add the new combined column back to the Seurat object
  seurat_object <- AddMetaData(object = seurat_object, metadata = metadata$subcluster, col.name = "subcluster")
  
  
  return(seurat_object)
}

subset_seurat_merge <- subcluster_clusters(subset_seurat_merge)

# Plot UMAPs
plot_UMAPs <- function(seurat_object, cluster_annotation) {
  UMAP_theme_split <- theme(axis.text = element_text(size = 5),
                            strip.text = element_text(size = 8),
                            legend.text = element_text(size = 7),
                            legend.title = element_text(size = 8))
  
  pdf("plots_excl_h55/UMAP/UMAP_subclusters_condition.pdf", width = 9, height = 3)
  
  for (i in unique(cluster_annotation)) {
    if(i == "T-cells") {
      next
    }
    group <- paste0("sub.", i)
    
    print(paste("Creating UMAP plot of", group))
    umap <- DimPlot(subset(seurat_object, subset = seurat_annotations == i), reduction = "umap", 
                    group.by = `group`, split.by = "Condition", shape.by = "Condition", pt.size = 0.75, label = TRUE,
                    label.size = 2, repel = TRUE, alpha = 0.5) + scale_color_npg() + 
      UMAP_theme_split
    print(umap)
  }
  dev.off()
}

plot_UMAPs(subset_seurat_merge, cluster_annotation)



# create a new column to annotate sample-condition-celltype in the single-cell dataset
subset_seurat_merge$subcluster.condition <- paste0(subset_seurat_merge$Condition, "-", subset_seurat_merge$subcluster)


pseudo_CM <- AggregateExpression(subset(subset_seurat_merge, subset = seurat_annotations == "CM"), normalization.method = "LogNormalize", return.seurat = TRUE, group.by = c("Condition", "subcluster"))
pseudo_CM$subcluster.condition <- paste0(pseudo_CM$Condition, "-", pseudo_CM$subcluster)
Idents(pseudo_CM) <- "subcluster.condition"

matrix_CM <- FetchData(pseudo_CM, layer = "scale.data", vars = c("PKP2", PPRE_geneset) )

pheatmap(matrix_CM)


Idents(pseudo_CM) <- "Condition"
DE_pseudo_CM <- FindMarkers(pseudo_CM,
                         max.cells.per.ident = 1000,
                         ident.1 = "PKP2", ident.2 = "Control", 
                         test.use = "DESeq2",
                         verbose = FALSE)

top30_pseudo_CM <- DE_pseudo_CM %>%
  arrange(p_val_adj) %>%
  head(50) %>%
  rownames()

DE_pseudo_CM <- DE_pseudo_CM %>%
  rownames_to_column(var = "gene_symbol") %>%
  mutate(
    diff_exp = case_when(avg_log2FC > 0.6 & p_val < 0.05 ~ "UP",
                         avg_log2FC < -0.6 & p_val < 0.05 ~ "DOWN",
                         avg_log2FC < -0.5 & avg_log2FC >= -0.6 & p_val < 0.05 ~ "semi DOWN",
                         avg_log2FC > 0.5 & avg_log2FC <= 0.6 & p_val < 0.05 ~ "semi UP",
                         TRUE ~ "NO"),
    sig_label = if_else(gene_symbol %in% top30_pseudo, gene_symbol, NA_character_),
    PPRE_genes = if_else(gene_symbol %in% PPRE_geneset, gene_symbol, NA_character_)
  )

pdf("plots_excl_h55/volcanoplot/pseudobulk_DE_PPRE_labels.pdf", width = 12, height = 12)
ggplot(DE_pseudo_CM, aes(x = avg_log2FC, y = -log10(p_val), colour = diff_exp, label = PPRE_genes)) +
  geom_point(size = 6, alpha = 0.5) + theme_minimal() +
  labs(
    y = expression("-log"[10]*"p-value"),
    x = expression("log"[2]*"FC"),
    color = "Differentially expressed"
  ) +
  theme(
    axis.line = element_line(color = "black"),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    plot.title = element_text(size = 18),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14)
  ) +
  scale_color_manual(values = c("#DC0000FF", "darkgrey","orange",  "orange", "#00A087FF"), 
                     labels = c("Downregulated", "Not significant", "Semi-downregulated", "Semi-upregulated", "Upregulated")) +
  geom_label_repel(size = 5, max.overlaps = 500, force = 10, label.padding = 0.15, max.time = 1.5, show.legend = FALSE, max.iter = 15000)

dev.off()




##### STEP 5d: subcluster analysis #####
seurat_CM <- subset(subset_seurat_merge, subset = seurat_annotations == "CM")
Idents(seurat_CM) <- "subcluster"

markers_CM_0_CM_4 <- FindMarkers(seurat_CM,
                                 ident.1 = "CM_0",
                                 ident.2 = "CM_4",
                                 test.use = "wilcox"
                                 )

top30_markers_CM_0_CM_4 <- markers_CM_0_CM_4 %>%
  arrange(p_val_adj) %>%
  head(50) %>%
  rownames()

markers_CM_0_CM_4 <- markers_CM_0_CM_4 %>%
  rownames_to_column(var = "gene_symbol") %>%
  mutate(
    diff_exp = case_when(avg_log2FC > 0.6 & p_val < 0.05 ~ "UP",
                         avg_log2FC < -0.6 & p_val < 0.05 ~ "DOWN",
                         avg_log2FC < -0.5 & avg_log2FC >= -0.6 & p_val < 0.05 ~ "semi DOWN",
                         avg_log2FC > 0.5 & avg_log2FC <= 0.6 & p_val < 0.05 ~ "semi UP",
                         TRUE ~ "NO"),
    sig_label = if_else(gene_symbol %in% top30_markers_CM_0_CM_4, gene_symbol, NA_character_),
    PPRE_genes = if_else(gene_symbol %in% PPRE_geneset, gene_symbol, NA_character_)
  )


pdf("plots_excl_h55/volcanoplot/CM0vsCM4_DE_PPRE_labels.pdf", width = 12, height = 12)
ggplot(markers_CM_0_CM_4, aes(x = avg_log2FC, y = -log10(p_val_adj), colour = diff_exp, label = sig_label)) +
  geom_point(size = 6, alpha = 0.5) + theme_minimal() +
  labs(
    y = expression("-log"[10]*"p-value"),
    x = expression("log"[2]*"FC"),
    color = "Differentially expressed"
  ) +
  theme(
    axis.line = element_line(color = "black"),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    plot.title = element_text(size = 18),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14)
  ) +
  scale_color_manual(values = c("#DC0000FF", "darkgrey","orange",  "orange", "#00A087FF"), 
                     labels = c("Downregulated", "Not significant", "Semi-downregulated", "Semi-upregulated", "Upregulated")) +
  geom_label_repel(size = 3, max.overlaps = 500, force = 10, label.padding = 0.15, max.time = 1.5, show.legend = FALSE, max.iter = 15000)

dev.off()


# GSEA per subcluster
organism = "org.Hs.eg.db"
library(organism, character.only = TRUE)
gene_list_CM0vsCM4 <- markers_CM_0_CM_4$avg_log2FC
names(gene_list_CM0vsCM4) <- markers_CM_0_CM_4$gene
gene_list_CM0vsCM4 <- na.omit(gene_list_CM0vsCM4)
gene_list_CM0vsCM4 <- sort(gene_list_CM0vsCM4, decreasing = TRUE)
gse_CM0vsCM4 <- gseGO(geneList = gene_list_CM0vsCM4,
                  ont = "ALL",
                  keyType = "SYMBOL",
                  minGSSize = 3,
                  maxGSSize = 10000,
                  pvalueCutoff = 0.05,
                  verbose = TRUE,
                  eps = 0,
                  OrgDb = organism,
                  pAdjustMethod = "bonferroni")


require(DOSE)
pdf("plots_excl/dotplot_GSE_pathways_CM_0.pdf", height = 10)
dotplot(gse_CM0vsCM4, showCategory = 10, split = ".sign", title = "CM0vsCM4") + 
  facet_grid(.~.sign)
dev.off()

##### STEP 5e: PPRE Expression #####
plot_PKP2_vs_PPRE <- FeatureScatter(seurat_CM, feature1 = "PKP2", 
                                    feature2 = "HADHA", split.by = "subcluster", 
                                    group.by = "Condition",
                                    cols = c("#008080", "#DC143C"))

strip_color_clusters <- strip_themed(background_y = elem_list_rect(fill = my_pal))
strip_color_CM <- strip_themed(background_y = elem_list_rect(fill = pal_cm))
strip_color_FB <- strip_themed(background_y = elem_list_rect(fill = pal_fb))
strip_color_AD <- strip_themed(background_y = elem_list_rect(fill = pal_ad))


pdf("plots_excl_h55/Fig4_expression_Vln.pdf", width = 30, height = 40)
VlnPlot(subset_seurat_merge, features = PPRE_geneset, split.by = "Condition", group.by = "seurat_annotations")
dev.off()


# Plotting PKP2 expression per sample
expression_pkp2_PPRE <- FetchData(subset_seurat_merge, vars = c("PKP2", PPRE_geneset, "Condition", 
                                                                "subcluster", "seurat_annotations"), layer = "data") 

expression_pkp2_PPRE_scale <- FetchData(subset_seurat_merge, vars = c("PKP2", PPRE_geneset, "Condition", 
                                                                "subcluster", "seurat_annotations"), layer = "scale.data") 


expression_pkp2_PPRE_logN <- as.data.frame(expression_pkp2_PPRE) %>%
  rownames_to_column("cell_id") %>%
  pivot_longer(cols = 2:28, names_to = "gene", values_to = "expression") %>%
  mutate(
    sample = str_extract(cell_id, "H\\d{2}")
  ) 

pdf("plots_excl_h55/Fig3_correct_data.pdf", height = 50, width = 35)
Fig4 <- ggplot(expression_pkp2_PPRE_logN, aes(x = gene, y = expression, fill = Condition)) +
  geom_point(position = position_jitterdodge(dodge.width = 0.85, jitter.width = 0.70), aes(color = Condition), size = 1, alpha = 0.5) + 
  scale_color_manual(values =  c("#008080", "#DC143C")) +
  facet_grid2(seurat_annotations~., scales = "free", strip = strip_color_clusters) +
  theme_pubr() +
  labs(
    x = "PPRE responsive genes",
    y = "LogNormalized expression"
  ) + 
  stat_compare_means(method = "wilcox.test", hide.ns = TRUE, label = "p.signif" , size = 8, vjust = 1.5)


Fig4 + theme(
    axis.text.x = element_text(angle = 90, hjust = 1, face = "italic", size = 35),
    axis.text.y = element_text(size = 35),
    axis.title.y = element_text(size = 40),
    axis.title.x = element_text(size = 40),
    legend.text = element_text(size = 35),
    legend.title = element_blank(),
    strip.text.y = element_text(size = 35),
    panel.spacing = unit(1.5, "lines")
  ) +
  guides(color = guide_legend(override.aes = list(size = 16)))

dev.off()

pdf("plots_excl_h55/Fig3_CM_FB_Mu_EC_un.pdf", height = 25, width = 35)
ggplot(expression_pkp2_PPRE_scale %>% filter(seurat_annotations %in% c("CM", "FB", "Mural", "EC", "Unknown")), aes(x = gene, y = expression, fill = Condition)) +
  geom_point(position = position_jitterdodge(dodge.width = 0.85, jitter.width = 0.70), aes(color = Condition), size = 1, alpha = 0.5) + 
  scale_color_manual(values = c("darkgreen", "red2")) +
  facet_grid2(seurat_annotations~., scales = "free", strip = strip_color) +
  theme_pubr() +
  ylim(-1, 12) +
  labs(
    x = "PPRE responsive genes"
  ) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, face = "italic", size = 35),
    axis.text.y = element_text(size = 35),
    axis.title.y = element_text(size = 40),
    axis.title.x = element_text(size = 40),
    strip.text.y = element_text(size = 35)
  ) +
  stat_compare_means(method = "wilcox.test", label = "p.format", size = 5)
dev.off()



# Subcluster CMs
cell_counts_subcluster <- as.data.frame(table(Idents(subset_seurat_merge), subset_seurat_merge$orig.ident, subset_seurat_merge$Condition, subset_seurat_merge$subcluster)) %>%
  rename(cluster = "Var1", sample = "Var2", Condition = "Var3", subcluster = "Var4") %>%
  pivot_wider(
    names_from = Condition,
    values_from = Freq
  ) %>%
  pivot_longer(cols = c(PKP2, Control), values_to = "cell_count", names_to = "Condition") %>%
  filter(cell_count > 0) %>%
  group_by(sample) %>%
  mutate(
    total_cells = sum(cell_count),
    ratio_celltype = cell_count/total_cells
  )

pdf("plots_excl_h55/Fig4_CM_subcluster.pdf", height = 15, width = 15)
ggplot(expression_pkp2_PPRE_logN %>% filter(seurat_annotations == "CM"), aes(x = gene, y = expression, fill = Condition)) +
  geom_point(position = position_jitterdodge(dodge.width = 0.85, jitter.width = 0.70), aes(color = Condition), size = 1, alpha = 0.5) +
  scale_color_manual(values =  c("#008080", "#DC143C")) +
  facet_grid2(subcluster~., scales = "free", strip = strip_color_CM) +
  theme_pubr() +
  labs(
    x = "PPRE responsive genes",
    y = "LogNormalized Expression"
  ) +  
  ylim(0,5) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, face = "italic", size = 16),
    axis.text.y = element_text(size = 16),
    axis.title.y = element_text(size = 18),
    axis.title.x = element_text(size = 18),
    strip.text.y = element_text(size = 16),
    legend.position = "right",
    panel.spacing = unit(1.5, "lines")
  ) +
  stat_compare_means(method = "wilcox.test", hide.ns = TRUE, label = "p.signif" , size = 5, vjust = 1) +
  guides(color = guide_legend(override.aes = list(size = 5)))
dev.off()

pdf("plots_excl_h55/Fig4_FB_subcluster.pdf", height = 25, width = 35)
ggplot(expression_pkp2_PPRE_scale %>% filter(seurat_annotations == "FB"), aes(x = gene, y = expression, fill = Condition)) +
  geom_point(position = position_jitterdodge(dodge.width = 0.85, jitter.width = 0.70), aes(color = Condition), size = 1, alpha = 0.5) +
  scale_color_manual(values =  c("#008080", "#DC143C")) +
  facet_grid2(subcluster~., scales = "free", strip = strip_color_FB) +
  theme_pubr() +
  labs(
    x = "PPRE responsive genes",
    y = "LogNormalized Expression"
  ) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, face = "italic", size = 35),
    axis.text.y = element_text(size = 35),
    axis.title.y = element_text(size = 40),
    axis.title.x = element_text(size = 40),
    strip.text.y = element_text(size = 35)
  ) +
  stat_compare_means(method = "wilcox.test", hide.ns = TRUE, label = "p.signif" , size = 3, symnum.args = symnum.args)
dev.off()

pdf("plots_excl_h55/Fig4_AD_subcluster.pdf", height = 10, width = 20)
ggplot(expression_pkp2_PPRE_scale %>% filter(seurat_annotations == "AD"), aes(x = gene, y = expression, fill = Condition)) +
  geom_point(position = position_jitterdodge(dodge.width = 0.85, jitter.width = 0.70), aes(color = Condition), size = 1, alpha = 0.5) +
  scale_color_manual(values =  c("#008080", "#DC143C")) +
  facet_grid2(subcluster~., scales = "free", strip = strip_color_AD) +
  theme_pubr() +
  ylim(0, 6.5) +
  labs(
    x = "PPRE responsive genes",
    y = "LogNormalized Expression"
  ) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, face = "italic", size = 16),
    axis.text.y = element_text(size = 16),
    axis.title.y = element_text(size = 18),
    axis.title.x = element_text(size = 18),
    strip.text.y = element_text(size = 16),
    legend.position = "right"
  ) +
  stat_compare_means(method = "wilcox.test", hide.ns = TRUE, label = "p.signif" , size = 2, symnum.args = symnum.args) +
  guides(color = guide_legend(override.aes = list(size = 5)))
dev.off()

ggplot(expression_pkp2_PPRE_scale %>% filter(seurat_annotations == "AD"), aes(y = gene, x = expression, fill = Condition)) +
  geom_density_ridges(alpha = 0.3, aes(point_color = Condition)) + 
  scale_fill_manual(values =  c("#008080", "#DC143C")) +
  scale_color_manual(values =  c("#008080", "#DC143C"))+
  coord_flip() + 
  facet_grid2(subcluster~.)

# Wilcoxon test between gene expression of PKP2 and control per PPRE gene per cluster
statistical_analysis_gene_expression <- function(seurat_object, PPRE_geneset) {
  pkp2 <- subset(subset_seurat_merge, subset = Condition == "PKP2")
  control <- subset(subset_seurat_merge, subset = Condition == "Control")
  
  wilcox_testresults <- list()
  for (gene in PPRE_geneset) {
    print(gene)
    for (cluster in unique(cluster_annotation)) {
      print(cluster)
      # Subset cluster data
      pkp2_cluster <- subset(pkp2, subset = seurat_annotations == cluster)
      control_cluster <- subset(control, subset = seurat_annotations == cluster)
      
      # Fetch expression data
      gene_expr_pkp2 <- as.numeric(unlist(FetchData(pkp2_cluster, vars = gene)))
      gene_expr_control <- as.numeric(unlist(FetchData(control_cluster, vars = gene)))
      
      print(head(gene_expr_pkp2))
      print(head(gene_expr_control))
    
      result <- wilcox.test(gene_expr_pkp2, gene_expr_control)
      
      wilcox_testresults[[paste0(gene, cluster)]] <- result

    }
  }
  
  return(wilcox_testresults)
}

PPRE_expression_percluster <- statistical_analysis_gene_expression(subset_seurat_merge, PPRE_geneset)

expression_genes <- expression_pkp2_PPRE_scale %>%
  spread(expression, gene)

expression_PKP2 <- expression_pkp2_PPRE %>%
  filter(gene == "PKP2") %>%
  select(expression, cell_id, Condition, subcluster) %>%
  rename(PKP_expression = "expression")


extract_expression_gene <- function(expression_data) {
  expression_list <- list()
  unique_genes <- unique(expression_data$gene)
  
  final_expression <- expression_data %>%
    select(cell_id, Condition, subcluster, seurat_annotations) %>%
    distinct()
  
  for (i in unique(expression_data$gene)) {
    print(i)
    column_name <- paste0(i, "_expression")
    gene_data <- expression_data %>%
      filter(gene == i) %>%
      select(expression, cell_id, Condition, subcluster, seurat_annotations) %>%
      rename(!!i := expression)
    
    expression_list[[i]] <- gene_data
    
    final_expression <- full_join(final_expression, gene_data, by = c("cell_id", "Condition", "subcluster", "seurat_annotations"))
  }
  
  # This does not work, try to join all the expression columns of all datasets together in one file
  
  
  return(final_expression)
}




# Include all genes and run spearman correlation with non zero log normalized counts and get top 10 upper and lower bound per subcluster of CMs
final_expression <- extract_expression_gene(expression_pkp2_PPRE)
final_expression_scale <- extract_expression_gene(expression_pkp2_PPRE_scale)

final_expression_data_scale <- final_expression_scale %>%
  pivot_longer(cols = -c(cell_id, Condition, subcluster, PKP2, seurat_annotations), names_to = "gene", values_to = "expression") 

pdf("plots_excl_h55/expression_pkp2_vs_PPRE.pdf", width = 25, height = 23)
ggplot(final_expression_data_scale %>% filter(seurat_annotations == "CM"), aes(x = PKP2, y = expression, color = Condition)) +
  geom_point(alpha = 0.5, size = 0.5) +
  scale_color_manual(values = c("#008080", "#DC143C")) +
  facet_grid(gene ~ subcluster + Condition, scales = "free_y") +
  theme_minimal() +
  stat_cor()
dev.off()

pdf("plots_excl_h55/featureplot/feature_plot_cell_cardiomyocytes_coexp_pkp2_hadhab_decr1_acacb.pdf", width = 15, height = 30)
pkp2_hadhb <- FeaturePlot(seurat_CM, features = c("HADHB", "PKP2"), blend = TRUE, split.by = "Condition")
pkp2_hadha <- FeaturePlot(seurat_CM, features = c("HADHA", "PKP2"), blend = TRUE, split.by = "Condition")
pkp2_decr1 <- FeaturePlot(seurat_CM, features = c("DECR1", "PKP2"), blend = TRUE, split.by = "Condition")
pkp2_ACACB <- FeaturePlot(seurat_CM, features = c("ACACB", "PKP2"), blend = TRUE, split.by = "Condition")

pkp2_hadha / pkp2_hadhb / pkp2_decr1 / pkp2_ACACB
dev.off()
ggplot(final_expression_data_scale %>% filter(subcluster == "CM_3"), aes(x = PKP2, y = expression, color = Condition)) +
  geom_point(alpha = 0.5, size = 0.5) +
  scale_color_manual(values = c("#008080", "#DC143C")) +
  facet_grid(Condition ~ gene , scales = "free_y") +
  stat_cor(method = "spearman", size = 2) +
  theme_minimal()

final_expression_no_pkp2 <- final_expression %>% select(!PKP2_expression)

final_expression_hm <- final_expression_scale %>%
  column_to_rownames("cell_id") %>%
  select(!c(subcluster, Condition, seurat_annotations)) %>%
  as.matrix()

pdf("plots_excl_h55/heatmap_PPRE_genes_PKP2_CM_subclusters.pdf", width = 15, height = 12)
DoHeatmap(subset(subset_seurat_merge, subset = seurat_annotations == 'CM'),
          features = c("PKP2", PPRE_geneset), group.by = c("subcluster.condition"), draw.lines = TRUE,
          group.colors = pal_cm) + scale_fill_gradientn(colors = c("blue", "white", "red"))
dev.off()

# make a pseudobulk per subcluster and plot the expression of pkp2 against the expression of the ppre genes
pseudo_seurat_subcluster <- AggregateExpression(subset_seurat_merge, assays = "RNA",
                                                return.seurat = TRUE, group.by = c("Condition", "orig.ident", "subcluster"))
Idents(pseudo_seurat_subcluster) <- "subcluster.condition"
DGE_subclusters_all <- FindAllMarkers(pseudo_seurat_subcluster,
                                            max.cells.per.ident = 200,
                                            test.use = "DESeq2",
                                            verbose = FALSE)

ggplot(DE_CM_sub_pseudo %>% filter(gene %in% c(PPRE_geneset, "PKP2")), aes(x = gene, y = avg_log2FC)) +
  geom_point() +
  geom_text(aes(label = gene)) +
  facet_grid(~comparison) +
  theme_pubr() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1)
  )

volcano_pseudo_plot(DE_CM_sub_pseudo, PPRE_geneset)
  

FeaturePlot(subset(subset_seurat_merge, subset = seurat_annotations == 'CM'),
             features = c("PKP2", "HADHA"), blend = TRUE, split.by = "Condition")

ggplot(mean_expr_subCM, aes(Identity, Feature, fill = mean_exp)) +
  geom_dotplot()

DimPlot(subset_seurat_merge, reduction = "umap", group.by = c("sub.AD", "sub.CM", "sub.FB", "sub.SMC"), label = TRUE)



# Create color panel
x1 <- c("pink", "red")
x2 <- c("pink", "blue")
panel <- data.frame(x1, x2)

pdf("plots_excl_h55/dotplot/PPRE_responsive_genes_bysubcluster_split_condition.pdf", width = 8, height = 25)
DotPlot(subset_seurat_merge, split.by = "Condition", group.by = c("subcluster"),
        features = PPRE_geneset, cols = panel, dot.scale = 3) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.1),
        axis.text.y = element_text(size = 6),
        axis.title.x = element_blank())
dev.off()

percent_stats_sub <- Percent_Expressing(seurat_CM, split_by = "Condition",
                                    features = c(PPRE_geneset, "PKP2"))
pheatmap(percent_stats_sub,
         )
pdf("plots_excl_h55/heatmap/heatmap_percentage_expressed_CM_subclusters.pdf", width = 10, height = 10)


list_of_clusters <- c("CM_CM_0")
Heatmap(percent_stats_sub %>% select(contains("CM")), column_split = 3, row_split = 3)
dev.off()

pdf("plots_excl_h55/heatmap/heatmap_percentage_expressed_FB_subclusters.pdf", width = 10, height = 10)
Heatmap(percent_stats %>% select(contains("FB")))
dev.off()



percent_stats_longer <- percent_stats %>%
  rownames_to_column("gene") %>%
  pivot_longer(
    cols = 2:139,
    values_to = "percentage_expression",
    names_to = "condition_cluster_subcluster"
  ) %>%
  separate(condition_cluster_subcluster, 
           into = c("Condition", "cluster", "subcluster"), sep = "_") %>% 
  mutate(
    comparison = paste0(cluster, "_", subcluster) # it is not a comparison, but the columnname is then the same with the DE dataset
  )


DE_CM_subclusters <- function(seurat_object) {
  CM_sub_sc_DE <- list()
  for (i in 0:6) {
    name <- paste0("CM_", i)
    
    id.1 <- paste0("PKP2-", name)
    id.2 <- paste0("Control-", name)
    print(paste("Comparing", id.1, "with", id.2))
    CM_sub_sc_DE[[name]] <- FindMarkers(subset_seurat_merge, slot = "data", 
                                        min.pct = -Inf,
                                        logfc.threshold = -Inf,
                                        min.cells.feature = 1,
                                        min.cells.group = 1,
                                        ident.1 = id.1, ident.2 = id.2, 
                                        verbose = FALSE) %>% 
      mutate(comparison = name) %>%
      rownames_to_column("gene")
  }
  DE_CM_sub <- bind_rows(CM_sub_sc_DE)
  return(DE_CM_sub)
}

Idents(subset_seurat_merge) <- "subcluster.condition"
DE_CM_sub <- DE_CM_subclusters(subset_seurat_merge)

get_topx_genes <- function(seurat_object, number) {
  CM_seurat <- subset(seurat_object, subset = seurat_annotations == "CM")
  top30_list <- list()
  for (subcluster in unique(CM_seurat$subcluster)){
    topx <- DE_CM_sub %>%
      filter(comparison == subcluster, avg_log2FC > 1) %>%
      arrange(p_val_adj) %>%
      head(number) %>%
      select(gene, avg_log2FC)
    
    top30_list[[subcluster]] <- topx
    
  }
  
  return(top30_list)
}

get_bottomx_genes <- function(seurat_object, number) {
  top30_list <- list()
  for (subcluster in unique(seurat_object$subcluster)){
    topx <- DE_CM_sub %>%
      filter(comparison == subcluster, avg_log2FC < -1) %>%
      arrange(p_val_adj) %>%
      head(number) %>%
      select(gene, avg_log2FC)
    
    top30_list[[subcluster]] <- topx
    
  }
  
  return(top30_list)
}

list_CM_top30 <- get_topx_genes(subset_seurat_merge, 30)
list_CM_top10_up <- get_topx_genes(seurat_CM, 10)
list_CM_top10_down <- get_bottomx_genes(seurat_CM, 10)
list_CM_top10_unlist <- unlist(list_CM_top10_up)
list_CM_top10_down_unlist <- unlist(list_CM_top10_down)
DoHeatmap(seurat_CM, features = list_CM_top30, group.by = "subcluster.condition", slot = "data")
DoHeatmap(seurat_CM, features = list_CM_top10_down_unlist, group.by = "subcluster.condition")


# Subclusters of CM
CM_0_genes_top30 <- as.list(list_CM_top30$CM_0$gene)
CM_2_genes_top30 <- as.list(list_CM_top30$CM_2$gene)
CM_4_genes_top30 <- as.list(list_CM_top30$CM_4$gene)
CM_4_genes <- list_CM_top30$CM_4$gene
 
DE_CM_2_sub <- DE_CM_sub %>%
  filter(comparison == "CM_2") %>%
  mutate(
    diff_exp = case_when(avg_log2FC > 0.6 & p_val < 0.05 ~ "UP",
                         avg_log2FC < -0.6 & p_val < 0.05 ~ "DOWN",
                         avg_log2FC < -0.5 & avg_log2FC >= -0.6 & p_val < 0.05 ~ "semi DOWN",
                         avg_log2FC > 0.5 & avg_log2FC <= 0.6 & p_val < 0.05 ~ "semi UP",
                         TRUE ~ "NO"),
    sig_label = if_else(gene %in% CM_2_genes_top30, gene, NA_character_),
    PPRE_genes = if_else(gene %in% PPRE_geneset, gene, NA_character_)
  )

DE_CM_3_sub <- DE_CM_sub %>%
  filter(comparison == "CM_3") %>%
  mutate(
    diff_exp = case_when(avg_log2FC > 0.6 & p_val < 0.05 ~ "UP",
                         avg_log2FC < -0.6 & p_val < 0.05 ~ "DOWN",
                         avg_log2FC < -0.5 & avg_log2FC >= -0.6 & p_val < 0.05 ~ "semi DOWN",
                         avg_log2FC > 0.5 & avg_log2FC <= 0.6 & p_val < 0.05 ~ "semi UP",
                         TRUE ~ "NO"),
    sig_label = if_else(gene %in% CM_0_genes_top30, gene, NA_character_),
    PPRE_genes = if_else(gene %in% PPRE_geneset, gene, NA_character_)
  )

DE_CM_4_sub <- DE_CM_sub %>%
  filter(comparison == "CM_4") %>%
  mutate(
    diff_exp = case_when(avg_log2FC > 0.6 & p_val < 0.05 ~ "UP",
                         avg_log2FC < -0.6 & p_val < 0.05 ~ "DOWN",
                         avg_log2FC < -0.5 & avg_log2FC >= -0.6 & p_val < 0.05 ~ "semi DOWN",
                         avg_log2FC > 0.5 & avg_log2FC <= 0.6 & p_val < 0.05 ~ "semi UP",
                         TRUE ~ "NO"),
    sig_label = if_else(gene %in% CM_4_genes_top30, gene, NA_character_),
    PPRE_genes = if_else(gene %in% PPRE_geneset, gene, NA_character_)
  )

volcano_plot <- ggplot(DE_CM_2_sub, aes(x = avg_log2FC, y = -log10(p_val), colour = diff_exp, label = PPRE_genes)) +
  geom_point(size = 5, alpha = 0.5) + theme_minimal() +
  labs(
    y = expression("-log"[10]*"p-value"),
    x = expression("log"[2]*"FC"),
    color = "Differentially expressed"
  ) +
  scale_color_manual(values = c("#DC0000FF", "darkgrey","orange",  "orange", "#00A087FF"), 
                     labels = c("Downregulated", "Not significant", "Semi-downregulated", "Semi-upregulated", "Upregulated")) +
  geom_label_repel(size = 5, max.overlaps = 200, label.padding = 0.15, max.time = 1)

pdf("plots/volcano_plots_CM_subcluster_PPRE_low_FAO.pdf", width = 18, height = 5) # Check error message!!
ggplot(DE_CM_sub %>% filter(comparison %in% c("CM_0", "CM_1", "CM_2", "CM_5")), aes(x = avg_log2FC, y = -log10(p_val), colour = diff_exp, label = sig_label)) +
  geom_point(size = 2, alpha = 0.5) + theme_bw() +
  facet_grid(~comparison, scales = "free") +
  labs(
    y = expression("-log"[10]*"p-value"),
    x = expression("log"[2]*"FC"),
    color = "Differentially expressed"
  ) +
  scale_color_manual(values = c("#DC0000FF", "darkgrey","orange",  "orange", "#00A087FF"), 
                     labels = c("Downregulated", "Not significant", "Semi-downregulated", "Semi-upregulated", "Upregulated")) +
  geom_label_repel(size = 2, max.overlaps = 200, label.padding = 0.15, max.time = 1)
dev.off()

pdf("plots_excl_h55/volcanoplot/CM_3_volcanoplot_top_labeled.pdf", width = 8, height = 10)
volcano_plot
dev.off()

FeaturePlot(seurat_CM, features = c("HADHA", "PKP2"), blend = TRUE, blend.threshold = 0.1)
expr_cor <- FeatureScatter(seurat_CM, feature1 = "PKP2", "HADHA", group.by = "Condition", slot = "scale.data")

pdf("plots_excl_h55/test_expr_cor_PKP2_HADHA.pdf", width = 15)
expr_cor
dev.off()

VlnPlot(seurat_CM, group.by = c("Condition"), split.by = "subcluster", 
        features = "PDK4", layer = "data", cols = pal_cm,
        alpha = 0.5)

VlnPlot(seurat_CM, group.by = c("Condition"), split.by = "orig.ident", 
        features = "XIST", layer = "data", cols = pal_cm,
        alpha = 0.5)

vp <- volcano_pseudo_plot(DE_CM_sub, genes_to_label = PPRE_geneset) +
  labs(title = paste("DE", "in", "subcluster CM_0", "PKP2 vs Control")) +
  theme_minimal() +
  theme(axis.text = element_text(size = 14),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 18),
        axis.title = element_text(size = 16),
        axis.line.x = element_line(color = "black"),
        axis.line.y = element_line(color = "black"))
print(vp)  




# Count cells per cluster
Idents(subset_seurat_merge) <- "seurat_annotations"
cell_counts_subclusters <- as.data.frame(table(Idents(subset_seurat_merge), subset_seurat_merge$orig.ident, subset_seurat_merge$Condition, subset_seurat_merge$subcluster)) %>%
  rename(cluster = "Var1", sample = "Var2", Condition = "Var3", subcluster = "Var4") %>%
  pivot_wider(
    names_from = Condition,
    values_from = Freq
  ) %>%
  pivot_longer(cols = c(PKP2, Control), values_to = "cell_count", names_to = "Condition") %>%
  filter(cell_count > 0) %>%
  group_by(sample) %>%
  mutate(
    total_cells = sum(cell_count),
    ratio_celltype = cell_count/total_cells
  )                                    
                                      
                                      
##### STEP 5f: Top DE genes per subclusters ####
Idents(subset_seurat_merge) <- "subcluster.condition"
gettop10_logfc <- function(seurat_object, range, cluster) {
  list_top10_genes <- list()
  for (condition in c("PKP2", "Control")) {
    for (i in range) {
      
      print(i)
      name <- paste0(condition, "-", cluster, "_", i)
      print(name)
      if (name == "Control-FB_7") {
        print("Skipping this cluster as it contains to little cells")
        next
      }
      top10_genes <- FindMarkers(subset(subset_seurat_merge, subset = seurat_annotations == cluster),
                                 ident.1 = name,
                                 logfc.threshold = -Inf
      ) %>%
        rownames_to_column(var = "gene") %>%
        filter(avg_log2FC > 0) %>%
        arrange(p_val_adj, avg_log2FC) %>%
        head(10) %>%
        select(gene, avg_log2FC, p_val_adj) %>%
        mutate(
          subcluster = name,
          cluster = cluster
        ) 
      
      top10_genes_wide <- top10_genes %>%
        select(gene, subcluster, avg_log2FC) %>%
        pivot_wider(names_from = subcluster,
                    values_from = avg_log2FC)
      
      list_top10_genes[[name]] <- top10_genes_wide
    }
  }
  return(list_top10_genes)
}
# Run gettop10 function for AD, CM and FB
list_top10_genes_AD <- gettop10_logfc(subset_seurat_merge, (0:5), "AD")
list_top10_genes_CM <- gettop10_logfc(subset_seurat_merge, (0:7), "CM")
list_top10_genes_FB <- gettop10_logfc(subset_seurat_merge, (0:9), "FB")

# combine dataframes of list_top10_genes
list_top10_genes_AD <- list_top10_genes_AD[-4] # no top 10 genes for AD_3
top10_genes_AD <- power_full_join(list_top10_genes_AD, by = "gene") 
top10_genes_FB <- power_full_join(list_top10_genes_FB, by = "gene")
top10_genes_CM <- power_full_join(list_top10_genes_CM, by = "gene")

top10_genes_AD_FB_CM <- full_join(top10_genes_CM, top10_genes_FB, by = "gene") %>%
  full_join(top10_genes_AD, by = "gene") %>%
  replace(is.na(.), 0) 

genes <- top10_genes_AD_FB_CM$gene 


top10_genes_AD_FB_CM <- top10_genes_AD_FB_CM %>% select(-gene)

rownames(top10_genes_AD_FB_CM) <- genes
top10_genes_matrix <- as.matrix(top10_genes_AD_FB_CM)

# Filter FC > 5
top10_genes_matrix_5 <- top10_genes_matrix[ rowSums(top10_genes_matrix != 0) >= 2, ]

# Filter FC < 2
top10_genes_matrix_2 <- top10_genes_matrix[ rowSums(top10_genes_matrix !=0) == 2, ]

pheatmap(top10_genes_matrix_2)

RidgePlot(subset(subset_seurat_merge, subset = seurat_annotations == "CM"), features = "PKP2", group.by = "subcluster.condition", cols = alpha(colour = pal_cm2, alpha = 0.5))

# Create metadata
## Extract condition and cluster from column names
condition <- sapply(strsplit(colnames(top10_genes_matrix_5), "-|_"), `[`, 1)
cluster <- sapply(strsplit(colnames(top10_genes_matrix_5), "-|_"), `[`, 2)

## Extract condition and cluster from column names
condition_2 <- sapply(strsplit(colnames(top10_genes_matrix_2), "-|_"), `[`, 1)
cluster_2 <- sapply(strsplit(colnames(top10_genes_matrix_2), "-|_"), `[`, 2)

presence <- as.data.frame(rowSums(top10_genes_matrix_5 != 0)) %>%
  rename(sum = "rowSums(top10_genes_matrix_5 != 0)")

# Create annotation data frame
annotation <- data.frame(
  Condition = condition,
  Cluster = cluster
)
rownames(annotation) <- colnames(top10_genes_matrix_5)

annotation_2 <- data.frame(
  Condition = condition_2,
  Cluster = cluster_2
)
rownames(annotation_2) <- colnames(top10_genes_matrix_2)





# Add clustering
DoHeatmap(subset(seurat_CM, subset = subcluster == "CM_4"), features = CM_4_genes, group.by = c("Condition"))

PPRE_cells <- WhichCells(seurat_CM, expression = RXRA > 1 & PPARA > 1, slot = "data")
# MEF2A is a TF speciic for cardiac development, use to proof data is correct
# Compare (FC between MEF2A and PPARA)
PPARA_MEF2 <- WhichCells(seurat_CM, expression = MEF2A > 1 & PPARA > 1, slot = "data")


plot_density(seurat_CM, features = c("MEF2A", "PPARA"), reduction = "umap")
plot_density(seurat_CM, features = "PPARA", reduction = "umap")

# PKP2 expression vs all other genes --> check correlation --> top 10 positive and negative correlation
# PPRE geneset vs all other genes
  # Non zero expression


top30_DE_ps_CM <- DE_CM_sub_pseudo %>%
  arrange(p_val_adj) %>%
  head(50) %>%
  rownames()

DE_CM_sub_pseudo <- DE_CM_sub_pseudo %>%
  group_by(comparison) %>%
  mutate(
    fdr = p.adjust(p_val_adj, method = "fdr"),
    diff_exp = case_when(avg_log2FC > 0.6 & p_val < 0.05 ~ "UP",
                         avg_log2FC < -0.6 & p_val < 0.05 ~ "DOWN",
                         avg_log2FC < -0.5 & avg_log2FC >= -0.6 & p_val < 0.05 ~ "semi DOWN",
                         avg_log2FC > 0.5 & avg_log2FC <= 0.6 & p_val < 0.05 ~ "semi UP",
                         TRUE ~ "NO"),
    sig_label = if_else(gene %in% top30_DE_ps_CM, gene, NA_character_),
    PPRE_genes = if_else(gene %in% PPRE_geneset, gene, NA_character_)
  )

DE_CM_sub_PPRE <- na.omit(DE_CM_sub, PPRE_genes)
DE_CM_sub_PPRE_pseudo <- DE_CM_sub_pseudo %>%
  filter(!is.na(PPRE_genes)) %>%
  mutate(
    gene = factor(gene, levels = sort(unique(gene))),
    comparison = factor(comparison),
    not_signf = if_else(-log10(fdr) == 0, "not significant", "significant")
  )


DE_CM_sub_PPRE_percent_pseudo <- percent_stats_longer %>%
  filter(cluster == "CM") %>%
  right_join(DE_CM_sub_PPRE_pseudo, by = c("gene", "comparison"))


Percent_Expressing()



pdf("plots_excl_h55/UMAP/UMAP_subclusters_whole.pdf", width = 7, height = 6)
DimPlot(subset(subset_seurat_merge, subset = seurat_annotations == "CM"), reduction = "umap", 
        group.by = c("sub.CM"), shape.by = "Condition", pt.size = 2, label = TRUE,
        label.size = 2, repel = TRUE, alpha = 0.5) + scale_color_npg()  
  
DimPlot(subset(subset_seurat_merge, subset = seurat_annotations == "FB"), reduction = "umap", 
        group.by = c("sub.FB"), shape.by = "Condition", pt.size = 2, label = TRUE,
        label.size = 2, repel = TRUE, alpha = 0.5) + scale_color_npg() 
DimPlot(subset(subset_seurat_merge, subset = seurat_annotations == "AD"), reduction = "umap", 
        group.by = c("sub.AD"), shape.by = "Condition", pt.size = 2, label = TRUE,
        label.size = 2, repel = TRUE, alpha = 0.5) + scale_color_npg() 
DimPlot(subset(subset_seurat_merge, subset = seurat_annotations == "Monocytes"), reduction = "umap",
        group.by = c("sub.Monocytes"), shape.by = "Condition",
        label = TRUE, label.size = 2, repel = TRUE, alpha = 0.5) + scale_color_npg() 
dev.off()





# perfrom DE analysis across subclusters, plot dotplot with FC and FDR



plots_dot <- list()
# Loop through each cluster type
for (cluster_type in unique(cluster_annotation)) {
  print(cluster_type)
  # Create DotPlot for the current cluster type
  dotplot <- DotPlot(subset(subset_seurat_merge, subset = seurat_annotations == "FB"), 
                     split.by = "Condition", group.by = "subcluster", features = PPRE_geneset,
                     assay = "RNA",
                     dot.scale = 3) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.1),
          axis.text.y = element_text(size = 9))
  print(dotplot)
  # Add the plot to the list
  plot_list[[cluster_type]] <- dotplot
}




# PPARG - RXRA
pdf("plots_excl_h55/featureplot/feature_plot_cell_cardiomyocytes_coexp_pkp2_hadhab_decr1_acacb.pdf", width = 15, height = 30)
pkp2_hadhb <- FeaturePlot(seurat_CM, features = c("HADHB", "PKP2"), blend = TRUE, split.by = "Condition")
pkp2_hadha <- FeaturePlot(seurat_CM, features = c("HADHA", "PKP2"), blend = TRUE, split.by = "Condition")
pkp2_decr1 <- FeaturePlot(seurat_CM, features = c("DECR1", "PKP2"), blend = TRUE, split.by = "Condition")
pkp2_ACACB <- FeaturePlot(seurat_CM, features = c("ACACB", "PKP2"), blend = TRUE, split.by = "Condition")

pkp2_hadha / pkp2_hadhb / pkp2_decr1 / pkp2_ACACB
dev.off()








# Checking counts per cell per subcluster
metadata <- subset_seurat_merge@meta.data %>% 
  rownames_to_column("cell_id") %>% 
  select(cell_id, nCount_RNA, subcluster, Condition, seurat_annotations) %>%
  group_by(cell_id, Condition, subcluster, seurat_annotations) %>%
  summarize(
    counts = sum(nCount_RNA)
    
  )

cell_counts_subcluster <- metadata %>%
  filter(seurat_annotations ==  "CM") %>%
  mutate(
    cell_id = factor(cell_id)
  ) %>%
  group_by(Condition, subcluster) %>%
  summarize(
    cell_counts = n()
  )

ggplot(metadata %>% filter(seurat_annotations == "CM"), aes(x = Condition, y = log(counts))) +
  geom_point(position = "jitter", size = 0.5, color = "lightgrey", alpha = 0.3) +
  geom_boxplot(aes(fill = Condition), alpha = 0.5, outlier.size = 0.5, outlier.color = "darkgrey") +
  theme_bw() +
  scale_fill_npg() +
  facet_grid(~subcluster) +
  stat_compare_means(method = "wilcox.test", label = "p.format", size = 2) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.2)
    
  )
  
dev.off()

##### STEP 5g: Violin plots of specific genes (co expression) #####

vp_case1 <- function(gene_signature, file_name, test_sign){
  plot_case1 <- function(signature, y_max = NULL){
    VlnPlot(seurat_CM, features = signature,
            pt.size = 0.1, 
            split.by = "Condition",
            y.max = y_max # add the y-axis maximum value - otherwise p-value hidden
    ) + stat_compare_means(comparisons = test_sign, label = "p.signif")
  }
  plot_list <- list()
  y_max_list <- list()
  
  for (gene in gene_signature) {
    plot_list[[gene]] <- plot_case1(gene)
    y_max_list[[gene]] <- max(plot_list[[gene]]$data[[gene]]) # get the max no. for each gene
    plot_list[[gene]] <- plot_case1(gene, y_max = (y_max_list[[gene]] + 1) )
  }
  cowplot::plot_grid(plotlist = plot_list)
  file_name <- paste0(file_name, "_r.png")
  ggsave(file_name, width = 14, height = 8)
}

gene_sig <- c("HADHB", "HADH", "HADHA")
comparisons <- list(c("Control-CM_0", "PKP2-CM_0", "Control-CM_1", "PKP2-CM_1", "Control-CM_2", "PKP2-CM_2",
  "Control-CM_3", "PKP2-CM_3", "Control-CM_4", "PKP2-CM_4", "Control-CM_5", "PKP2-CM_5",
  "Control-CM_6", "PKP2-CM_6", "Control-CM_7", "PKP2-CM_7"))

vp_case1(gene_signature = gene_sig, file_name = "plots_excl_h55/gene_sig_2", test_sign = comparisons)


vlnplot_d <- VlnPlot(subset(subset_seurat_merge, subset = seurat_annotations == 'CM'), features = c("HADH", "HADHA", "HADHB", "ACACB"), 
        cols = c("#008080", "#DC143C", "grey"), split.by = "Condition", 
        layer = "data", split.plot = FALSE,
        ncol = 1,
        alpha = 0.5) + stat_compare_means(method = "wilcox.test", label = "p.signif")

vlnplot_d[[1]] <- vlnplot_d[[1]] + stat_compare_means(method = "wilcox.test", label = "p.signif") + RestoreLegend() + theme(legend.position = "top")
vlnplot_d[[2]] <- vlnplot_d[[2]] + stat_compare_means(method = "wilcox.test", label = "p.signif")
vlnplot_d[[3]] <- vlnplot_d[[3]] + stat_compare_means(method = "wilcox.test", label = "p.signif")

vlnplot_d
RidgePlot(subset(subset_seurat_merge, subset = seurat_annotations == "CM"), features = c("HADH", "HADHA", "HADHB", "ACACB"), group.by = "subcluster.condition", ncol = 2, cols = alpha(colour = pal_cm2, alpha = 0.5))

vlnplot_p <- VlnPlot(seurat_CM, features = c("PPARA", "PPARG", "PPARD"), 
                     cols = c("#008080", "#DC143C", "grey"), split.by = "Condition", 
                     layer = "data", split.plot = FALSE,
                     ncol = 1,
                     alpha = 0.5) + stat_compare_means(method = "wilcox.test", label = "p.signif")

vlnplot_p[[1]] <- vlnplot_p[[1]] + stat_compare_means(method = "wilcox.test", label = "p.signif") + RestoreLegend() + theme(legend.position = "top")
vlnplot_p[[2]] <- vlnplot_p[[2]] + stat_compare_means(method = "wilcox.test", label = "p.signif")
vlnplot_p[[3]] <- vlnplot_p[[3]] + stat_compare_means(method = "wilcox.test", label = "p.signif")

vlnplot_r <- VlnPlot(seurat_CM, features = c("RXRA", "RXRB", "RXRG"), 
                     cols = c("#008080", "#DC143C", "grey"), split.by = "Condition", 
                     layer = "data", split.plot = FALSE,
                     ncol = 1,
                     alpha = 0.5) + stat_compare_means(method = "wilcox.test", label = "p.signif")

vlnplot_r[[1]] <- vlnplot_r[[1]] + stat_compare_means(method = "wilcox.test", label = "p.signif") + RestoreLegend() + theme(legend.position = "top")
vlnplot_r[[2]] <- vlnplot_r[[2]] + stat_compare_means(method = "wilcox.test", label = "p.signif")
vlnplot_r[[3]] <- vlnplot_r[[3]] + stat_compare_means(method = "wilcox.test", label = "p.signif")

png("plots_excl_h55/expression_PPRE_vlnplots.png", width = 750, height = 900)
vlnplot_p | vlnplot_r
dev.off()

VlnPlot(seurat_CM, features = c("ACACB", "HADHA"))


blend_PKP2_HADHA <- FeaturePlot(seurat_CM, features = c("PKP2", "HADHA"), split.by = "Condition", blend = TRUE)
blend_PKP2_HADH <- FeaturePlot(seurat_CM, features = c("PKP2", "HADH"), split.by = "Condition", blend = TRUE)

# First plot
data_blend_PKP2_HADHA <- blend_PKP2_HADHA[[7]][["data"]] %>%
  mutate(
    co_expr = case_when(PKP2_HADHA %in% c(90, 9) ~ "co_expressed")
  ) %>%
  filter(co_expr == "co_expressed") %>%
  rownames_to_column(var = "cell_id") %>%
  mutate(
    Condition = "PKP2",
    sample = str_extract(cell_id, "H[0-9]{2}")
  ) %>%
  group_by(ident, sample, Condition, PKP2_HADHA) %>%
  summarize(
    n = n()
  )

data_blend_PKP2_HADH <- blend_PKP2_HADH[[7]][["data"]] %>%
  mutate(
    co_expr = case_when(PKP2_HADH %in% c(90, 9) ~ "co_expressed")
  ) %>%
  filter(co_expr == "co_expressed") %>%
  rownames_to_column(var = "cell_id") %>%
  mutate(
    Condition = "PKP2",
    sample = str_extract(cell_id, "H[0-9]{2}")
  ) %>%
  group_by(ident, sample, Condition) %>%
  summarize(
    n = n()
  )

# Second plot
data_blend_PKP2_HADHA_c <- blend_PKP2_HADHA[[3]][["data"]] %>%
  mutate(
    co_expr = case_when(PKP2_HADHA %in% c(90, 9) ~ "co_expressed")
  ) %>%
  filter(co_expr == "co_expressed") %>%
  rownames_to_column(var = "cell_id") %>%
  mutate(
    Condition = "Control",
    sample = str_extract(cell_id, "H[0-9]{2}")
  ) %>%
  group_by(ident, sample, Condition) %>%
  summarize(
    n = n()
  )

PKP2_HADHA_blend <- rbind(data_blend_PKP2_HADHA, data_blend_PKP2_HADHA_c)

HADHA_PKP2 <- ggplot(PKP2_HADHA_blend, aes(x = Condition, y = n, fill = Condition)) +
  geom_boxplot() +
  geom_point(position = "jitter") +
  facet_grid(~ident, scales = "free") +
  stat_compare_means(method = "wilcox.test") +
  theme_pubr() +
  scale_fill_manual(values = c("#008080", "#DC143C"))

# PPRE
Idents(subset_seurat_merge) <- "subcluster"
get_coexpression_count <- function(seurat_object, cluster, feature1, feature2, cell_counts_subclusters = cell_counts_subclusters) {
  # Create the combined expression feature name
  co_expression <- paste0(feature1, "_", feature2)
  print(co_expression)
  
  # Generate the blend plot
  blend_plot <- FeaturePlot(subset(seurat_object, subset = seurat_annotations == cluster), blend.threshold = 0.25,
                            features = c(feature1, feature2), split.by = "Condition", blend = TRUE)
  
  co_expression_plot <- FeatureScatter(subset(seurat_object, subset = seurat_annotations == cluster), feature1 = feature1, 
                                       feature2 = feature2, split.by = "subcluster", group.by = "Condition", cols = colors_control_PKP2)
  
  # # Helper function to process data for a given condition
  # process_blend_data <- function(data, condition) {
  #   data %>%
  #     rownames_to_column(var = "cell_id") %>%
  #     mutate(
  #       co_expr = case_when(
  #         .data[[co_expression]] == 99 ~ "co_expressed",
  #         .data[[co_expression]] %in% c(0, 9, 90) ~ "not_co_expressed",
  #         TRUE ~ NA_character_
  #       ),
  #       Condition = condition,
  #       sample = str_extract(cell_id, "H[0-9]{2}")
  #     ) %>%
  #     group_by(ident, sample, Condition, co_expr) %>%
  #     summarize(n = n(), .groups = 'drop')
  # }
  # 
  # # Process data for the first condition (PKP2)
  # print("getting data from co expression of condition 1")
  # data_blend <- process_blend_data(blend_plot[[7]][["data"]], "PKP2")
  # 
  # # Process data for the second condition (Control)
  # print("getting data from co expression of condition 2")
  # data_blend_control <- process_blend_data(blend_plot[[3]][["data"]], "Control")
  # 
  # # Combine data from both conditions
  # data_blend_all <- bind_rows(data_blend, data_blend_control)
  # 
  # # Add total cell counts to the combined data
  # print("adding cell count data to blend data")
  # cell_counts_subcluster_cluster <- cell_counts_subclusters %>%
  #   filter(cluster == cluster) %>%
  #   select(sample, cell_count, subcluster, Condition)
  # 
  # data_blend_all_count <- data_blend_all %>%
  #   rename(subcluster = ident) %>%
  #   left_join(cell_counts_subcluster_cluster, by = c("subcluster", "sample", "Condition")) %>%
  #   group_by(co_expr) %>%
  #   mutate(ratio_co_exp = n / cell_count,
  #          cluster = cluster)
  
  return(co_expression_plot)
}

## PPARA and RXRA
CM_coexpression_PPARA_RXRA <- get_coexpression_count(subset_seurat_merge, "CM", "RXRA", "PPARA", cell_counts_subclusters)
AD_coexpression_PPARA_RXRA <- get_coexpression_count(subset_seurat_merge, "AD", "RXRA", "PPARA", cell_counts_subclusters)
FB_coexpression_PPARA_RXRA <- get_coexpression_count(subset_seurat_merge, "FB", "RXRA", "PPARA", cell_counts_subclusters)

coexpression_PPARA_RXRA <- rbind(CM_coexpression_PPARA_RXRA, AD_coexpression_PPARA_RXRA, FB_coexpression_PPARA_RXRA) %>%
  mutate(dimer = "PPARA_RXRA")

## PPARG and RXRA
CM_coexpression_PPARG_RXRA <- get_coexpression_count(subset_seurat_merge, "CM", "RXRA", "PPARG", cell_counts_subclusters)
AD_coexpression_PPARG_RXRA <- get_coexpression_count(subset_seurat_merge, "AD", "RXRA", "PPARG", cell_counts_subclusters)
FB_coexpression_PPARG_RXRA <- get_coexpression_count(subset_seurat_merge, "FB", "RXRA", "PPARG", cell_counts_subclusters)

coexpression_PPARG_RXRA <- rbind(CM_coexpression_PPARG_RXRA, AD_coexpression_PPARG_RXRA, FB_coexpression_PPARG_RXRA) %>%
  mutate(dimer = "PPARG_RXRA")

## PPARD and RXRA
CM_coexpression_PPARD_RXRA <- get_coexpression_count(subset_seurat_merge, "CM", "RXRA", "PPARD", cell_counts_subclusters)
AD_coexpression_PPARD_RXRA <- get_coexpression_count(subset_seurat_merge, "AD", "RXRA", "PPARD", cell_counts_subclusters)
FB_coexpression_PPARD_RXRA <- get_coexpression_count(subset_seurat_merge, "FB", "RXRA", "PPARD", cell_counts_subclusters)

coexpression_PPARD_RXRA <- rbind(CM_coexpression_PPARD_RXRA, AD_coexpression_PPARD_RXRA, FB_coexpression_PPARD_RXRA) %>%
  mutate(dimer = "PPARD_RXRA")


## PPARA and RXRB
CM_coexpression_PPARA_RXRB <- get_coexpression_count(subset_seurat_merge, "CM", "RXRB", "PPARA", cell_counts_subclusters)
AD_coexpression_PPARA_RXRB <- get_coexpression_count(subset_seurat_merge, "AD", "RXRB", "PPARA", cell_counts_subclusters)
FB_coexpression_PPARA_RXRB <- get_coexpression_count(subset_seurat_merge, "FB", "RXRB", "PPARA", cell_counts_subclusters)

coexpression_PPARA_RXRB <- rbind(CM_coexpression_PPARA_RXRB, AD_coexpression_PPARA_RXRB, FB_coexpression_PPARA_RXRB) %>%
  mutate(dimer = "PPARA_RXRB")

## PPARG and RXRB
CM_coexpression_PPARG_RXRB <- get_coexpression_count(subset_seurat_merge, "CM", "RXRB", "PPARG", cell_counts_subclusters)
AD_coexpression_PPARG_RXRB <- get_coexpression_count(subset_seurat_merge, "AD", "RXRB", "PPARG", cell_counts_subclusters)
FB_coexpression_PPARG_RXRB <- get_coexpression_count(subset_seurat_merge, "FB", "RXRB", "PPARG", cell_counts_subclusters)

coexpression_PPARG_RXRB <- rbind(CM_coexpression_PPARG_RXRB, AD_coexpression_PPARG_RXRB, FB_coexpression_PPARG_RXRB) %>%
  mutate(dimer = "PPARG_RXRB")

## PPARD and RXRB
CM_coexpression_PPARD_RXRB <- get_coexpression_count(subset_seurat_merge, "CM", "RXRB", "PPARD", cell_counts_subclusters)
AD_coexpression_PPARD_RXRB <- get_coexpression_count(subset_seurat_merge, "AD", "RXRB", "PPARD", cell_counts_subclusters)
FB_coexpression_PPARD_RXRB <- get_coexpression_count(subset_seurat_merge, "FB", "RXRB", "PPARD", cell_counts_subclusters)

coexpression_PPARD_RXRB <- rbind(CM_coexpression_PPARD_RXRB, AD_coexpression_PPARD_RXRB, FB_coexpression_PPARD_RXRB) %>%
  mutate(dimer = "PPARD_RXRB")

## PPARA and RXRG
CM_coexpression_PPARA_RXRG <- get_coexpression_count(subset_seurat_merge, "CM", "RXRG", "PPARA", cell_counts_subclusters)
AD_coexpression_PPARA_RXRG <- get_coexpression_count(subset_seurat_merge, "AD", "RXRG", "PPARA", cell_counts_subclusters)
FB_coexpression_PPARA_RXRG <- get_coexpression_count(subset_seurat_merge, "FB", "RXRG", "PPARA", cell_counts_subclusters)

coexpression_PPARA_RXRG <- rbind(CM_coexpression_PPARA_RXRG, AD_coexpression_PPARA_RXRG, FB_coexpression_PPARA_RXRG) %>%
  mutate(dimer = "PPARA_RXRG")

## PPARG and RXRG
CM_coexpression_PPARG_RXRG <- get_coexpression_count(subset_seurat_merge, "CM", "RXRG", "PPARG", cell_counts_subclusters)
AD_coexpression_PPARG_RXRG <- get_coexpression_count(subset_seurat_merge, "AD", "RXRG", "PPARG", cell_counts_subclusters)
FB_coexpression_PPARG_RXRG <- get_coexpression_count(subset_seurat_merge, "FB", "RXRG", "PPARG", cell_counts_subclusters)

coexpression_PPARG_RXRG <- rbind(CM_coexpression_PPARG_RXRG, AD_coexpression_PPARG_RXRG, FB_coexpression_PPARG_RXRG) %>%
  mutate(dimer = "PPARG_RXRG")

## PPARD and RXRG
CM_coexpression_PPARD_RXRG <- get_coexpression_count(subset_seurat_merge, "CM", "RXRG", "PPARD", cell_counts_subclusters)
AD_coexpression_PPARD_RXRG <- get_coexpression_count(subset_seurat_merge, "AD", "RXRG", "PPARD", cell_counts_subclusters)
FB_coexpression_PPARD_RXRG <- get_coexpression_count(subset_seurat_merge, "FB", "RXRG", "PPARD", cell_counts_subclusters)

coexpression_PPARD_RXRG <- rbind(CM_coexpression_PPARD_RXRG, AD_coexpression_PPARD_RXRG, FB_coexpression_PPARD_RXRG) %>%
  mutate(dimer = "PPARD_RXRG")


coexpression_PPREs <- rbind(coexpression_PPARA_RXRA, coexpression_PPARA_RXRB, coexpression_PPARA_RXRG,
                            coexpression_PPARD_RXRA, coexpression_PPARD_RXRB, coexpression_PPARD_RXRG,
                            coexpression_PPARG_RXRA, coexpression_PPARG_RXRB, coexpression_PPARG_RXRG)

PPRE <- ggplot(coexpression_PPREs %>% filter(co_expr == "co_expressed"), aes(x = subcluster, y = n, fill = Condition)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_grid(cluster~dimer, scales = "free_x") +
  scale_fill_manual(values = colors_control_PKP2) +
  theme_pubr() +
  theme(panel.border = element_rect(color = "black", fill = NA),
        panel.grid.major = element_line(color = "lightgrey", linetype = 3),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    y = "Cell count with co-expression of a PPAR and RXR"
  )


PPRE
RXRA <- FeatureScatter(seurat_CM, feature1 = "HADH", feature2 = "RXRA", group.by = "Condition", split.by = "subcluster", cols = c("#008080", "#DC143C"), plot.cor = FALSE)
PPARA <- FeatureScatter(seurat_CM, feature1 = "HADH", feature2 = "PPARA", group.by = "Condition", split.by = "subcluster", cols = c("#008080", "#DC143C"), plot.cor = FALSE)

RXRA + PPARA + plot_layout(nrow = 2)


PPARA <- (CM_coexpression_PPARA_RXRA / CM_coexpression_PPARA_RXRB / CM_coexpression_PPARA_RXRG) + plot_layout(axis_titles = "collect", nrow = 3) + theme(legend.position = "top")
PPARD <- (CM_coexpression_PPARD_RXRA / CM_coexpression_PPARD_RXRB / CM_coexpression_PPARD_RXRG) + plot_layout(axis_titles = "collect", nrow = 3) + NoLegend()
PPARG <- (CM_coexpression_PPARG_RXRA / CM_coexpression_PPARG_RXRB / CM_coexpression_PPARG_RXRG) + plot_layout(axis_titles = "collect", nrow = 3) + NoLegend()

pdf("plots_excl_h55/PPAR_RXR_coexpression.pdf", height = 10, width = 14)
PPARA
PPARD
PPARG
dev.off()

## Single cell

print(merge_dat[merge_dat$gene%in%common[1:10],c('gene','p_val_sc','p_val_bulk')])

print(merge_dat[merge_dat$gene%in%common_FB[1:10],c('gene','p_val_sc','p_val_bulk')])

VlnPlot(subset_seurat_merge, features = c("FKBP5", "CD36", "PPARG", "LRP1", "SLC8A1"), idents = c("CM_Control", "CM_PKP2"), group.by = "orig.ident.condition", layer = "scale.data", alpha = 0.5)
VlnPlot(subset_seurat_merge, features = c("LRP1", "FKBP5", "PPARG"), idents = c("CM_Control", "CM_PKP2"), group.by = "orig.ident.condition", layer = "scale.data", alpha = 0.5)
VlnPlot(subset_seurat_merge, features = c("ADD1", "PER1", "CEBPB", "CEBPD", "PPARG"), idents = c("CM_Control", "CM_PKP2"), split.by = "Condition", group.by = "sub.CM", layer = "scale.data", alpha = 0.5) + RestoreLegend() 
VlnPlot(subset_seurat_merge, features = c("ADD1", "PER1", "CEBPB", "CEBPD", "PPARG"), idents = c("FB_Control", "FB_PKP2"), split.by = "Condition", group.by = "sub.FB", layer = "scale.data", alpha = 0.5) + RestoreLegend()

# Co expression of PKP2 with PPRE-genes
list_co_expression_PPRE_PKP2 <- list()

for (i in PPRE_geneset) {
  co_expression_PKP2_gene <- get_coexpression_count(subset_seurat_merge, "CM", "PKP2", i, cell_counts_subclusters) %>% mutate(co_expression = paste0("PKP2_", i))
  
  list_co_expression_PPRE_PKP2[[i]] <- co_expression_PKP2_gene 
}

coexpression_PKP2_PPRE_CM <- bind_rows(list_co_expression_PPRE_PKP2)

PPRE_PKP2_CM <- ggplot(coexpression_PKP2_PPRE_CM %>% filter(co_expr == "co_expressed", co_expression %in% c("PKP2_ACACB", "PKP2_CPT1A", "PKP2_DECR1", "PKP2_ETFDH", 
                                                                                                            "PKP2_FABP3", "PKP2_HADH", "PKP2_HADHA",
                                                                                                            "PKP2_HADHB", "PKP2_PPARA")),
                       aes(x = subcluster, y = ratio_co_exp, fill = Condition)) +
  geom_col(position = "dodge") +
  facet_grid(~co_expression, scales = "free_x") +
  scale_fill_manual(values = colors_control_PKP2) +
  theme_pubr() +
  theme(panel.border = element_rect(color = "black", fill = NA),
        panel.grid.major = element_line(color = "lightgrey", linetype = 3),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    y = "Cell count with co-expression of a PPAR and RXR"
  )

png("plots_excl_h55/cell_count_coexpression_PPRE_vs_PKP2.png", width = 1000, height = 700)
PPRE_PKP2_CM
dev.off()


get_correlation_coexpression <- function(seurat_object, gene_list, cluster) {
  list_co_expression <- list()
  for (i in gene_list) {
    Idents(seurat_object) <- "subcluster"
    co_expression_PKP2_gene <- FeatureScatter(subset(seurat_object, subset = seurat_annotations == cluster), "PKP2", i, split.by = "Condition")
    
    list_co_expression[[i]] <- co_expression_PKP2_gene 
    
  }
  return(list_co_expression)
}

PPREs <- c("PPARA", "PPARD", "PPARG", "RXRA", "RXRB", "RXRG")
  
PPRE_coexpression_PKP2 <- get_correlation_coexpression(subset_seurat_merge, c("HADH", "HADHB", "HADHA", "ACACB"), "CM")


ACACB_PKP2 <- PPRE_coexpression_PKP2[["ACACB"]][["data"]]
HADHA_PKP2 <- PPRE_coexpression_PKP2[["HADHA"]][["data"]]
HADH_PKP2 <- PPRE_coexpression_PKP2[["HADH"]][["data"]]
HADHB_PKP2 <- PPRE_coexpression_PKP2[["HADHB"]][["data"]]


ACACB <- ggplot(ACACB_PKP2, aes(x = PKP2, y = ACACB, color = Condition)) +
  geom_point() +
  geom_smooth(data = ACACB_PKP2 %>% filter(ACACB > 0, PKP2 > 0), color = "black", method = "lm") +
  stat_cor(label.y = 4) +
  facet_grid(Condition~colors, scales = "free") +
  scale_fill_manual(values = colors_control_PKP2) +
  scale_color_manual(values = colors_control_PKP2) +
  theme_pubr()

HADH <- ggplot(HADH_PKP2, aes(x = PKP2, y = HADH, color = Condition)) +
  geom_point() +
  geom_smooth(data = HADH_PKP2 %>% filter(HADH > 0, PKP2 > 0), color = "black", method = "lm") +
  stat_cor(label.y = 4) +
  facet_grid(Condition~colors, scales = "free") +
  scale_fill_manual(values = colors_control_PKP2) +
  scale_color_manual(values = colors_control_PKP2) +
  theme_pubr() +
  theme(legend.position = "none")

HADHA <- ggplot(HADHA_PKP2, aes(x = PKP2, y = HADHA, color = Condition)) +
  geom_point() +
  geom_smooth(data = HADHA_PKP2 %>% filter(HADHA > 0, PKP2 > 0), color = "black", method = "lm") +
  stat_cor(label.y = 4) +
  facet_grid(Condition~colors, scales = "free") +
  scale_fill_manual(values = colors_control_PKP2) +
  scale_color_manual(values = colors_control_PKP2) +
  theme_pubr() +
  theme(legend.position = "none")

HADHB <- ggplot(HADHB_PKP2, aes(x = PKP2, y = HADHB, color = Condition)) +
  geom_point() +
  geom_smooth(data = HADHB_PKP2 %>% filter(HADHB > 0, PKP2 > 0), color = "black", method = "lm") +
  stat_cor(label.y = 4) +
  facet_grid(Condition~colors, scales = "free") +
  scale_fill_manual(values = colors_control_PKP2) +
  scale_color_manual(values = colors_control_PKP2) +
  theme_pubr() +
  theme(legend.position = "none")

png("plots_excl_h55/co_expr_pkp2_HADHs_ACACB.png", height = 1600, width = 1000)
ACACB / HADH / HADHA / HADHB
dev.off()



##### STEP 5a.5: GSEA with clusterProfiler #####
# Install pathview and XML via command line
# gene list for GSEA requires ENSEMBL ID
# Set organism for GSEA
organism = "org.Hs.eg.db"
library(organism, character.only = TRUE)

gsea_cluster <- function(cluster_annotation) {
  
  pdf("plots/GSEA_dotplots.pdf", height = 12, width = 10)
  for (i in unique(cluster_annotation)) {
    if(i == "AD") {
      next
    }
    print(paste("Running GSEA on", i))
    bulk_DE <- get(paste0("bulk_", i, "_DE"))
    gene_list <- bulk_DE$avg_log2FC
    names(gene_list) <- bulk_DE$gene_symbol
    gene_list <- na.omit(gene_list)
    gene_list <- sort(gene_list, decreasing = TRUE)
    gse <- gseGO(geneList = gene_list,
                    ont = "ALL",
                    keyType = "SYMBOL",
                    minGSSize = 3,
                    maxGSSize = 800,
                    pvalueCutoff = 0.05,
                    verbose = TRUE,
                    OrgDb = organism,
                    pAdjustMethod = "bonferroni")
  
  require(DOSE)
  dplt <- dotplot(gse, showCategory = 10, split = ".sign", title = i) + 
    facet_grid(.~.sign)
  print(dplt)
  assign(paste0("gse_", i), gse, envir = .GlobalEnv)
  }
  dev.off()
}

gsea_cluster(cluster_annotation)





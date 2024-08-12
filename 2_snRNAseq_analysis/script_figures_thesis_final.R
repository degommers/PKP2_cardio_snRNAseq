#### Thesis figures ####
# See previous script: snRNAseq_analysis_PKP2_Control.R
# required to have loaded in environment:
# - cell_counts
# - subset_seurat_merge


# Subset CM
seurat_CM <- subset(subset_seurat_merge, subset = seurat_annotations == "CM")



#### Results: Figure 1: PPRE responsive genelist ####




#### Results: Figure 2 ####
## PCA plot ##
bulk <- AggregateExpression(subset_seurat_merge, group.by = c("Condition", "orig.ident"), return.seurat = TRUE)
tail(Cells(bulk))
bulk <- NormalizeData()
bulk <- FindVariableFeatures(bulk, selection.method = "vst", nfeatures = 2000)
bulk <- ScaleData(bulk, features = rownames(bulk))

bulk <- RenameCells(bulk[["RNA"]], new.names = paste0("bulk_", colnames(x = bulk[['RNA']])))

num_samples <- ncol(bulk)
num_genes <- nrow(bulk)
npcs_bulk <- min(50, num_samples - 1, num_genes - 1)
bulk <- RunPCA(bulk, features = VariableFeatures(object = bulk), npcs = npcs_bulk)

pca_figure <- DimPlot(bulk, shape.by = "Condition", reduction = "pca", group.by = "orig.ident", cols = c("#008080", "#DC143C", "#000080", "#228B22", "#FF4500", "darkgrey"), pt.size = 4) 
pca_figure <- pca_figure + labs(title = "") + theme(legend.position = "right") 



## Abundance plot ##
abundance <- bar_plot_count(cell_counts, fill = "seurat_annotations")
abundance_cell_count <- abundance + cellcounts_total + plot_layout(widths = c(15, 1))

abundance + NoLegend()

## Dot plot marker genes ##
order_dotplot <- c("AD", "NC", "T-cells", "Monocytes", "EC", "FB", "Mural", "CM", "SMC", "Myeloid", "Lymphoid", "Unknown")
dotplot_marker_genes <- DotPlot(subset_seurat_merge, features =  genes_cluster_annotation_paper_combined, group.by = "seurat_annotations", cols = "Spectral") + 
  scale_y_discrete(limits = order_dotplot) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.2, face = "italic", size = 18),
        axis.text.y = element_text(size = 16),
        axis.title = element_text(size = 16),
        # legend.position = "top",
        legend.title = element_text(size = 16),
        legend.title.position = "top",
        legend.text = element_text(size = 14),
        legend.box = "vertical")

pdf("plots_excl_h55/presentation/marker_genes.pdf", height = 14, width = 7)
dotplot_marker_genes
dev.off()
##### Boxplot cell abundance #####
boxplot_cells <- boxplot_cell_counts(cell_counts) + 
  scale_y_break(c(0.50, 0.90), scales = 0.2, ticklabels = c(0.9, 1.0)) + theme(legend.position = "right")

##### UMAP #####
umap_labeled <- DimPlot(subset_seurat_merge, reduction = "umap") + scale_color_d3(palette = "category20") + ylim(-15, 15) +
  theme(legend.position = "right") + guides(color = guide_legend(ncol = 1, override.aes = list(size = 5)))


Idents(subset_seurat_merge) <- "seurat_annotations"
umap_labeled_pres <- umap_labeled + theme(legend.text = element_text(size = 16),
                                          axis.text = element_text(size = 16),
                                          axis.title = element_text(size = 18))

##### PDF: snRNAseq analysis #####
pdf("plots_excl_h55/Figure1_ABC_3.pdf", width = 15, height = 15)
topplot <- ((dotplot_marker_genes) / (umap_labeled | abundance | pca_figure)) + plot_layout(heights = c(2, 3))
topplot / boxplot_cells + plot_layout(nrow = 2, heights = c(2, 1)) # add colors to boxplot in illustrator and add tags as well
dev.off()


#### Results: Figure 4 ####
pdf("plots_excl_h55/Fig3_correct_data.pdf", height = 50, width = 35)
Fig4 <- ggplot(expression_pkp2_PPRE_logN, aes(x = gene, y = expression, fill = Condition)) +
  geom_point(position = position_jitterdodge(dodge.width = 0.85, jitter.width = 0.70), aes(color = Condition), size = 1, alpha = 0.5) + 
  scale_color_manual(values =  c("#008080", "#DC143C")) +
  facet_grid2(seurat_annotations~., scales = "free", strip = strip_color_clusters) +
  theme_pubr() +
  ylim(0, 6.5) +
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


Fig4_CM <- ggplot(expression_pkp2_PPRE_logN %>% filter(seurat_annotations == "CM", gene != "PKP2"), aes(x = gene, y = expression, fill = Condition)) +
  geom_point(position = position_jitterdodge(dodge.width = 0.85, jitter.width = 0.70), aes(color = Condition), size = 1, alpha = 0.5) + 
  scale_color_manual(values =  c("#008080", "#DC143C")) +
  facet_grid2(seurat_annotations~., scales = "free", strip = strip_color_CM) +
  theme_pubr() +
  ylim(0, 6.5) +
  labs(
    x = "PPRE responsive genes",
    y = "LogNormalized expression"
  ) + 
  stat_compare_means(method = "wilcox.test", hide.ns = TRUE, label = "p.signif" , size = 8, vjust = 1.5)

pdf("plots_excl_h55/presentation/Fig3_CM_data.pdf", height = 10, width = 20)
Fig4_CM + theme(
  axis.text.x = element_text(angle = 90, hjust = 1, face = "italic", size = 15),
  axis.text.y = element_text(size = 15),
  axis.title.y = element_text(size = 18),
  axis.title.x = element_text(size = 18),
  legend.text = element_text(size = 18),
  legend.title = element_blank(),
  strip.text.y = element_text(size = 18),
  panel.spacing = unit(1.5, "lines")
) +
  guides(color = guide_legend(override.aes = list(size = 12)))
dev.off()



#### Results: Figure 6: Volcano Plot ####


#### Results: PKP2 expression
pal_cm2 <- c("#D3D3D3", "#D3D3D3", "#D9C7B6", "#D9C7B6", "#DFBB9A", "#DFBB9A",
             "#E5AF7E", "#E5AF7E", "#ECA362", "#ECA362", "#F29746", "#F29746",
             "#F88B2A", "#F88B2A", "#FF7F0E", "#FF7F0E")
print(pal_cm2)

VlnPlot(subset_seurat_merge, group.by = "orig.ident.condition", features = "PKP2", cols = c("#008080", "#DC143C", "#000080", "#228B22", "#FF4500", "darkgrey"))

Idents(seurat_CM) <- "subcluster.condition"
levels(seurat_CM)

# Relevel the seurat object
levels_sub_con <- c("PKP2-CM_0", "Control-CM_0", "PKP2-CM_1", "Control-CM_1", "PKP2-CM_2", "Control-CM_2",
                    "PKP2-CM_3", "Control-CM_3", "PKP2-CM_4", "Control-CM_4", "PKP2-CM_5", "Control-CM_5",
                    "PKP2-CM_6", "Control-CM_6", "PKP2-CM_7", "Control-CM_7")

subset_seurat_merge$subcluster.condition <- factor(subset_seurat_merge$subcluster.condition, levels = levels_sub_con)

PPRE_geneset_ordered <- c("ABCB11", "ACAD11", "FABP1", "LEP", "SLC27A2", "RXRG", "DECR2",
                          "ADIPOQ", "SESN2", "PPARG", "PPARD", "EHHADH", "SLC25A17", "CPT1A", "RXRA",
                          "RXRB", "SMIM37", "ACADS", "HADH", "ETFDH", "PPARA", "ETFDH", "PPARA",
                          "FABP3", "ACACB", "HADHA", "HADHB", "DECR1", "PKP2")


# comment: for glycolysis look at:  "SLC2A4", 
hm2 <- DoHeatmap(subset(subset_seurat_merge, subset = seurat_annotations == "CM", downsample = 1000), group.by = "subcluster.condition", 
                group.colors = pal_cm2, features = c(PPRE_geneset_ordered),
                raster = FALSE,
                size = 3) +
  scale_fill_gradientn(colors = c("darkblue", "white", "red"), limits = c(-3, 3))
pdf("plots_excl_h55/hm_PPRE_PKP2_sd.pdf", width = 20)
hm
dev.off()

#### Results Figure ?: Co expression PPRE dimer ####

rxra_ppara <- FeaturePlot(subset(subset_seurat_merge, subset = seurat_annotations == "CM"), features = c("PPARA", "RXRA"), 
                          slot = "data", blend = TRUE, split.by = "Condition",
                          pt.size = 0.5, alpha = 0.5)

rxrg_ppara <- FeaturePlot(subset(subset_seurat_merge, subset = seurat_annotations == "CM"), features = c("PPARA", "RXRG"), blend = TRUE, split.by = "Condition") + scale_fill_manual(values = "white")
rxrb_ppara <- FeaturePlot(seurat_CM, features = c("PPARA", "RXRB"), blend = TRUE, split.by = "Condition",
                          alpha = 0.5) 


rxra_pparg <- FeaturePlot(seurat_CM, features = c("PPARG", "RXRA"), blend = TRUE, split.by = "Condition",
                          pt.size = 0.5, alpha = 0.5)
rxrg_pparg <- FeaturePlot(seurat_CM, features = c("PPARG", "RXRG"), blend = TRUE, split.by = "Condition",
                          pt.size = 0.5, alpha = 0.5)
rxrb_pparg <- FeaturePlot(seurat_CM, features = c("PPARG", "RXRB"), blend = TRUE, split.by = "Condition",
                          pt.size = 0.5, alpha = 0.5)

rxra_ppard <- FeaturePlot(seurat_CM, features = c("PPARD", "RXRA"), blend = TRUE, split.by = "Condition",
                          pt.size = 0.5, alpha = 0.5)
rxrg_ppard <- FeaturePlot(seurat_CM, features = c("PPARD", "RXRG"), blend = TRUE, split.by = "Condition",
                          pt.size = 0.5, alpha = 0.5)
rxrb_ppard <- FeaturePlot(seurat_CM, features = c("PPARD", "RXRB"), blend = TRUE, split.by = "Condition",
                          pt.size = 0.5, alpha = 0.5)
pdf("plots_excl_h55/PPRE_co_expr.pdf", height = 35, width = 20)
rxra_ppara / rxrb_ppara/ rxra_pparg / rxrb_pparg / rxra_ppard / rxrb_ppard + plot_layout(nrow = 6)
dev.off()

#### Results: Figure 7: CM specific analysis ####
# A: UMAP
umap_CM <- DimPlot(subset(subset_seurat_merge, subset = seurat_annotations == "CM"), reduction = "umap", 
                   group.by = c("sub.CM"), shape.by = "Condition", pt.size = 2, label = TRUE,
                   label.size = 4 , repel = TRUE, alpha = 0.5, ) + scale_color_manual(values = pal_cm)
umap_CM + theme(legend.text = element_text(size = 16),
                legend.title = element_text(size = 18),
                axis.text = element_text(size = 16),
                axis.title = element_text(size = 18))

# B: Abundance
abundance_CM <- cell_counts_subcluster %>% filter(cluster == "CM") %>%
  bar_plot_count(fill = "subcluster") +
  scale_fill_manual(values = pal_cm)
  
dev.off()



# D: expression of HADH, HADHA, HADHB
CM_subcl_levels <- c("CM_0", "CM_1", "CM_2", "CM_3", "CM_4", "CM_5", "CM_6", "CM_7")

subset_seurat_merge$subcluster <- factor(x = subset_seurat_merge$subcluster, levels = CM_subcl_levels)

VlnPlot(subset(subset_seurat_merge, subset = seurat_annotations == "CM"), features = c("HADHB", "HADH", "HADHA"), 
        cols = c("#008080", "#DC143C", "grey"), split.by = "Condition", group.by = "subcluster",
        layer = "data", split.plot = FALSE,
        ncol = 1,
        alpha = 0.5) + RestoreLegend()

# Coexpression
FeatureScatter(seurat_CM, feature1 = "RXRA", feature2 = "PPARA", group.by = "subcluster", split.by = "Condition")

HADH_PPARA <- FeatureScatter(subset(subset_seurat_merge, subset = seurat_annotations == "CM"), feature1 = "PPARA", 
               feature2 = "HADH", split.by = "subcluster", group.by = "Condition",
               cols = alpha(c("#008080", "#DC143C"), 0.66)) + ggtitle("PPARA:HADH") +
  guides(color = guide_legend(override.aes = list(size = 5)))
RXRA_HADH <- FeatureScatter(subset(subset_seurat_merge, subset = seurat_annotations == "CM"), feature1 = "RXRA", 
               feature2 = "HADH", split.by = "subcluster", group.by = "Condition", 
               cols = alpha(c("#008080", "#DC143C"), 0.66)) + ggtitle("RXRA:HADH") + NoLegend()
PPARA_RXRA <- FeatureScatter(subset(subset_seurat_merge, subset = seurat_annotations == "CM"), feature1 = "PPARA", 
               feature2 = "RXRA", split.by = "subcluster", group.by = "Condition",
               cols = alpha(c("#008080", "#DC143C"), 0.66)) + ggtitle("PPARA:RXRA") + NoLegend()


pdf("plots_excl_h55/scatter_coexpr_PPARRXRHADH.pdf", width = 12, height = 8)
HADH_PPARA / RXRA_HADH / PPARA_RXRA
dev.off()


# E: PPRE expression (scaled expression, )
PPRE_expression_CM <- ggplot(expression_pkp2_PPRE_logN %>% filter(seurat_annotations == "CM"), aes(x = gene, y = expression, fill = Condition)) +
  geom_point(position = position_jitterdodge(dodge.width = 0.85, jitter.width = 0.70), aes(color = Condition), size = 1, alpha = 0.5) +
  scale_color_manual(values =  c("#008080", "#DC143C")) +
  facet_grid2(subcluster~., scales = "free", strip = strip_color_CM) +
  theme_pubr() +
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
    legend.position = "top",
    panel.spacing = unit(1.5, "lines")
  ) +
  stat_compare_means(method = "wilcox.test", hide.ns = TRUE, label = "p.signif" , size = 5, vjust = 1) +
  guides(color = guide_legend(override.aes = list(size = 5)))

figure7 <- (umap_CM | abundance_CM) / PPRE_expression_CM + plot_layout(heights = c(1, 4)) + plot_annotation(tag_levels = "A")

ggsave("plots_excl_h55/Figure7_CM_abundance_expression.png", plot = figure7, height = 16, width = 12)


# E: PPRE expression: Zoomed in expression of HADH/A/B
PPRE_expression_CM_H <- ggplot(expression_pkp2_PPRE_scale %>% filter(seurat_annotations == "CM", gene %in% c("HADHA", "HADH", "HADHB")), aes(x = gene, y = expression, fill = sample)) +
  geom_point(position = position_jitterdodge(dodge.width = 1, jitter.width = 0.2), aes(color = sample), size = 1, alpha = 0.5) +
  scale_color_manual(values = c("#008080", "#DC143C", "#000080", "#228B22", "#FF4500", "darkgrey")) +
  facet_grid2(subcluster~., scales = "free", strip = strip_color_CM) +
  theme_pubr() +
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
    legend.position = "right",
    panel.spacing = unit(1.5, "lines")
  ) +
  stat_compare_means(method = "wilcox.test", hide.ns = TRUE, label = "p.signif" , size = 5, vjust = 1) +
  guides(color = guide_legend(override.aes = list(size = 5)))

#### Results: Figure 8: Fibroblasts ####
# A: UMAP
umap_FB <- DimPlot(subset(subset_seurat_merge, subset = seurat_annotations == "FB"), reduction = "umap", 
                   group.by = c("sub.FB"), shape.by = "Condition", pt.size = 2, label = TRUE,
                   label.size = 4 , repel = TRUE, alpha = 0.5, ) + scale_color_manual(values = pal_fb)
umap_FB + theme(legend.text = element_text(size = 16),
                legend.title = element_text(size = 18),
                axis.text = element_text(size = 16),
                axis.title = element_text(size = 18))

# B: Abundance
Idents(subset_seurat_merge)
abundance_FB <- cell_counts_subclusters %>% filter(cluster == "FB") %>%
  bar_plot_count(fill = "subcluster") +
  scale_fill_manual(values = pal_fb)


# C: Top expressed genes
Idents(subset_seurat_merge) <- "subcluster"


# E: Expression of PPRE-rG in FB
PPRE_expression_FB <- ggplot(expression_pkp2_PPRE_logN %>% filter(seurat_annotations == "FB"), aes(x = gene, y = expression, fill = Condition)) +
  geom_point(position = position_jitterdodge(dodge.width = 0.85, jitter.width = 0.70), aes(color = Condition), size = 1, alpha = 0.5) +
  scale_color_manual(values =  c("#008080", "#DC143C")) +
  facet_grid2(subcluster~., scales = "free", strip = strip_color_FB) +
  theme_pubr() +
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
    legend.position = "top",
    panel.spacing = unit(1.5, "lines")
  ) +
  stat_compare_means(method = "wilcox.test", hide.ns = TRUE, label = "p.signif" , size = 5, vjust = 1) +
  guides(color = guide_legend(override.aes = list(size = 5)))

# Combined plot Figure 8

figure8 <- (umap_FB + abundance_FB) / PPRE_expression_FB + plot_layout(heights = c(1, 4)) + plot_annotation(tag_levels = "A")
ggsave("plots_excl_h55/Figure8_FB_abundance_expression.png", plot = figure8, height = 18, width = 12)

#### Results: Figure 9: AD specific analysis ####
# A: UMAP
umap_AD <- DimPlot(subset(subset_seurat_merge, subset = seurat_annotations == "AD"), reduction = "umap", 
                   group.by = c("sub.AD"), shape.by = "Condition", pt.size = 2, label = TRUE,
                   label.size = 4 , repel = TRUE, alpha = 0.5, ) + scale_color_manual(values = pal_ad)
umap_AD + theme(legend.text = element_text(size = 16),
                legend.title = element_text(size = 18),
                axis.text = element_text(size = 16),
                axis.title = element_text(size = 18))

# B: Abundance
Idents(subset_seurat_merge)
abundance_AD <- cell_counts_subclusters %>% filter(cluster == "AD") %>%
  bar_plot_count(fill = "subcluster") +
  scale_fill_manual(values = pal_ad)

# C: Top expressed genes
Idents(subset_seurat_merge) <- "subcluster"

# E: Expression of PPRE-rG in FB
PPRE_expression_AD <- ggplot(expression_pkp2_PPRE_logN %>% filter(seurat_annotations == "AD"), aes(x = gene, y = expression, fill = Condition)) +
  geom_point(position = position_jitterdodge(dodge.width = 0.85, jitter.width = 0.70), aes(color = Condition), size = 1, alpha = 0.5) +
  scale_color_manual(values =  c("#008080", "#DC143C")) +
  facet_grid2(subcluster~., scales = "free", strip = strip_color_AD) +
  theme_pubr() +
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
    legend.position = "top",
    panel.spacing = unit(1.5, "lines")
  ) +
  stat_compare_means(method = "wilcox.test", hide.ns = TRUE, label = "p.signif" , size = 5, vjust = 1) +
  guides(color = guide_legend(override.aes = list(size = 5)))

# Combined plot Figure 8

figure9 <- (umap_AD + abundance_AD) / PPRE_expression_AD + plot_layout(heights = c(1, 4)) + plot_annotation(tag_levels = "A")
ggsave("plots_excl_h55/Figure8_AD_abundance_expression.png", plot = figure9, height = 16, width = 12)




#### Results: Figure x: CM AD FB

cms_plot <- umap_CM + abundance_CM
fbs_plot <- umap_FB + abundance_FB
ads_plot <- umap_AD + abundance_AD

figurex <- cms_plot / fbs_plot / ads_plot + plot_annotation(tag_levels = "A")

pdf("plots_excl_h55/figurex_umap_abudance_all.pdf", height = 20, width = 15)
figurex
dev.off()


png("plots_excl_h55/figurex_umap_abudance_all.png", height = 1000, width = 800)
figurex
dev.off()

#### Results: Figure 10: top10 expressed genes per subcluster ####
# C: Top expressed genes
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
top10_genes_matrix_2 <- top10_genes_matrix[ rowSums(top10_genes_matrix !=0) < 2, ]



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

# Colors
custom_colors <- c("white", colorRampPalette(c("blue", "red"))(10))
ann_colors <- list(
  Condition = c(PKP2 = "#DC143C", Control = "#008080"),
  Cluster = c(CM = "#FF7F0EFF", FB = "#1F77B4FF", AD = "#7F7F7FFF")
)

# Heatmap with top genes that overlap in at least 2 subclusters
heatmap_top10_genes <- pheatmap(top10_genes_matrix_5, 
                                color = custom_colors,
         fontsize_row = 12,
         fontsize_col = 12,
         annotation_row = presence,
         annotation_col = annotation,
         annotation_colors = ann_colors,
         border_color = "lightgrey"
         )

pdf("plots_excl_h55/heatmap_top10_genes_overlap.pdf", height = 16, width = 18)
heatmap_top10_genes
dev.off()

png("plots_excl_h55/heatmap_top10_genes_overlap.png", height = 900, width = 900)
heatmap_top10_genes
dev.off()


# Heatmap with top genes that are unique per cluster
heatmap_top10_genes_2 <- pheatmap(top10_genes_matrix_2,
         fontsize_row = 3,
         annotation_col = annotation_2,
         annotation_colors = ann_colors,
         border_color = "lightgrey"
         )


pdf("plots_excl_h55/heatmap_top10_genes_unique.pdf", height = 16, width = 18)
heatmap_top10_genes_2
dev.off()

#### Results: Figure 11: Percent expressed per subcluster ####
perc_stats_AD_FB_CM <- percent_stats %>%
  select(contains(c("CM", "AD", "FB"))) %>%
  as.matrix()

## Extract condition and cluster from column names
condition <- sapply(strsplit(colnames(perc_stats_AD_FB_CM), "-|_"), `[`, 1)
cluster <- sapply(strsplit(colnames(perc_stats_AD_FB_CM), "-|_"), `[`, 2)

# Create annotation data frame
annotation <- data.frame(
  Condition = condition,
  Cluster = cluster
)
rownames(annotation) <- colnames(top10_genes_matrix_5)

percent_expressed_plot <- pheatmap(perc_stats_AD_FB_CM,
         annotation_col = annotation,
         annotation_colors = ann_colors,
         color = c("lightyellow", muted("purple")),
         border_color = "lightgrey")

pdf("plots_excl_h55/percent_expressed_plot_AD_FB_CM.pdf", height = 16, width = 21)
percent_expressed_plot
dev.off()

#### Supplementary ####
##### Absolute cell counts #####
pdf("plots_excl_h55/Supplementary_fig_absolute_cell_counts.pdf", width = 7, height = 8)
ggplot(cell_counts, aes(x = sample, y = cell_count, fill = cluster)) +
  geom_col(position = "dodge") +
  scale_fill_d3(palette = "category20") +
  theme_pubclean() +
  labs(
    y = "cell count"
  )
dev.off()

table_cell_counts <- cell_counts %>% group_by(cluster, Condition) %>% summarize(cell_count = sum(cell_count)), n=30)

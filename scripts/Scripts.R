library(Seurat)
library(ggplot2)
library(SingleR)
library(dplyr)
library(celldex)
library(RColorBrewer)
library(SingleCellExperiment)
library(tidyverse)
library(gridExtra)
library(SeuratObject)
library(remotes)
library(SeuratData)
library(scales)
library(scater)
library(SingleCellExperiment)
library(patchwork)
library(harmony)
library(scDblFinder)
library(garnett)
library(monocle3)

inf_csf = readRDS(file = '/Users/tee/ms_rds/tot_harmony_csf.Rds')
inf_pbmc = readRDS(file = '/Users/tee/ms_rds/tot_harmony_pbmc.Rds')
inf_pbmc = readRDS(file = '/Users/tee/ms/inf_harmony_SCT_pbmc1')
SetIdent(inf_csf, value = 'orig.ident')
Idents(inf_csf)
Idents(inf_csf) = 'orig.ident'
View(inf_csf@meta.data)
meta = inf_csf@meta.data
dim(meta)
summary(meta$nCount_RNA)
summary(meta$nFeature_RNA)
summary(meta$percent.mt)
inf_csf[["percent.rb"]] <- PercentageFeatureSet(inf_csf, pattern = "^RP[SL]")
FeatureScatter(inf_csf, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(inf_csf, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
FeatureScatter(inf_csf, feature1 = "nCount_RNA", feature2 = "percent.rb")
FeatureScatter(inf_csf, feature1 = "percent.rb", feature2 = "percent.mt")

inf_csf2 = subset(inf_csf, subset = nFeature_RNA > 500 & percent.mt < 15)
meta2 = inf_csf2@meta.data
dim(meta2)
summary(meta2$nCount_RNA)
summary(meta2$nFeature_RNA)
summary(meta2$percent.mt)
FeatureScatter(inf_csf2, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(inf_csf2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
FeatureScatter(inf_csf2, feature1 = "nCount_RNA", feature2 = "percent.rb")
FeatureScatter(inf_csf2, feature1 = "percent.rb", feature2 = "percent.mt")
inf_csf2 = NormalizeData(inf_csf2)
inf_csf2 = FindVariableFeatures(inf_csf2, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(inf_csf2), 10)
top10 
plot1 <- VariableFeaturePlot(inf_csf2)
LabelPoints(plot = plot1, points = top10, repel = TRUE, xnudge = 0, ynudge = 0)
all.genes <- rownames(inf_csf2)
inf_csf2 <- ScaleData(inf_csf2, features = all.genes)
inf_csf2 <- RunPCA(inf_csf2)
VizDimLoadings(inf_csf2, dims = 1:9, reduction = "pca") & 
  theme(axis.text=element_text(size=5), axis.title=element_text(size=8,face="bold"))
DimHeatmap(inf_csf2, dims = 1:6, nfeatures = 20, cells = 500, balanced = T)
DimPlot(inf_csf2, reduction = "pca")
ElbowPlot(inf_csf2)



############################
inf_csf <- inf_csf %>% 
  RunUMAP(reduction = "harmony", dims = 1:20, verbose = F) %>% 
  FindNeighbors(reduction = "harmony", k.param = 10, dims = 1:20) %>% 
  FindClusters() %>% 
  identity()
inf_csf <- SetIdent(inf_csf,value = "orig.ident")
DimPlot(inf_csf,reduction = "harmony") + 
  plot_annotation(title = "after integration (Harmony)")
DimPlot(inf_csf,reduction = "umap") + 
  plot_annotation(title = "after integration (Harmony)")

s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes

inf_csf <- JoinLayers(inf_csf, overwrite = TRUE)
DefaultAssay(inf_csf) = "RNA"
inf_csf <- CellCycleScoring(inf_csf, s.features = s.genes, g2m.features = g2m.genes)
table(inf_csf[[]]$Phase)

options(future.globals.maxSize = 8000 * 1024^2)
inf_csf <- SCTransform(inf_csf, ncells = 8824, method = "glmGamPoi",
  vars.to.regress = c("percent.mt", "nFeature_RNA"), verbose = F)
inf_csf

inf_csf <- RunPCA(inf_csf, verbose = F)
inf_csf <- RunUMAP(inf_csf, dims = 1:30, verbose = F)
inf_csf <- FindNeighbors(inf_csf, dims = 1:30, verbose = F)
inf_csf <- FindClusters(inf_csf, verbose = F)
table(inf_csf[[]]$seurat_clusters)


# note that you can chain multiple commands together with %>%
inf_csf <- inf_csf %>%
  RunPCA() %>%
  FindNeighbors(dims = 1:30) %>%
  FindClusters() %>%
  RunUMAP(dims = 1:30)
DimPlot(object = inf_csf, reduction = "umap", label = T, split.by = 'orig.ident')
DimPlot(object = inf_csf, reduction = "umap", label = T)

DefaultAssay(inf_csf) = 'RNA'

######TOP FEATURES CSF
VizDimLoadings(inf_csf, dims = 1:5, reduction = "pca", ncol=3)
LabelPoints(
  VariableFeaturePlot(inf_csf),
  points = head(VariableFeatures(inf_csf),10),
  repel = TRUE) + 
  theme(legend.position="bottom")

#######TOP FEATURES PBMC
VizDimLoadings(inf_pbmc, dims = 1:3, reduction = "pca", ncol=3)
LabelPoints(
  VariableFeaturePlot(inf_pbmc),
  points = head(VariableFeatures(inf_pbmc),10),
  repel = TRUE) + 
  theme(legend.position="bottom")

######FINAL_AUTOMATIC ANNOTATION
inf_csf_new <- JoinLayers(inf_csf, overwrite = TRUE)
DefaultAssay(inf_csf) = "RNA"
inf_cts = GetAssayData(inf_csf_new, assay = "RNA", slot = "data")
View(inf_csf@meta.data)
Seurat::FeaturePlot(inf_csf, "GNLY") + 
scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral"))) + 
  ggtitle("GNLY")
Seurat::FeaturePlot(inf_csf, "CD74") + 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral"))) + 
  ggtitle("CD74")
FeaturePlot(inf_csf,"LYZ") + 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral"))) + 
  ggtitle("LYZ: monocytes")
FeaturePlot(inf_csf,"NKG7") + 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral"))) + 
  ggtitle("NKG7: natural killers")
FeaturePlot(inf_csf,"HBA1") + 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral"))) + 
  ggtitle("NKG7: natural killers")
tcell_genes <- c("IL7R", "LTB", "TRAC", "CD3D")
monocyte_genes <- c("CD14", "CST3", "CD68", "CTSS")
Seurat::FeaturePlot(inf_csf, tcell_genes, ncol=2)
Seurat::VlnPlot(inf_csf,
                features = tcell_genes,
                ncol = 2)
Seurat::FeaturePlot(inf_csf, monocyte_genes, ncol=2)
Seurat::VlnPlot(inf_csf,
                features = monocyte_genes,
                ncol = 2)
inf_csf <- Seurat::AddModuleScore(inf_csf,
          features = list(tcell_genes),
          name = "tcell_genes")
Seurat::DefaultAssay(inf_csf) <- "RNA"
inf_csf <- Seurat::SetIdent(inf_csf, value = inf_csf$SCT_snn_res.0.8)
Seurat::FeaturePlot(inf_csf, "tcell_genes1")
Seurat::VlnPlot(inf_csf, "tcell_genes1")
hpca = celldex::HumanPrimaryCellAtlasData()
imce = celldex::DatabaseImmuneCellExpressionData()
nov <- celldex::NovershternHematopoieticData()
com = SingleR(test = inf_cts, ref = list(hpca, imce, nov), 
              labels = list(hpca$label.main, imce$label.main, nov$label.main))
table(com$labels)
plotScoreHeatmap(com)
SingleR::plotDeltaDistribution(com)
singleR_labels <- com$labels
t <- table(singleR_labels)
other <- names(t)[t < 10]
singleR_labels[singleR_labels %in% other] <- "none"
inf_csf$SingleR_annot <- singleR_labels
dittoSeq::dittoDimPlot(inf_csf, "SingleR_annot", size = 0.7, split.by = 'orig.ident')
dittoSeq::dittoBarPlot(inf_csf, var = "SingleR_annot", group.by = "orig.ident")
dittoSeq::dittoBarPlot(inf_csf, 
                       var = "SingleR_annot", 
                       group.by = "SCT_snn_res.0.8")
inf_csf@meta.data$com.label = com$labels
inf_csf@meta.data$com.pruned = com$pruned.labels
inf_csf <- SetIdent(inf_csf, value = "SingleR_annot")
a = DimPlot(inf_csf, label = T , repel = T, label.size = 4) + NoLegend()
b = DimPlot(inf_csf, label = T , repel = T, label.size = 4, split.by = 'orig.ident') +
  NoLegend()
DimPlot(inf_csf, label = T , repel = T, label.size = 4, split.by = 'orig.ident') +
NoLegend()
ggsave(a, filename = 'csf_annotation.jpeg', device = 'jpeg', 
       width = 15, height = 10, units = 'in')
ggsave(b, filename = 'csf_annotation_split.jpeg', device = 'jpeg', 
       width = 15, height = 10, units = 'in')
saveRDS(inf_csf, file = "inf_harmony_SCT")

######PBMCS
########final
inf_pbmc[["percent.rb"]] <- PercentageFeatureSet(inf_pbmc, pattern = "^RP[SL]")
options(future.globals.maxSize = 8000 * 1024^2)
inf_pbmc <- SCTransform(inf_pbmc, ncells = 8824, method = "glmGamPoi",
    vars.to.regress = c("percent.mt", "nFeature_RNA"), verbose = F)
inf_pbmc
# note that you can chain multiple commands together with %>%
inf_pbmc <- inf_pbmc |> 
  RunPCA() |>
  FindNeighbors(dims = 1:30) |>
  FindClusters() |>
  RunUMAP(dims = 1:30)
DefaultAssay(inf_pbmc) = "RNA"
inf_pbmc_new <- JoinLayers(inf_pbmc, overwrite = TRUE)
inf_pbmc_cts = GetAssayData(inf_pbmc_new, assay = "RNA", slot = "data")
hpca = celldex::HumanPrimaryCellAtlasData()
imce = celldex::DatabaseImmuneCellExpressionData()
nov <- celldex::NovershternHematopoieticData()
com2 = SingleR(test = inf_pbmc_cts, ref = list(hpca, imce, nov), 
               labels = list(hpca$label.main, imce$label.main, nov$label.main))
saveRDS(com2, file = 'com2.Rds')
table(com2$labels)
plotScoreHeatmap(com2)  
inf_pbmc@meta.data$com.label = com2$labels
inf_pbmc@meta.data$com.pruned = com2$pruned.labels
inf_pbmc <- SetIdent(inf_pbmc, value = "com.pruned")
c = DimPlot(inf_pbmc, label = T , repel = T, label.size = 4, raster = F) +
  NoLegend()
d = DimPlot(inf_pbmc, label = T , repel = T, label.size = 4, raster = F,
            split.by = 'orig.ident') + NoLegend()
ggsave(c, filename = 'pbmc_cellular_annotation.jpeg', device = 'jpeg', width = 15,
       height = 10, units = 'in')
ggsave(d, filename = 'pbmc_cellular_annotation_split.jpeg', device = 'jpeg', 
       width = 15, height = 10, units = 'in')
table(com$labels)
plotScoreHeatmap(com)
View(inf_pbmc@meta.data)
SingleR::plotDeltaDistribution(com2)
singleR_labels <- com2$labels
t2 <- table(singleR_labels)
other <- names(t2)[t2 < 10]
singleR_labels[singleR_labels %in% other] <- "none"
inf_pbmc$SingleR_annot <- singleR_labels
dittoSeq::dittoDimPlot(inf_pbmc, "SingleR_annot", size = 0.7, split.by = 'orig.ident')
dittoSeq::dittoBarPlot(inf_pbmc, var = "SingleR_annot", group.by = "orig.ident")
dittoSeq::dittoBarPlot(inf_pbmc, 
                       var = "SingleR_annot", 
                       group.by = "SCT_snn_res.0.8")
inf_pbmc@meta.data$com.label = com2$labels
inf_pbmc@meta.data$com.pruned = com2$pruned.labels
inf_pbmc <- SetIdent(inf_pbmc, value = "SingleR_annot")
c = DimPlot(inf_pbmc, label = T , repel = T, label.size = 4,raster = F) + NoLegend()
d = DimPlot(inf_pbmc, label = T , repel = T, label.size = 4, split.by = 'orig.ident', 
  raster = F) +
  NoLegend()
DimPlot(inf_csf, label = T , repel = T, label.size = 4, split.by = 'orig.ident') +
  NoLegend()
ggsave(c, filename = 'pbmc_annotation.jpeg', device = 'jpeg', 
       width = 15, height = 10, units = 'in')
ggsave(d, filename = 'pbmc_annotation_split.jpeg', device = 'jpeg', 
       width = 15, height = 10, units = 'in')
saveRDS(inf_pbmc, file = "inf_harmony_SCT_pbmc2")
inf_pbmc = readRDS("inf_harmony_SCT_pbmc2.Rds")
library(monocle3)
inf_pbmc = readRDS("inf_harmony_SCT_pbmc2.Rds")
com2 = readRDS('/Users/tee/ms/com2.Rds')
View(inf_pbmc@meta.data)
inf_pbmc_new <- JoinLayers(inf_pbmc, overwrite = TRUE)
rm(inf_pbmc)
######################
cds_obj = readRDS('/Users/tee/ms/cds.Rds')
cds_obj <- SeuratWrappers::as.cell_data_set(inf_pbmc_new)
cds_obj = monocle3::preprocess_cds(cds_obj)
monocle3::plot_pc_variance_explained(cds_obj)
cds_obj<- monocle3::reduce_dimension(cds_obj, reduction_method = "UMAP")
cds_obj<- monocle3::reduce_dimension(cds_obj, reduction_method = "tSNE")
monocle3::plot_cells(cds_obj, label_cell_groups = T, 
      color_cells_by="SingleR_annot")
monocle3::plot_cells(cds_obj, label_cell_groups = T, reduction_method='UMAP',
    color_cells_by="SingleR_annot", group_label_size = 3)
cds_obj <- cluster_cells(cds_obj, resolution=1e-5)
monocle3::plot_cells(cds_obj, label_cell_groups = F, show_trajectory_graph = F, 
                     color_cells_by = "partition")
ws = plot_cells(cds_obj, color_cells_by="SingleR_annot", group_cells_by="partition", 
        labels_per_group = 2, label_cell_groups = T, group_label_size = 4)
plot_cells(cds_obj, color_cells_by="SingleR_annot", group_cells_by="partition", 
           group_label_size = 4)
plot_cells(cds_obj, color_cells_by="SingleR_annot", 
           label_groups_by_cluster=FALSE)
ggsave(ws, filename = 'pbmc_annotation_monocle.jpeg', device = 'jpeg', 
       width = 15, height = 10, units = 'in')

## Step 5: Learn a graph
cds_obj <- learn_graph(cds_obj, use_partition=FALSE, close_loop=FALSE)
#cds_obj <- learn_graph(cds_obj)
p1 <- plot_cells(cds_obj, color_cells_by="SingleR_annot",
  group_label_size=4, graph_label_size=3,
  label_cell_groups=FALSE, label_principal_points=TRUE,
  label_groups_by_cluster=FALSE)
p1
ggsave(p1, filename = 'pbmc_graph.jpeg', device = 'jpeg', 
       width = 15, height = 10, units = 'in')

cds_obj <- order_cells(cds_obj, root_pr_nodes = "Y_1")
## Step 6: Order cells
cds_obj <- order_cells(cds_obj)
plot_cells(cds_obj, color_cells_by = "pseudotime", label_branch_points = F, label_leaves = F)
rowData(cds_obj)$gene_name = rownames(cds_obj)
rowData(cds_obj)$gene_short_name = rowData(cds_obj)$gene_name
saveRDS(cds_obj, file = "cds.Rds")
View(cds_obj@colData@metadata)
p2 = plot_cells(cds_obj, genes = c("ISG15", "TNFRSF18", "TNFRSF4", "M1B2"), 
           label_cell_groups = F, 
           show_trajectory_graph = F)
ggsave(p2, filename = 'pbmc_graph_genes.jpeg', device = 'jpeg', 
       width = 15, height = 10, units = 'in')
cd_graph = graph_test(cds_obj, neighbor_graph = 'principal_graph', cores = 8)

cds = readRDS(file = "cds.Rds")
cdg = readRDS(file = "cdgraph.Rds")
colData(cds)
tcell_cds <- cds[,grepl("T_cells", colData(cds)$ident, ignore.case=TRUE)]
plot_cells(tcell_cds, color_cells_by="partition", label_branch_points = F, label_leaves = F)
tcell_res <- graph_test(tcell_cds, neighbor_graph="knn", cores=8)
tcell_res = na.omit(tcell_res)
tcell_res = tcell_res[order(-tcell_res$morans_test_statistic),]
tcell_ids <- row.names(subset(tcell_res, q_value < 0.05))
gene_module_df <- find_gene_modules(tcell_cds[tcell_ids,], resolution=1e-2)
unique(gene_module_df$supermodule)
cell_group_df <- tibble::tibble(cell=row.names(colData(tcell_cds)), 
                                cell_group=partitions(cds)[colnames(tcell_cds)])
agg_mat <- aggregate_gene_expression(tcell_cds, gene_module_df, cell_group_df)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
colnames(agg_mat) <- stringr::str_c("Partition ", colnames(agg_mat))
pheatmap::pheatmap(agg_mat, cluster_rows=TRUE, cluster_cols=TRUE,
                   scale="column", clustering_method="ward.D2",
                   fontsize=6)
rowData(tcell_cds)$gene_name = rownames(tcell_cds)
rowData(tcell_cds)$gene_short_name = rowData(tcell_cds)$gene_name
plot_cells(tcell_cds, 
           genes=gene_module_df %>% filter(module %in% c(8, 28, 33, 37)),
           group_cells_by="partition",
           color_cells_by="partition",
           show_trajectory_graph=FALSE)
saveRDS(tcell_cds, file = "tcell_res.Rds")

############

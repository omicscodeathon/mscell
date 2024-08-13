# script to integrate scRNA-Seq datasets


# load libraries
library(Seurat)
library(ggplot2)
library(tidyverse)
library(gridExtra)


# get data location
dir_prj <- list.dirs(path = 'dat/CSF', recursive = F, full.names = F)


for(x in dir_prj){
  name <- gsub(".*_", "", x)
  cts <- ReadMtx(mtx = paste0('dat/CSF/',x,'/matrix.mtx'),
           features = paste0('dat/CSF/',x,'/genes.tsv'),
           cells = paste0('dat/CSF/',x,'/barcodes.tsv'))
  assign(name, CreateSeuratObject(counts = cts))
}

ls()
rm(list = "name")
# merge datasets

merged_seurat <- merge(`MS-MA`, y = c(`MS-MB`, `MS-MC`, `MS-MD`, `MS-ME`, 
   `NMS-NMSA1`, `NMS-NMSB2`, `NMS-NMSC1`, `NMS-NMSD2`, `NMS-NMSE2`), 
   add.cell.ids = ls()[3:12])


merged_seurat

# QC & filtering -----------------------

View(merged_seurat@meta.data)
# create a sample column
merged_seurat$sample <- rownames(merged_seurat@meta.data)

# split sample column
merged_seurat@meta.data <- separate(merged_seurat@meta.data, col = 'sample', into = c('Patient', 'Barcode'), 
                                    sep = '_')
merged_seurat@meta.data <- separate(merged_seurat@meta.data, col = 'Patient', into = c('Type', 'Number'), 
                                    sep = '-')
unique(merged_seurat$Number)
unique(merged_seurat$Type)
# calculate mitochondrial percentage
merged_seurat$mitoPercent <- PercentageFeatureSet(merged_seurat, pattern='^MT-')

# explore QC


# filtering
merged_seurat_filtered <- subset(merged_seurat, subset = nCount_RNA > 500 &
                                   nFeature_RNA > 500 &
                                   mitoPercent < 10)

 merged_seurat_filtered

merged_seurat


# perform standard workflow steps to figure out if we see any batch effects --------
merged_seurat_filtered <- NormalizeData(object = merged_seurat_filtered)
merged_seurat_filtered <- FindVariableFeatures(object = merged_seurat_filtered)
merged_seurat_filtered <- ScaleData(object = merged_seurat_filtered)
merged_seurat_filtered <- RunPCA(object = merged_seurat_filtered)
ElbowPlot(merged_seurat_filtered)
merged_seurat_filtered <- FindNeighbors(object = merged_seurat_filtered, dims = 1:20)
merged_seurat_filtered <- FindClusters(object = merged_seurat_filtered)
merged_seurat_filtered <- RunUMAP(object = merged_seurat_filtered, dims = 1:20)


# plot
p1 <- DimPlot(merged_seurat_filtered, reduction = 'umap', group.by = 'Patient')
p2 <- DimPlot(merged_seurat_filtered, reduction = 'umap', group.by = 'Type',
              cols = c('red','green','blue'))

grid.arrange(p1, p2, ncol = 2, nrow = 2)


# perform integration to correct for batch effects ------
obj.list <- SplitObject(merged_seurat_filtered, split.by = 'Patient')
for(i in 1:length(obj.list)){
  obj.list[[i]] <- NormalizeData(object = obj.list[[i]])
  obj.list[[i]] <- FindVariableFeatures(object = obj.list[[i]])
}


# select integration features
features <- SelectIntegrationFeatures(object.list = obj.list)

# find integration anchors (CCA)
anchors <- FindIntegrationAnchors(object.list = obj.list,
                                  anchor.features = features)

# integrate data
seurat.integrated <- IntegrateData(anchorset = anchors)


# Scale data, run PCA and UMAP and visualize integrated data
seurat.integrated <- ScaleData(object = seurat.integrated)
seurat.integrated <- RunPCA(object = seurat.integrated)
seurat.integrated <- RunUMAP(object = seurat.integrated, dims = 1:50)


p3 <- DimPlot(seurat.integrated, reduction = 'umap', group.by = 'Patient')
p4 <- DimPlot(seurat.integrated, reduction = 'umap', group.by = 'Type',
              cols = c('red','green','blue'))


grid.arrange(p1, p2, p3, p4, ncol = 2, nrow = 2)









library(SeuratObject)
library(remotes)
library(SeuratData)
library(scales)
library(scater)
library(patchwork)
library(SingleCellExperiment)
library(scDblFinder)

####load in files 
fin1 = list.files(path = "/Users/tee/Downloads/scr/GSE133028/csf/5/H/", full.names = F, 
                  recursive = F)
for (i in fin1){
  names = gsub(pattern = "_t[0-9].csv", replacement = "", i)
  print(names)
  cts = read.csv(paste0("/Users/tee/Downloads/scr/GSE133028/csf/5/H/", i))
  assign(names, CreateSeuratObject(counts = cts, min.cells = 3, min.features = 200, 
                                   project = "NMS"))
}

fin2 = list.files(path = "/Users/tee/Downloads/scr/GSE133028/csf/5/M/", full.names = F,
                  recursive = F)
for (i in fin2){
  names = gsub(pattern = "_t[0-9].csv", replacement = "", i)
  print(names)
  cts = read.csv(paste0("/Users/tee/Downloads/scr/GSE133028/csf/5/M/", i))
  assign(names, CreateSeuratObject(counts = cts, min.cells = 3, min.features = 200, 
                                   project = "MS"))
}

msa = Read10X(data.dir = '/Users/tee/ms/dat/CSF/MS-MA')
msb = Read10X(data.dir = '/Users/tee/ms/dat/CSF/MS-MB')
msc = Read10X(data.dir = '/Users/tee/ms/dat/CSF/MS-MC')
msd = Read10X(data.dir = '/Users/tee/ms/dat/CSF/MS-MD')
mse = Read10X(data.dir = '/Users/tee/ms/dat/CSF/MS-ME')
nmsa = Read10X(data.dir = '/Users/tee/ms/dat/CSF/NMS-NMSA1')
nmsb = Read10X(data.dir = '/Users/tee/ms/dat/CSF/NMS-NMSB2')
nmsc = Read10X(data.dir = '/Users/tee/ms/dat/CSF/NMS-NMSC1')
nmsd = Read10X(data.dir = '/Users/tee/ms/dat/CSF/NMS-NMSD2')
nmse = Read10X(data.dir = '/Users/tee/ms/dat/CSF/NMS-NMSE2')
msa = CreateSeuratObject(counts=msa, min.cells = 3, min.features = 200, project = 'MS')
msb = CreateSeuratObject(counts=msb, min.cells = 3, min.features = 200, project = 'MS')
msc = CreateSeuratObject(counts=msc, min.cells = 3, min.features = 200, project = 'MS')
msd = CreateSeuratObject(counts=msd, min.cells = 3, min.features = 200, project = 'MS')
mse = CreateSeuratObject(counts=mse, min.cells = 3, min.features = 200, project = 'MS')
nmsa = CreateSeuratObject(counts=nmsa, min.cells = 3, min.features = 200, project = 'NMS')
nmsb = CreateSeuratObject(counts=nmsb, min.cells = 3, min.features = 200, project = 'NMS')
nmsc = CreateSeuratObject(counts=nmsc, min.cells = 3, min.features = 200, project = 'NMS')
nmsd = CreateSeuratObject(counts=nmsd, min.cells = 3, min.features = 200, project = 'NMS')
nmse = CreateSeuratObject(counts=nmse, min.cells = 3, min.features = 200, project = 'NMS')

#####QC
## Mitochondrion Percentage
msa = PercentageFeatureSet(object = msa, pattern = "^MT-", col.name = 'percent.mt')
msb = PercentageFeatureSet(object = msb, pattern = "^MT-", col.name = 'percent.mt')
msc = PercentageFeatureSet(object = msc, pattern = "^MT-", col.name = 'percent.mt')
msd = PercentageFeatureSet(object = msd, pattern = "^MT-", col.name = 'percent.mt')
mse = PercentageFeatureSet(object = mse, pattern = "^MT-", col.name = 'percent.mt')
MS6 = PercentageFeatureSet(object = MS6, pattern = "^MT-", col.name = 'percent.mt')
MS7 = PercentageFeatureSet(object = MS7, pattern = "^MT-", col.name = 'percent.mt')
MS8 = PercentageFeatureSet(object = MS8, pattern = "^MT-", col.name = 'percent.mt')
MS9 = PercentageFeatureSet(object = MS9, pattern = "^MT-", col.name = 'percent.mt')
MS10 = PercentageFeatureSet(object = MS10, pattern = "^MT-", col.name = 'percent.mt')
MS24 = PercentageFeatureSet(object = MS24, pattern = "^MT-", col.name = 'percent.mt')
MS25 = PercentageFeatureSet(object = MS25, pattern = "^MT-", col.name = 'percent.mt')
MS27 = PercentageFeatureSet(object = MS27, pattern = "^MT-", col.name = 'percent.mt')
MS28 = PercentageFeatureSet(object = MS28, pattern = "^MT-", col.name = 'percent.mt')
MS29 = PercentageFeatureSet(object = MS29, pattern = "^MT-", col.name = 'percent.mt')
nmsa = PercentageFeatureSet(object = nmsa, pattern = "^MT-", col.name = 'percent.mt')
nmsb = PercentageFeatureSet(object = nmsb, pattern = "^MT-", col.name = 'percent.mt')
nmsc = PercentageFeatureSet(object = nmsc, pattern = "^MT-", col.name = 'percent.mt')
nmsd = PercentageFeatureSet(object = nmsd, pattern = "^MT-", col.name = 'percent.mt')
nmse = PercentageFeatureSet(object = nmse, pattern = "^MT-", col.name = 'percent.mt')
H2 = PercentageFeatureSet(object = H2, pattern = "^MT-", col.name = 'percent.mt')
H3 = PercentageFeatureSet(object = H3, pattern = "^MT-", col.name = 'percent.mt')
H32 = PercentageFeatureSet(object = H32, pattern = "^MT-", col.name = 'percent.mt')

###Exploration and Filtering
VlnPlot(msa, features = c('nCount_RNA', 'nFeature_RNA', 'percent.mt'), ncol = 3, pt.size = 0)
msa = subset(msa, subset = nCount_RNA < 12000 & nFeature_RNA < 3000 & percent.mt < 5)
VlnPlot(msb, features = c('nCount_RNA', 'nFeature_RNA', 'percent.mt'), ncol = 3, pt.size = 0)
msb = subset(msb, subset = nCount_RNA < 11000 & nFeature_RNA < 2500 & percent.mt < 5)
VlnPlot(msc, features = c('nCount_RNA', 'nFeature_RNA', 'percent.mt'), ncol = 3, pt.size = 0)
msc = subset(msc, subset = nCount_RNA < 12000 & nFeature_RNA < 3000 & percent.mt < 7)
VlnPlot(msd, features = c('nCount_RNA', 'nFeature_RNA', 'percent.mt'), ncol = 3, pt.size = 0)
msd = subset(msd, subset = nCount_RNA < 18000 & nFeature_RNA < 4200 & percent.mt < 7)
VlnPlot(mse, features = c('nCount_RNA', 'nFeature_RNA', 'percent.mt'), ncol = 3, pt.size = 0)
mse = subset(mse, subset = nCount_RNA < 12000 & nFeature_RNA < 2600 & percent.mt < 10)
VlnPlot(MS6, features = c('nCount_RNA', 'nFeature_RNA', 'percent.mt'),ncol = 3, pt.size = 0)
MS6 = subset(MS6, subset = nCount_RNA < 12000 & nFeature_RNA < 2500 & percent.mt < 15)
VlnPlot(MS7, features = c('nCount_RNA', 'nFeature_RNA', 'percent.mt'), ncol = 3, pt.size = 0)
MS7 = subset(MS7, subset = nCount_RNA < 9000 & nFeature_RNA < 2500 & percent.mt < 13)
VlnPlot(MS8, features = c('nCount_RNA', 'nFeature_RNA', 'percent.mt'), ncol = 3, pt.size = 0)
MS8 = subset(MS8, subset = nCount_RNA < 12000 & nFeature_RNA < 2500 & percent.mt < 12)
VlnPlot(MS9, features = c('nCount_RNA', 'nFeature_RNA', 'percent.mt'), ncol = 3, pt.size = 0)
MS9 = subset(MS9, subset = nCount_RNA < 9000 & nFeature_RNA < 2500 & percent.mt < 15)
VlnPlot(MS10, features = c('nCount_RNA', 'nFeature_RNA', 'percent.mt'), ncol = 3, pt.size = 0)
MS10 = subset(MS10, subset = nCount_RNA < 9000 & nFeature_RNA < 2500 & percent.mt < 10)
VlnPlot(MS24, features = c('nCount_RNA', 'nFeature_RNA', 'percent.mt'),ncol = 3, pt.size = 0)
MS24 = subset(MS24, subset = nCount_RNA < 7500 & nFeature_RNA < 2500 & percent.mt < 8)
VlnPlot(MS27, features = c('nCount_RNA', 'nFeature_RNA', 'percent.mt'), ncol = 3, pt.size = 0)
MS27 = subset(MS27, subset = nCount_RNA < 7500 & nFeature_RNA < 2500 & percent.mt < 15)
VlnPlot(MS28, features = c('nCount_RNA', 'nFeature_RNA', 'percent.mt'), ncol = 3, pt.size = 0)
MS28 = subset(MS28, subset = nCount_RNA < 5000 & nFeature_RNA < 2000 & percent.mt < 10)
VlnPlot(MS29, features = c('nCount_RNA', 'nFeature_RNA', 'percent.mt'), ncol = 3, pt.size = 0)
MS29 = subset(MS29, subset = nCount_RNA < 9000 & nFeature_RNA < 2300 & percent.mt < 10)
VlnPlot(H2, features = c('nCount_RNA', 'nFeature_RNA', 'percent.mt'),ncol = 3, pt.size = 0)
H2 = subset(H2, subset = nCount_RNA < 10000 & nFeature_RNA < 2300 & percent.mt < 10)
VlnPlot(H3, features = c('nCount_RNA', 'nFeature_RNA', 'percent.mt'), ncol = 3, pt.size = 0)
H3 = subset(H3, subset = nCount_RNA < 8000 & nFeature_RNA < 2500 & percent.mt < 10)
VlnPlot(H32, features = c('nCount_RNA', 'nFeature_RNA', 'percent.mt'), ncol = 3, pt.size = 0)
H32 = subset(H32, subset = nCount_RNA < 9000 & nFeature_RNA < 2500 & percent.mt < 12)
VlnPlot(nmsa, features = c('nCount_RNA', 'nFeature_RNA', 'percent.mt'), ncol = 3, pt.size = 0)
nmsa = subset(nmsa, subset = nCount_RNA < 18000 & nFeature_RNA < 3900 & percent.mt < 6)
VlnPlot(nmsb, features = c('nCount_RNA', 'nFeature_RNA', 'percent.mt'), ncol = 3, pt.size = 0)
nmsb = subset(nmsb, subset = nCount_RNA < 14000 & nFeature_RNA < 3500 & percent.mt < 7)
VlnPlot(nmsc, features = c('nCount_RNA', 'nFeature_RNA', 'percent.mt'), ncol = 3, pt.size = 0)
nmsc = subset(nmsc, subset = nCount_RNA < 15000 & nFeature_RNA < 3000 & percent.mt < 10)
VlnPlot(nmsd, features = c('nCount_RNA', 'nFeature_RNA', 'percent.mt'), ncol = 3, pt.size = 0)
nmsd = subset(nmsa, subset = nCount_RNA < 13000 & nFeature_RNA < 2500 & percent.mt < 7)
VlnPlot(nmse, features = c('nCount_RNA', 'nFeature_RNA', 'percent.mt'),ncol = 3, pt.size = 0)
nmse = subset(nmse, subset = nCount_RNA < 10000 & nFeature_RNA < 2500 & percent.mt < 7)

####Merge the datasets for MS and NON-MS(NMS)
merged_int = list(msa = msa, msb = msb, msc = msc, msd = msd , mse = mse, MS6 = MS6, 
                  MS7 = MS7, MS8 = MS8, MS9 = MS9, MS10 = MS10, MS24 = MS24, 
                  MS25 = MS25, MS27 = MS27, MS28 = MS28, MS29 = MS29, 
                  nmsa = nmsa, nmsb = nmsb, nmsc = nmsc, nmsd = nmsd, nmse = nmse, 
                  H2 = H2, H3 = H3, H32 = H32)

#remove doublets
merged_int = lapply(X = merged_int, FUN = function(X){
  X = NormalizeData(X)
  X = FindVariableFeatures(X, selection.method = "vst", nFeatures = 2000)
  X = ScaleData(X)
  X = RunPCA(X)
  X = FindNeighbors(X, dims = 1:20)
  X = FindClusters(X)
  X = RunUMAP(X, dims = 1:20)
})

set.seed(42)

for (i in 1:23) {
  sce <- as.SingleCellExperiment(merged_int[[i]], graphs=c("pca","umap"))
  sce <- scDblFinder(sce)
  merged_int[[i]]$db_score <- sce$scDblFinder.score
  merged_int[[i]]$db_type <- factor( sce$scDblFinder.class,
  levels=c("singlet", "doublet") )
}

merged = unlist(merged_int)
msa = subset(merged$msa, subset = db_type == 'singlet')
msb = subset(merged$msb, subset = db_type == 'singlet')
msc = subset(merged$msc, subset = db_type == 'singlet')
msd = subset(merged$msd, subset = db_type == 'singlet')
mse = subset(merged$mse, subset = db_type == 'singlet')
MS6 = subset(merged$MS6, subset = db_type == 'singlet')
MS7 = subset(merged$MS7, subset = db_type == 'singlet')
MS8 = subset(merged$MS8, subset = db_type == 'singlet')
MS9 = subset(merged$MS9, subset = db_type == 'singlet')
MS10 = subset(merged$MS10, subset = db_type == 'singlet')
MS24 = subset(merged$MS24, subset = db_type == 'singlet')
MS25 = subset(merged$MS25, subset = db_type == 'singlet')
MS27 = subset(merged$MS27, subset = db_type == 'singlet')
MS28 = subset(merged$MS28, subset = db_type == 'singlet')
MS29 = subset(merged$MS29, subset = db_type == 'singlet')
nmsa = subset(merged$nmsa, subset = db_type == 'singlet')
nmsb = subset(merged$nmsb, subset = db_type == 'singlet')
nmsc = subset(merged$nmsc, subset = db_type == 'singlet')
nmsd = subset(merged$nmsd, subset = db_type == 'singlet')
nmse = subset(merged$nmse, subset = db_type == 'singlet')
H2 = subset(merged$H2, subset = db_type == 'singlet')
H3 = subset(merged$H3, subset = db_type == 'singlet')
H32 = subset(merged$H32, subset = db_type == 'singlet')


#####mMerge for HARMONY or MS and NON-MS(NMS)
inf_int = merge(x = msa, y=c(msb, msc, msd, mse, nmsa, nmsb,H2, H3, H32, MS6, 
MS7, MS8, MS9, MS10, nmsc, nmsd, nmse, MS24, MS25, MS27, MS28, MS29))

####More exploration of merged dataset for optimal QC
Idents(inf_int) = 'orig.ident'
plot1 <- FeatureScatter(inf_int, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(inf_int, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Standard workflow for merged datasets
inf_int = inf_int |> NormalizeData() |> 
  FindVariableFeatures(nfeatures = 3000) |>
  ScaleData() |> 
  RunPCA(npcs = 50)
inf_int = RunUMAP(inf_int, reduction = "pca", dims = 1:30, verbose = F)
a = DimPlot(inf_int,reduction = "umap", group.by = "orig.ident") + 
plot_annotation(title =  "before integration")
inf_int$patient = c('a', "b", "c", "d", "e", '6', "7", "8", "9", "10", '24', '25', '27',
                    '28', '29', 'sa', 'sb','sc', 'sd', 'se', '2', '3', '32')
View(inf_int@meta.data)
plot_px = DimPlot(inf_int,reduction = "umap", group.by = "patient") + 
  plot_annotation(title =  "before integration") + NoLegend()
plot_type = DimPlot(inf_int,reduction = "umap", group.by = "orig.ident") + 
  plot_annotation(title =  "before integration") + NoLegend()
plot_px + plot_type

####integration with Harmony
inf_int<- RunHarmony(inf_int, group.by.vars = "orig.ident", dims.use = 1:20, 
      max.iter.harmony = 50, plot_convergence = T)

###check embeddings
harmony_embeddings <- Embeddings(inf_int, 'harmony')
harmony_embeddings[1:5, 1:5]

###Visualize integration
inf_int <- inf_int  %>% 
  RunUMAP(reduction = "harmony", dims = 1:20, verbose = F) %>% 
  FindNeighbors(reduction = "harmony", k.param = 10, dims = 1:20) %>% 
  FindClusters() %>% 
  identity()
inf_int =  SetIdent(inf_int, value = "seurat_clusters")
e = DimPlot(inf_int, group.by = "orig.ident", reduction = 'harmony') + 
  plot_annotation(title = "after integration (Harmony)") + 
  NoLegend()
f = DimPlot(inf_int, group.by = "patient", reduction = 'harmony') + 
  plot_annotation(title = "after integration (Harmony)") + 
  NoLegend()
c = DimPlot(inf_int, label = T, reduction = "harmony") + 
  plot_annotation(title = "cluster identity (Harmony)") + 
  NoLegend()
d = DimPlot(inf_int, label = T, reduction = "harmony", split.by = 'orig.ident') + 
  plot_annotation(title = "cluster identity (Harmony)") + 
  NoLegend()
e + f
c + d
saveRDS(inf_int, file = "tot_harmony_csf.Rds")

########PERIPHERAL BLOOD MONONUCLEAR CELLS

####load in files 
fin3 = list.files(path = "/Users/tee/Downloads/scr/GSE133028/pbmc/5/H/", full.names = F, 
                  recursive = F)
for (i in fin3){
  names = gsub(pattern = "_t[0-9].csv", replacement = "", i)
  print(names)
  cts = read.csv(paste0("/Users/tee/Downloads/scr/GSE133028/pbmc/5/H/", i))
  assign(names, CreateSeuratObject(counts = cts, min.cells = 3, min.features = 200, 
                                   project = "NMS"))
}

fin4 = list.files(path = "/Users/tee/Downloads/scr/GSE133028/pbmc/5/M/", full.names = F,
                  recursive = F)
for (i in fin4){
  names = gsub(pattern = "_t[0-9].csv", replacement = "", i)
  print(names)
  cts = read.csv(paste0("/Users/tee/Downloads/scr/GSE133028/pbmc/5/M/", i))
  assign(names, CreateSeuratObject(counts = cts, min.cells = 3, min.features = 200, 
                                   project = "MS"))
}

msa = Read10X(data.dir = '/Users/tee/ms/dat/PBMC/MS-MA')
msb = Read10X(data.dir = '/Users/tee/ms/dat/PBMC/MS-MB')
msc = Read10X(data.dir = '/Users/tee/ms/dat/PBMC/MS-MC')
msd = Read10X(data.dir = '/Users/tee/ms/dat/PBMC/MS-MD')
mse = Read10X(data.dir = '/Users/tee/ms/dat/PBMC/MS-ME')
nmsa = Read10X(data.dir = '/Users/tee/ms/dat/PBMC/NMS-NMSA1')
nmsb = Read10X(data.dir = '/Users/tee/ms/dat/PBMC/NMS-NMSB2')
nmsc = Read10X(data.dir = '/Users/tee/ms/dat/PBMC/NMS-NMSC1')
nmsd = Read10X(data.dir = '/Users/tee/ms/dat/PBMC/NMS-NMSD2')
nmse = Read10X(data.dir = '/Users/tee/ms/dat/PBMC/NMS-NMSE2')
msa = CreateSeuratObject(counts=msa, min.cells = 3, min.features = 200, project = 'MS')
msb = CreateSeuratObject(counts=msb, min.cells = 3, min.features = 200, project = 'MS')
msc = CreateSeuratObject(counts=msc, min.cells = 3, min.features = 200, project = 'MS')
msd = CreateSeuratObject(counts=msd, min.cells = 3, min.features = 200, project = 'MS')
mse = CreateSeuratObject(counts=mse, min.cells = 3, min.features = 200, project = 'MS')
nmsa = CreateSeuratObject(counts=nmsa, min.cells = 3, min.features = 200, project = 'NMS')
nmsb = CreateSeuratObject(counts=nmsb, min.cells = 3, min.features = 200, project = 'NMS')
nmsc = CreateSeuratObject(counts=nmsc, min.cells = 3, min.features = 200, project = 'NMS')
nmsd = CreateSeuratObject(counts=nmsd, min.cells = 3, min.features = 200, project = 'NMS')
nmse = CreateSeuratObject(counts=nmse, min.cells = 3, min.features = 200, project = 'NMS')

#####QC
## Mitochondrion Percentage
msa = PercentageFeatureSet(object = msa, pattern = "^MT-", col.name = 'percent.mt')
msb = PercentageFeatureSet(object = msb, pattern = "^MT-", col.name = 'percent.mt')
msc = PercentageFeatureSet(object = msc, pattern = "^MT-", col.name = 'percent.mt')
msd = PercentageFeatureSet(object = msd, pattern = "^MT-", col.name = 'percent.mt')
mse = PercentageFeatureSet(object = mse, pattern = "^MT-", col.name = 'percent.mt')
MS6 = PercentageFeatureSet(object = MS6, pattern = "^MT-", col.name = 'percent.mt')
MS7 = PercentageFeatureSet(object = MS7, pattern = "^MT-", col.name = 'percent.mt')
MS8 = PercentageFeatureSet(object = MS8, pattern = "^MT-", col.name = 'percent.mt')
MS9 = PercentageFeatureSet(object = MS9, pattern = "^MT-", col.name = 'percent.mt')
MS10 = PercentageFeatureSet(object = MS10, pattern = "^MT-", col.name = 'percent.mt')
MS24 = PercentageFeatureSet(object = MS24, pattern = "^MT-", col.name = 'percent.mt')
MS25 = PercentageFeatureSet(object = MS25, pattern = "^MT-", col.name = 'percent.mt')
MS27 = PercentageFeatureSet(object = MS27, pattern = "^MT-", col.name = 'percent.mt')
MS28 = PercentageFeatureSet(object = MS28, pattern = "^MT-", col.name = 'percent.mt')
MS29 = PercentageFeatureSet(object = MS29, pattern = "^MT-", col.name = 'percent.mt')
nmsa = PercentageFeatureSet(object = nmsa, pattern = "^MT-", col.name = 'percent.mt')
nmsb = PercentageFeatureSet(object = nmsb, pattern = "^MT-", col.name = 'percent.mt')
nmsc = PercentageFeatureSet(object = nmsc, pattern = "^MT-", col.name = 'percent.mt')
nmsd = PercentageFeatureSet(object = nmsd, pattern = "^MT-", col.name = 'percent.mt')
nmse = PercentageFeatureSet(object = nmse, pattern = "^MT-", col.name = 'percent.mt')
H2 = PercentageFeatureSet(object = H2, pattern = "^MT-", col.name = 'percent.mt')
H3 = PercentageFeatureSet(object = H3, pattern = "^MT-", col.name = 'percent.mt')
H32 = PercentageFeatureSet(object = H32, pattern = "^MT-", col.name = 'percent.mt')

###Exploration and Filtering
VlnPlot(msa, features = c('nCount_RNA', 'nFeature_RNA', 'percent.mt'), ncol = 3, pt.size = 0)
msa = subset(msa, subset = nCount_RNA < 10000 & nFeature_RNA < 2000 & percent.mt < 7)
VlnPlot(msb, features = c('nCount_RNA', 'nFeature_RNA', 'percent.mt'), ncol = 3, pt.size = 0)
msb = subset(msb, subset = nCount_RNA < 11000 & nFeature_RNA < 2500 & percent.mt < 10)
VlnPlot(msc, features = c('nCount_RNA', 'nFeature_RNA', 'percent.mt'), ncol = 3, pt.size = 0)
msc = subset(msc, subset = nCount_RNA < 10000 & nFeature_RNA < 2500 & percent.mt < 9)
VlnPlot(msd, features = c('nCount_RNA', 'nFeature_RNA', 'percent.mt'), ncol = 3, pt.size = 0)
msd = subset(msd, subset = nCount_RNA < 7000 & nFeature_RNA < 2000 & percent.mt < 8)
VlnPlot(mse, features = c('nCount_RNA', 'nFeature_RNA', 'percent.mt'), ncol = 3, pt.size = 0)
mse = subset(mse, subset = nCount_RNA < 12000 & nFeature_RNA < 2600 & percent.mt < 10)
VlnPlot(MS6, features = c('nCount_RNA', 'nFeature_RNA', 'percent.mt'),ncol = 3, pt.size = 0)
MS6 = subset(MS6, subset = nCount_RNA < 10000 & nFeature_RNA < 2500 & percent.mt < 12)
VlnPlot(MS7, features = c('nCount_RNA', 'nFeature_RNA', 'percent.mt'), ncol = 3, pt.size = 0)
MS7 = subset(MS7, subset = nCount_RNA < 10000 & nFeature_RNA < 2500 & percent.mt < 15)
VlnPlot(MS8, features = c('nCount_RNA', 'nFeature_RNA', 'percent.mt'), ncol = 3, pt.size = 0)
MS8 = subset(MS8, subset = nCount_RNA < 9000 & nFeature_RNA < 2500 & percent.mt < 11)
VlnPlot(MS9, features = c('nCount_RNA', 'nFeature_RNA', 'percent.mt'), ncol = 3, pt.size = 0)
MS9 = subset(MS9, subset = nCount_RNA < 7000 & nFeature_RNA < 2500 & percent.mt < 15)
VlnPlot(MS10, features = c('nCount_RNA', 'nFeature_RNA', 'percent.mt'), ncol = 3, pt.size = 0)
MS10 = subset(MS10, subset = nCount_RNA < 9000 & nFeature_RNA < 2500 & percent.mt < 10)
VlnPlot(MS24, features = c('nCount_RNA', 'nFeature_RNA', 'percent.mt'),ncol = 3, pt.size = 0)
MS24 = subset(MS24, subset = nCount_RNA < 9500 & nFeature_RNA < 2500 & percent.mt < 10)
VlnPlot(MS27, features = c('nCount_RNA', 'nFeature_RNA', 'percent.mt'), ncol = 3, pt.size = 0)
MS27 = subset(MS27, subset = nCount_RNA < 10000 & nFeature_RNA < 2500 & percent.mt < 15)
VlnPlot(MS28, features = c('nCount_RNA', 'nFeature_RNA', 'percent.mt'), ncol = 3, pt.size = 0)
MS28 = subset(MS28, subset = nCount_RNA < 5000 & nFeature_RNA < 2000 & percent.mt < 12)
VlnPlot(MS29, features = c('nCount_RNA', 'nFeature_RNA', 'percent.mt'), ncol = 3, pt.size = 0)
MS29 = subset(MS29, subset = nCount_RNA < 9000 & nFeature_RNA < 2300 & percent.mt < 10)
VlnPlot(H2, features = c('nCount_RNA', 'nFeature_RNA', 'percent.mt'),ncol = 3, pt.size = 0)
H2 = subset(H2, subset = nCount_RNA < 10000 & nFeature_RNA < 2300 & percent.mt < 12)
VlnPlot(H3, features = c('nCount_RNA', 'nFeature_RNA', 'percent.mt'), ncol = 3, pt.size = 0)
H3 = subset(H3, subset = nCount_RNA < 10000 & nFeature_RNA < 2500 & percent.mt < 10)
VlnPlot(H32, features = c('nCount_RNA', 'nFeature_RNA', 'percent.mt'), ncol = 3, pt.size = 0)
H32 = subset(H32, subset = nCount_RNA < 9500 & nFeature_RNA < 2500 & percent.mt < 14)
VlnPlot(nmsa, features = c('nCount_RNA', 'nFeature_RNA', 'percent.mt'), ncol = 3, pt.size = 0)
nmsa = subset(nmsa, subset = nCount_RNA < 12000 & nFeature_RNA < 2700 & percent.mt < 7)
VlnPlot(nmsb, features = c('nCount_RNA', 'nFeature_RNA', 'percent.mt'), ncol = 3, pt.size = 0)
nmsb = subset(nmsb, subset = nCount_RNA < 8000 & nFeature_RNA < 1500 & percent.mt < 7.5)
VlnPlot(nmsc, features = c('nCount_RNA', 'nFeature_RNA', 'percent.mt'), ncol = 3, pt.size = 0)
nmsc = subset(nmsc, subset = nCount_RNA < 8000 & nFeature_RNA < 1600 & percent.mt < 10)
VlnPlot(nmsd, features = c('nCount_RNA', 'nFeature_RNA', 'percent.mt'), ncol = 3, pt.size = 0)
nmsd = subset(nmsa, subset = nCount_RNA < 11000 & nFeature_RNA < 2000 & percent.mt < 9)
VlnPlot(nmse, features = c('nCount_RNA', 'nFeature_RNA', 'percent.mt'),ncol = 3, pt.size = 0)
nmse = subset(nmse, subset = nCount_RNA < 8000 & nFeature_RNA < 1800 & percent.mt < 8)

####Merge the datasets for MS and NON-MS(NMS)
merged_int = list(msa = msa, msb = msb, msc = msc, msd = msd , mse = mse, MS6 = MS6, 
                  MS7 = MS7, MS8 = MS8, MS9 = MS9, MS10 = MS10, MS24 = MS24, 
                  MS25 = MS25, MS27 = MS27, MS28 = MS28, MS29 = MS29, 
                  nmsa = nmsa, nmsb = nmsb, nmsc = nmsc, nmsd = nmsd, nmse = nmse, 
                  H2 = H2, H3 = H3, H32 = H32)

#remove doublets
merged_int = lapply(X = merged_int, FUN = function(X){
  X = NormalizeData(X)
  X = FindVariableFeatures(X, selection.method = "vst", nFeatures = 2000)
  X = ScaleData(X)
  X = RunPCA(X)
  X = FindNeighbors(X, dims = 1:20)
  X = FindClusters(X)
  X = RunUMAP(X, dims = 1:20)
})

set.seed(42)

for (i in 1:23) {
  sce <- as.SingleCellExperiment(merged_int[[i]], graphs=c("pca","umap"))
  sce <- scDblFinder(sce)
  merged_int[[i]]$db_score <- sce$scDblFinder.score
  merged_int[[i]]$db_type <- factor( sce$scDblFinder.class,
                                     levels=c("singlet", "doublet") )
}

merged = unlist(merged_int)
msa = subset(merged$msa, subset = db_type == 'singlet')
msb = subset(merged$msb, subset = db_type == 'singlet')
msc = subset(merged$msc, subset = db_type == 'singlet')
msd = subset(merged$msd, subset = db_type == 'singlet')
mse = subset(merged$mse, subset = db_type == 'singlet')
MS6 = subset(merged$MS6, subset = db_type == 'singlet')
MS7 = subset(merged$MS7, subset = db_type == 'singlet')
MS8 = subset(merged$MS8, subset = db_type == 'singlet')
MS9 = subset(merged$MS9, subset = db_type == 'singlet')
MS10 = subset(merged$MS10, subset = db_type == 'singlet')
MS24 = subset(merged$MS24, subset = db_type == 'singlet')
MS25 = subset(merged$MS25, subset = db_type == 'singlet')
MS27 = subset(merged$MS27, subset = db_type == 'singlet')
MS28 = subset(merged$MS28, subset = db_type == 'singlet')
MS29 = subset(merged$MS29, subset = db_type == 'singlet')
nmsa = subset(merged$nmsa, subset = db_type == 'singlet')
nmsb = subset(merged$nmsb, subset = db_type == 'singlet')
nmsc = subset(merged$nmsc, subset = db_type == 'singlet')
nmsd = subset(merged$nmsd, subset = db_type == 'singlet')
nmse = subset(merged$nmse, subset = db_type == 'singlet')
H2 = subset(merged$H2, subset = db_type == 'singlet')
H3 = subset(merged$H3, subset = db_type == 'singlet')
H32 = subset(merged$H32, subset = db_type == 'singlet')


#####mMerge for HARMONY or MS and NON-MS(NMS)
inf_int = merge(x = msa, y=c(msb, msc, msd, mse, nmsa, nmsb,H2, H3, H32, MS6, 
  MS7, MS8, MS9, MS10, nmsc, nmsd, nmse, MS24, MS25, MS27, MS28, MS29))

####More exploration of merged dataset for optimal QC
Idents(inf_int) = 'orig.ident'
plot1 <- FeatureScatter(inf_int, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(inf_int, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Standard workflow for merged datasets
inf_int = inf_int |> NormalizeData() |> 
  FindVariableFeatures(nfeatures = 3000) |>
  ScaleData() |> 
  RunPCA(npcs = 50)
inf_int = RunUMAP(inf_int, reduction = "pca", dims = 1:30, verbose = F)
a = DimPlot(inf_int,reduction = "umap", group.by = "orig.ident", raster = F) + 
  plot_annotation(title =  "before integration")
inf_int$patient = c('a', "b", "c", "d", "e", '6', "7", "8", "9", "10", '24', '25', '27',
  '28', '29', 'sa', 'sb','sc', 'sd', 'se', '2', '3', '32')
View(inf_int@meta.data)
plot_px = DimPlot(inf_int,reduction = "umap", group.by = "patient", raster = F) + 
  plot_annotation(title =  "before integration") + NoLegend()
plot_type = DimPlot(inf_int,reduction = "umap", group.by = "orig.ident", raster = F) + 
  plot_annotation(title =  "before integration") + NoLegend()
plot_px + plot_type

####integration with Harmony
inf_int<- RunHarmony(inf_int, group.by.vars = "orig.ident", dims.use = 1:20, 
                     max.iter.harmony = 50, plot_convergence = T)

###check embeddings
harmony_embeddings <- Embeddings(inf_int, 'harmony')
harmony_embeddings[1:5, 1:5]

###Visualize integration
inf_int <- inf_int  %>% 
  RunUMAP(reduction = "harmony", dims = 1:20, verbose = F) %>% 
  FindNeighbors(reduction = "harmony", k.param = 10, dims = 1:20) %>% 
  FindClusters() %>% 
  identity()
inf_int =  SetIdent(inf_int, value = "seurat_clusters")
e = DimPlot(inf_int, group.by = "orig.ident", reduction = 'harmony', raster = F) + 
  plot_annotation(title = "after integration (Harmony)") + 
  NoLegend()
f = DimPlot(inf_int, group.by = "patient", reduction = 'harmony', raster = F) + 
  plot_annotation(title = "after integration (Harmony)") + 
  NoLegend()
c = DimPlot(inf_int, label = T, reduction = "harmony", raster = F) + 
  plot_annotation(title = "cluster identity (Harmony)") + 
  NoLegend()
d = DimPlot(inf_int, label = T, reduction = "harmony", split.by = 'orig.ident', raster = F) + 
  plot_annotation(title = "cluster identity (Harmony)") + 
  NoLegend()
e + f
c + d
saveRDS(inf_int, file = "tot_harmony_pbmc.Rds")


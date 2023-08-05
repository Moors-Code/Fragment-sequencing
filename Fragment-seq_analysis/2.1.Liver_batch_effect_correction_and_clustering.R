########## Part 2.1: Batch effect correction and clustering of liver samples ##########
#This part merges all liver samples by applying batch effect correction following the vignette of mnn-correction 
#http://bioconductor.org/books/3.14/OSCA.multisample/integrating-datasets.html#mnn-correction
#after batch correction all samples undergo removal of low quality cells, normalization, scaling and clustering

########## Prepare environment ##########
###Setting the working directory 
setwd("/mnt/khandler/R_projects/Fragment-sequencing/Fragment-seq_analysis/")

###Load packages and functions 
source("./functions_and_packages/1.Packages.R")

########## Load Seurat objects after MULTI-seq barcode integration ##########
S1 <- readRDS(file = "./data_files_generated/SpS_LiverMets_sample1.Rda")
S2 <- readRDS(file = "./data_files_generated/SpS_LiverMets_sample2.Rda")
S3 <- readRDS(file = "./data_files_generated/SpS_LiverMets_sample3.Rda")
S4 <- readRDS(file = "./data_files_generated/SpS_LiverMets_sample4.Rda")
S5 <- readRDS(file = "./data_files_generated/SpS_LiverMets_sample5.Rda") 
S6 <- readRDS(file = "./data_files_generated/SpS_LiverMets_sample6.Rda")
S7 <- readRDS(file = "./data_files_generated/SpS_LiverMets_sample7.Rda")
S8 <- readRDS(file = "./data_files_generated/SpS_LiverMets_sample8.Rda")
S9 <- readRDS(file = "./data_files_generated/SpS_LiverMets_sample9.Rda")
S10 <- readRDS(file = "./data_files_generated/SpS_LiverWT_sample10.Rda")

#add prefix to the fragment ID 
S1$fragment <- paste(S1$fragment, "_1", sep="")
S2$fragment <- paste(S2$fragment, "_2", sep="")
S3$fragment <- paste(S3$fragment, "_3", sep="")
S4$fragment <- paste(S4$fragment, "_4", sep="")
S5$fragment <- paste(S5$fragment, "_5", sep="")
S6$fragment <- paste(S6$fragment, "_6", sep="")
S7$fragment <- paste(S7$fragment, "_7", sep="")
S8$fragment <- paste(S8$fragment, "_8", sep="")
S9$fragment <- paste(S9$fragment, "_9", sep="")
S10$fragment <- paste(S10$fragment, "_10", sep="")

########## Merging and Clustering of samples before batch effect correction ##########
#merge Seurat objects
merged <- merge(S1, y = c(S2,S3,S4,S5,S6,S7,S8,S9,S10), add.cell.ids = c("1","2","3","4","5","6","7","8","9","10"))
#add percentage of mitochondrial transcripts to Seurat object meta data 
mito.features <- grep(pattern = "^mt-", x = rownames(x = merged), value = TRUE)
percent.mito <- Matrix::colSums(x = GetAssayData(object = merged, slot = 'counts')[mito.features, ]) / Matrix::colSums(x = GetAssayData(object = merged, slot = 'counts'))
merged <- AddMetaData(object = merged, metadata = percent.mito, col.name = "percent.mito")
#normalize and scale data 
merged <- SCTransform(merged, vars.to.regress = c("nCount_RNA","percent.mito"), verbose = FALSE)
merged <- RunPCA(object = merged, features = VariableFeatures(object = merged), npcs = 30, verbose = FALSE)
#plot PCA space before batch effect correction --> you can see a batch effect 
DimPlot(merged, reduction = "pca", group.by = "orig.ident")+ 
  ggtitle("PCA space before batch correction")+ ggsave("./figures/2.1/pca_space_before_batch_correction.pdf")

########## Batch effect correction #########
#convert to single cell experiment object (SCE) 
merged_diet <- DietSeurat(merged)
merged_sce <- as.SingleCellExperiment(merged_diet)
#extract variable genes per sample and put them in a list 
batch.vars <- lapply(unique(merged_sce$orig.ident), function(BATCH) {
  batch.var <- modelGeneVar(merged_sce[,merged_sce$orig.ident==BATCH])})
#combine variable genes 
comb.var <- combineVar(batch.vars)
#take top 6000 genes 
hvgs <- getTopHVGs(comb.var,n=6000)
set.seed(42)
#run batch correction 
mnn.out <- fastMNN(merged_sce,
                   batch=merged_sce$orig.ident,
                   subset.row=hvgs,
                   d=50,
                   correct.all=TRUE)
stopifnot(identical(colnames(merged_sce),colnames(mnn.out)))
join.mnn.pca <- reducedDim(mnn.out,"corrected")
rownames(join.mnn.pca) <- colnames(mnn.out)
reducedDim(merged_sce,"PCA.MNN") <- join.mnn.pca
assays(merged_sce)[["reconstructed"]] <- assays(mnn.out)[["reconstructed"]]
set.seed(42)

#convert back to Seurat object
merged_S <- as.Seurat(merged_sce)
#plot PCA space after batch effect correction --> batch effect is gone 
DimPlot(merged_S, reduction = "PCA.MNN", group.by = "orig.ident") + 
  ggtitle("PCA space after batch correction")+ ggsave("./figures/2.1/pca_space_after_batch_correction.pdf")

########## Quality control and cutoff ##########
###check QC measurements 
p1 <- VlnPlot(merged_S, feature = "percent.mito")
p2 <- VlnPlot(merged_S, feature = "nCount_RNA")
p3 <- VlnPlot(merged_S, feature = "nFeature_RNA")
p4 <- FeatureScatter(merged_S, feature1 = "nCount_RNA", feature2 = "percent.mito")
p5 <- FeatureScatter(merged_S, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
median_percent.mito <- median(merged_S@meta.data$percent.mito)
median_nCount_RNA <- median(merged_S@meta.data$nFeature_RNA)
median_nFeature_RNA <- median(merged_S@meta.data$nFeature_RNA)
plot_grid(p1, p2,p3,p4,p5, labels = c(paste("median=",median_percent.mito, sep = ""), paste("median=",median_nCount_RNA, sep = ""),
                                      paste("median=",median_nFeature_RNA, sep = ""))) + ggtitle("Quality measurements")+
  ggsave("./figures/2.1/Quality_measures.pdf", width = 15, height = 10)

###Remove low quality cells 
#cells with higher percentage of mitochondrial genes than 0.2
#cells with higher than 7500 features 
#cells with lower than 200 features 
merged_S <- subset(merged_S, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mito < 0.2 )

########## Clustering of data ##########
merged_S <- subset(merged_S, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mito < 0.2 )
merged_S <- FindNeighbors(object = merged_S, dims = 1:10, reduction = "PCA.MNN")
merged_cl <- FindClusters(merged_S, resolution = 1, random.seed = 5, algorithm = 1, graph.name = "RNA_nn")
merged_cl <- RunUMAP(merged_cl, dims = 1:10, seed.use = 5, reduction = "PCA.MNN")
DimPlot(merged_cl, reduction = "umap", label = TRUE, pt.size = .5, group.by = "seurat_clusters")

########## Save Seurat object after batch correction, QC and clustering ##########
saveRDS(merged_cl, file = "./data_files_generated/LiverMerged_afterBC.Rda")






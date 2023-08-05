########## Part 8: Liver comp scRNAseq ##########
#This part compares fragment-seq data with conventional scRNA-seq
#For comparing cell types all objects will be annotated together 
#For comparison of UMI counts, gene counts and ratio mito/cytopl genes zUMI outputs were used with 50000 reads/cell downsampling 

########## Prepare environment ##########
###Setting the working directory 
setwd("/mnt/khandler/R_projects/Fragment-sequencing/Fragment-seq_analysis/")

###Load packages and functions 
source("./functions_and_packages/1.Packages.R")
source("./functions_and_packages/2.Functions_preprocessing.R")

########## Annotate zUMI outputs ##########
###load human and mouse ensemble ID/Gene ID dataframes 
##mouse gene/ensemble IDs 
ensembl<-useEnsembl(biomart="ensembl")
list<-listDatasets(ensembl)
mart <- useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl",version = 95)
attributes<-listAttributes(mart)
gene_ids <- getBM(attributes = c("ensembl_gene_id_version","external_gene_name"), mart = mart)
write.csv(gene_ids, "./data/gene_ids_ensemble_ids_mouse.csv")

##human gene/ensemble IDs 
mart1 <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl",version = 95)
attributes1<-listAttributes(mart1)
gene_ids1 <- getBM(attributes = c("ensembl_gene_id_version","external_gene_name"), mart = mart1)
write.csv(gene_ids1, "./data/gene_ids_ensemble_ids_human.csv")

gene_ids_human <- read.csv( "./data/gene_ids_ensemble_ids_human.csv")
gene_ids_mouse <- read.csv( "./data/gene_ids_ensemble_ids_mouse.csv")

#####annotate the zUMI output 
dge <- readRDS("/home/khandler/NAS/Kristina/Mets_sc1.dgecounts.rds")
Msc1 <- Annotation_mouse(dge)
Msc1_S <- CreateSeuratObject(Msc1, "sc1",min.cells = 3, min.features = 200)

dge <- readRDS("/home/khandler/NAS/Kristina/Mets_sc2.dgecounts.rds")
Msc2 <- Annotation_mouse(dge)
Msc2_S <- CreateSeuratObject(Msc2, "sc2",min.cells = 3, min.features = 200)

###save annotated objects 
saveRDS(Msc1_S, file = "./data_files_generated/scRNAseq_sample1.Rda")
saveRDS(Msc2_S, file = "./data_files_generated/scRNAseq_sample2.Rda")

########## Load Seurat objects merge and annotate together ##########
#only consider CRC injected samples 
S1 <- readRDS(file = "./data_files_generated/SpS_LiverMets_sample1.Rda")
S2 <- readRDS(file = "./data_files_generated/SpS_LiverMets_sample2.Rda")
S3 <- readRDS(file = "./data_files_generated/SpS_LiverMets_sample3.Rda")
S4 <- readRDS(file = "./data_files_generated/SpS_LiverMets_sample4.Rda")
S5 <- readRDS(file = "./data_files_generated/SpS_LiverMets_sample5.Rda") 
S6 <- readRDS(file = "./data_files_generated/SpS_LiverMets_sample6.Rda")
S7 <- readRDS(file = "./data_files_generated/SpS_LiverMets_sample7.Rda")
S8 <- readRDS(file = "./data_files_generated/SpS_LiverMets_sample8.Rda")
S9 <- readRDS(file = "./data_files_generated/SpS_LiverMets_sample9.Rda")
S11 <- readRDS(file = "./data_files_generated/scRNAseq_sample1.Rda")
S12 <- readRDS(file = "./data_files_generated/scRNAseq_sample2.Rda")

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

#####Merging and Clustering of samples before batch effect correction 
#merge Seurat objects
merged <- merge(S1, y = c(S2,S3,S4,S5,S6,S7,S8,S9,S11,S12), add.cell.ids = c("1","2","3","4","5","6","7","8","9","11","12"))
#add percentage of mitochondrial transcripts to Seurat object meta data 
mito.features <- grep(pattern = "^mt-", x = rownames(x = merged), value = TRUE)
percent.mito <- Matrix::colSums(x = GetAssayData(object = merged, slot = 'counts')[mito.features, ]) / Matrix::colSums(x = GetAssayData(object = merged, slot = 'counts'))
merged <- AddMetaData(object = merged, metadata = percent.mito, col.name = "percent.mito")
#normalize and scale data 
merged <- SCTransform(merged, vars.to.regress = c("nCount_RNA","percent.mito"), verbose = FALSE)
merged <- RunPCA(object = merged, features = VariableFeatures(object = merged), npcs = 30, verbose = FALSE)
#plot PCA space before batch effect correction --> you can see a batch effect 
DimPlot(merged, reduction = "pca", group.by = "orig.ident")

#####Batch effect correction 
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
DimPlot(merged_S, reduction = "PCA.MNN", group.by = "orig.ident") 

#####Remove low quality cells 
#cells with higher percentage of mitochondrial genes than 0.2
#cells with higher than 7500 features 
#cells with lower than 200 features 
merged_S <- subset(merged_S, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mito < 0.2 )

#####Clustering of data 
merged_S <- subset(merged_S, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mito < 0.2 )
merged_S <- FindNeighbors(object = merged_S, dims = 1:10, reduction = "PCA.MNN")
merged_cl <- FindClusters(merged_S, resolution = 1, random.seed = 5, algorithm = 1, graph.name = "RNA_nn")
merged_cl <- RunUMAP(merged_cl, dims = 1:10, seed.use = 5, reduction = "PCA.MNN")
DimPlot(merged_cl, reduction = "umap", label = TRUE, pt.size = .5, group.by = "seurat_clusters")

#####add protocol ID 
current.cluster.ids <- c("4M1","4M2","6M1","6M2","6M3","M1","M3","M4","M5","sc1","sc2")
new.cluster.ids <- c("SpS","SpS","SpS","SpS","SpS","SpS","SpS","SpS","SpS","sc","sc")
merged_cl$protocol <- plyr::mapvalues(x = merged_cl$orig.ident, from = current.cluster.ids, to = new.cluster.ids)

#####Broad annotation 
Idents(merged_cl) <- "seurat_clusters"
p <- DotPlot(merged_cl, features = c(
  "Ighm", "Igkc", #B cells: clusters 6,9,26
  "Trac","Cd3d" ,#T cells: clusters 8,12
  "Klra8","Cma1" ,#NK cells: clusters 8
  "Clec4f","Vsig4", #Kupffer cells: clusters 0,1,7,10,17,23
  "Bcl11a", "Ccr9", #DCs: cluster 21
  "S100a4", "Itgax", #general Monocytes: clusters 11,14,15,
  "Hba-a2", #Red blood cells: cluster 25
  "Csf3r", "S100a8",#Neutrophils: cluster 19
  "Siglech","Cox6a2", #pDCs: cluster 21
  "Fcer1a","Cpa3",#Basophils: cluster 16
  "Ttr", "Alb", "Cyp2e1", #Hepatocytes: cluster 20
  "Pecam1","Dll4", #LSECs: clusters 2,3,4,5,13,18,24
  "Lrat","Reln", #Stellate cells: cluster 22
  "Carmn","Nr1h5",  #Stromal cells: cluster 22
  "Svep1","Ncam1", #Fibroblasts: cluster 22
  "Krt19", "Epcam","Sox9", #Cholangiocytes: cluster  16
  "Sprr2a2","Pglyrp1","Gpx2" #Metastatic cells: cluster 16
)) 
p + theme(legend.title = element_text(size = 6), legend.text = element_text(size = 6)) + 
  theme(title = element_text(size = 10))+ theme(axis.text = element_text(size = 8)) +
  ggtitle("Marker genes for broad annotation") + theme(axis.text.x = element_text(angle = 90)) 

#Add idents of broad annotation 
Idents(merged_cl) <- "seurat_clusters"
current.cluster.ids <- c(0:26)
new.cluster.ids <- c("Kupffer","Kupffer","Endothelial","Endothelial","Endothelial","Endothelial","B","Kupffer","T_NK","B","Kupffer",
                     "Monocytes_DC","T_NK","Endothelial","Monocytes_DC","Monocytes_DC","Chol_Mets_Baso","Kupffer","Endothelial","Neutrophils",
                     "Hepatocytes","pDCs","Stellate_Stromal_Fib","Kupffer","Endothelial","RBC","B")
merged_cl$annotation <- plyr::mapvalues(x = merged_cl$seurat_clusters, from = current.cluster.ids, to = new.cluster.ids)

#####Fine grade annotation 
###Monocytes/DCs 
Idents(merged_cl) <- "annotation"
merged_cl_subCl <- FindSubCluster(merged_cl, "Monocytes_DC",graph.name = "RNA_snn", subcluster.name = "sub.cluster",resolution = 0.3, algorithm = 4)
DimPlot(merged_cl_subCl, reduction = "umap", label = TRUE, group.by = "sub.cluster")

Idents(merged_cl_subCl) <- "sub.cluster"
sub_celltype <- subset(merged_cl_subCl,idents = c( "Monocytes_DC_1","Monocytes_DC_2","Monocytes_DC_3",
                                                   "Monocytes_DC_4","Monocytes_DC_5",
                                                   "Monocytes_DC_6"))
DimPlot(sub_celltype, reduction = "umap", label = TRUE, group.by = "sub.cluster")

p <- DotPlot(sub_celltype, features = c("S100a4", "Itgax", "Crip1","Naaa","Adgre1","Axl","Mafb","Cx3cr1","C5ar1","Chil3","Sell","Gm9733", # general monocyte-derived cells
                                        "Pglyrp1","Spn","Trem3", #Patrolling monocytes: cluster  3
                                        "Clec9a","Cd24a", #cDC1: cluster 5
                                        "Cd209a","Mgl2","Clec10a","Cd7","Tnfsf9" , #cDC2: cluster 4
                                        "Siglech","Cox6a2",#pDCs: cluster 
                                        "Clec4f","Vsig4", #Kupffer cells: cluster
                                        "Ly6c2","Lyz2", #Ly6c+ macrophages: cluster 1,6
                                        "C1qc","C1qb","C1qa" #C1q+ macrophages: cluster 2
)) 

p + theme(legend.title = element_text(size = 6), legend.text = element_text(size = 6)) + 
  theme(title = element_text(size = 10))+ theme(axis.text = element_text(size = 6)) +
  ggtitle("Marker genes Monocytes/DC") + theme(axis.text.x = element_text(angle = 90)) 

#rename clusters that could be annotated 
Idents(merged_cl_subCl) <- "sub.cluster"
current.cluster.ids <- c("B","Chol_Mets_Baso","Endothelial","Hepatocytes","Kupffer","Monocytes_DC_1",
                         "Monocytes_DC_2","Monocytes_DC_3","Monocytes_DC_4",
                         "Monocytes_DC_5","Monocytes_DC_6",
                         "Neutrophils","pDCs", "RBC","Stellate_Stromal_Fib","T_NK")
new.cluster.ids <- c("B","Chol_Mets_Baso","Endothelial","Hepatocytes","Kupffer","Mac_Ly6c",
                     "Mac_C1q","Mono_patrolling","cDC2",
                     "cDC1","Mac_Ly6c",
                     "Neutrophils","pDC", "RBC","Stellate_Stromal_Fib","T_NK")
merged_cl_subCl$annotation <- plyr::mapvalues(x = merged_cl_subCl$sub.cluster, from = current.cluster.ids, to = new.cluster.ids)
DimPlot(merged_cl_subCl, reduction = "umap", label = TRUE, group.by = "annotation", label.size = 3)

###Metastatic cells/Cholangiocytes/Basophils
merged_cl <- merged_cl_subCl
Idents(merged_cl) <- "annotation"
merged_cl_subCl <- FindSubCluster(merged_cl, "Chol_Mets_Baso",graph.name = "RNA_snn", subcluster.name = "sub.cluster",resolution = 0.05, algorithm = 4)
DimPlot(merged_cl_subCl, reduction = "umap", label = TRUE, group.by = "sub.cluster")

Idents(merged_cl_subCl) <- "sub.cluster"
sub_celltype <- subset(merged_cl_subCl,idents = c( "Chol_Mets_Baso_1","Chol_Mets_Baso_2",
                                                   "Chol_Mets_Baso_3","Chol_Mets_Baso_4"))
DimPlot(sub_celltype, reduction = "umap", label = TRUE, group.by = "sub.cluster") 

p <- DotPlot(sub_celltype, features = c(  "Krt19", "Epcam","Sox9","Spp1", #Cholangiocytes: cluster  3
                                          "Sprr2a2","Pglyrp1","Gpx2", #Metastatic cells: clusters 1
                                          "Fcer1a","Cpa3"#Basophils: cluster 2
)) #4 mixed 
p + theme(legend.title = element_text(size = 6), legend.text = element_text(size = 6)) + 
  theme(title = element_text(size = 10))+ theme(axis.text = element_text(size = 8)) + theme(axis.text.x = element_text(angle = 90)) 

#rename clusters that could be annotated 
Idents(merged_cl) <- "annotation"
current.cluster.ids <- c("B","cDC1","cDC2", "Chol_Mets_Baso_1","Chol_Mets_Baso_2",
                         "Chol_Mets_Baso_3","Chol_Mets_Baso_4","Endothelial","Hepatocytes","Kupffer","Mac_C1q","Mac_Ly6c",
                         "Mono_patrolling",
                         "Neutrophils","pDC", "RBC","Stellate_Stromal_Fib", "T_NK")
new.cluster.ids <- c("B","cDC1","cDC2", "Metastasis","Basophils",
                     "Cholangiocytes","mixed","Endothelial","Hepatocytes","Kupffer","Mac_C1q","Mac_Ly6c",
                     "Mono_patrolling",
                     "Neutrophils","pDC", "RBC","Stellate_Stromal_Fib", "T_NK")
merged_cl_subCl$annotation <- plyr::mapvalues(x = merged_cl_subCl$sub.cluster, from = current.cluster.ids, to = new.cluster.ids)
DimPlot(merged_cl_subCl, reduction = "umap", label = TRUE, group.by = "annotation", label.size = 3) 

###Endothelial cells 
merged_cl <- merged_cl_subCl
Idents(merged_cl) <- "annotation"
merged_cl_subCl <- FindSubCluster(merged_cl, "Endothelial",graph.name = "RNA_snn", subcluster.name = "sub.cluster",resolution = 0.3, algorithm = 4)
DimPlot(merged_cl_subCl, reduction = "umap", label = TRUE, group.by = "sub.cluster")

Idents(merged_cl_subCl) <- "sub.cluster"
sub_celltype <- subset(merged_cl_subCl, idents = c( "Endothelial_1","Endothelial_2","Endothelial_3","Endothelial_4",
                                                    "Endothelial_5","Endothelial_6",
                                                    "Endothelial_7"))
DimPlot(sub_celltype, reduction = "umap", label = TRUE, group.by = "sub.cluster", label.size = 3) 

p <- DotPlot(sub_celltype, features = c("Lyve1","Flt4","Stab2", #LSECs: cluster 1,2,3,4,6
                                        "Il33","Pdgfb","Cd9","Timp3", #LVECs: cluster 5,7
                                        "Rspo3","Wnt9b", #Central vein ECs
                                        "Gja5", #Protain vein ECs
                                        "Ccl21a","Mmrn1","Thy1", #LECs (lymphatic ECs)
                                        "Clec4f","Vsig4" #Kupffer cells
)) 

p + theme(legend.title = element_text(size = 6), legend.text = element_text(size = 6)) + 
  theme(title = element_text(size = 10))+ theme(axis.text = element_text(size = 8)) + theme(axis.text.x = element_text(angle = 90)) 

#rename clusters that could be annotated 
Idents(merged_cl_subCl) <- "sub.cluster"
current.cluster.ids <- c("B","Basophils","cDC1","cDC2", "Cholangiocytes","Endothelial_1","Endothelial_2","Endothelial_3",
                         "Endothelial_4","Endothelial_5","Endothelial_6","Endothelial_7", "Hepatocytes","Kupffer","Mac_C1q", "Mac_Ly6c",
                         "Metastasis","mixed",
                         "Mono_patrolling",
                         "Neutrophils","RBC","Stellate_Stromal_Fib", "T_NK")
new.cluster.ids <- c("B","Basophils","cDC1","cDC2", "Cholangiocytes","LSECs","LSECs","LSECs",
                     "LSECs","LVECs","LSECs","LVECs", "Hepatocytes","Kupffer","Mac_C1q", "Mac_Ly6c",
                     "Metastasis","mixed",
                     "Mono_patrolling",
                     "Neutrophils","RBC","Stellate_Stromal_Fib", "T_NK")
merged_cl_subCl$annotation <- plyr::mapvalues(x = merged_cl_subCl$sub.cluster, from = current.cluster.ids, to = new.cluster.ids)
DimPlot(merged_cl_subCl, reduction = "umap", label = TRUE, group.by = "annotation", label.size = 3)

###T and NK cells
merged_cl <- merged_cl_subCl
Idents(merged_cl) <- "annotation"
merged_cl_subCl <- FindSubCluster(merged_cl, "T_NK",graph.name = "RNA_snn", subcluster.name = "sub.cluster",resolution = 0.4, algorithm = 4)
DimPlot(merged_cl_subCl, reduction = "umap", label = TRUE, group.by = "sub.cluster")

Idents(merged_cl_subCl) <- "sub.cluster"
sub_celltype <- subset(merged_cl_subCl,idents = c("T_NK_1","T_NK_2","T_NK_3","T_NK_4","T_NK_5","T_NK_6"))
DimPlot(sub_celltype, reduction = "umap", label = TRUE, group.by = "sub.cluster", label.size = 3) 

#Cd59b = Itga2, Cd11b = Itgam 
#Cd59b = Itga2, Cd11b = Itgam, Cd11c = Itgax, Cd11a = Itgal, Cd278 = Icos, Cd73 = Nt5e, Cd62l = Sell , Cd43 = Spn
p <- DotPlot(sub_celltype, features = c( "Trac", "Trbc1", "Trbc2","Cd3d","Themis" ,#T cells general 
                                           "Klra8","Cma1","Ncr1","Itga2","Itgam",#NK cells
                                           "Itgax", #NKT cells: cluster 4
                                           "Gzmc","Klrb1b", #ILC1s 
                                           "Cd4", #CD4+ T cells: cluster 3,1
                                           "Itgal", #Th1s
                                           "Icos","Nt5e","Itgb8", #Th17
                                           "Cd8a","Sell", #CD+8 T cells: cluster 2,5
                                           "Ly6c1","Spn" #TEMs
  ))
#6 mixed 
p + theme(legend.title = element_text(size = 6), legend.text = element_text(size = 6)) + 
  theme(title = element_text(size = 10))+ theme(axis.text = element_text(size = 8)) + theme(axis.text.x = element_text(angle = 90)) 

#rename clusters that could be annotated 
Idents(merged_cl_subCl) <- "sub.cluster"
current.cluster.ids <- c("B","Basophils","cDC1","cDC2", "Cholangiocytes","Hepatocytes","Kupffer", "LSECs","LVECs",
                         "Mac_C1q", "Mac_Ly6c","Metastasis","mixed",
                         "Mono_patrolling",
                         "Neutrophils","pDC", "RBC","Stellate_Stromal_Fib", "T_NK_1", "T_NK_2", "T_NK_3", "T_NK_4", "T_NK_5", 
                         "T_NK_6")
new.cluster.ids <- c("B","Basophils","cDC1","cDC2", "Cholangiocytes","Hepatocytes","Kupffer", "LSECs","LVECs",
                     "Mac_C1q", "Mac_Ly6c","Metastasis","mixed",
                     "Mono_patrolling",
                     "Neutrophils","pDC", "RBC","Stellate_Stromal_Fib", "T_CD4", "T_CD8", "T_CD4", "NKT", "T_CD8", 
                     "mixed")
merged_cl_subCl$annotation <- plyr::mapvalues(x = merged_cl_subCl$sub.cluster, from = current.cluster.ids, to = new.cluster.ids)
DimPlot(merged_cl_subCl, reduction = "umap", label = TRUE, group.by = "annotation", label.size = 3)

###B cells
merged_cl <- merged_cl_subCl
Idents(merged_cl) <- "annotation"
merged_cl_subCl <- FindSubCluster(merged_cl, "B",graph.name = "RNA_snn", subcluster.name = "sub.cluster",resolution = 0.1, algorithm = 4)
DimPlot(merged_cl_subCl, reduction = "umap", label = TRUE, group.by = "sub.cluster")

Idents(merged_cl_subCl) <- "sub.cluster"
sub_celltype <- subset(merged_cl_subCl,idents = c("B_1","B_2","B_3"))
DimPlot(sub_celltype, reduction = "umap", label = TRUE, group.by = "sub.cluster", label.size = 3) 

#https://www.abcam.com/primary-antibodies/b-cells-basic-immunophenotyping
p <- DotPlot(sub_celltype, features = c("Igkc", "Cd79b","Cd79a","Pax5",#B cells general
                                        "Cd19","Cd34","Cd38",#naiv B, (express no IgM)
                                        "Ighm","Ms4a1","Cd40",#memory B cells: clusters 1,2
                                        "Igha" ,"Cd27", "Sdc1"#Plasma cells: cluster 3
)) 
p + theme(legend.title = element_text(size = 6), legend.text = element_text(size = 6)) + 
  theme(title = element_text(size = 10))+ theme(axis.text = element_text(size = 8)) + theme(axis.text.x = element_text(angle = 90)) 

#rename clusters that could be annotated 
Idents(merged_cl_subCl) <- "sub.cluster"
current.cluster.ids <- c("B_1","B_2","B_3","Basophils","cDC1","cDC2", "Cholangiocytes","Hepatocytes","Kupffer", "LSECs","LVECs",
                         "Mac_C1q", "Mac_Ly6c","Metastasis","mixed",
                         "Mono_patrolling",
                         "Neutrophils","NKT","pDC", "RBC","Stellate_Stromal_Fib", "T_CD4","T_CD8")
new.cluster.ids <- c("B_mem","B_mem","B_plasma","Basophils","cDC1","cDC2", "Cholangiocytes","Hepatocytes","Kupffer", "LSECs","LVECs",
                     "Mac_C1q", "Mac_Ly6c","Metastasis","mixed",
                     "Mono_patrolling",
                     "Neutrophils","NKT","pDC", "RBC","Stellate_Stromal_Fib", "T_CD4","T_CD8")
merged_cl_subCl$annotation <- plyr::mapvalues(x = merged_cl_subCl$sub.cluster, from = current.cluster.ids, to = new.cluster.ids)
DimPlot(merged_cl_subCl, reduction = "umap", label = TRUE, group.by = "annotation", label.size = 3)

###Fibroblasts and Stellate cells
merged_cl <- merged_cl_subCl
Idents(merged_cl) <- "annotation"
merged_cl_subCl <- FindSubCluster(merged_cl, c("Stellate_Stromal_Fib"),graph.name = "RNA_snn", subcluster.name = "sub.cluster",resolution = 0.4, algorithm = 4)
DimPlot(merged_cl_subCl, reduction = "umap", label = TRUE, group.by = "sub.cluster")

Idents(merged_cl_subCl) <- "sub.cluster"
sub_celltype <- subset(merged_cl_subCl,idents = c("Stellate_Stromal_Fib_1","Stellate_Stromal_Fib_2","Stellate_Stromal_Fib_3","Stellate_Stromal_Fib_4"))
DimPlot(sub_celltype, reduction = "umap", label = TRUE, group.by = "sub.cluster", label.size = 3) 

p <- DotPlot(sub_celltype, features = c("Lrat", "Rbp1","Nt5e","Itgb1","Dach1","Nrxn1","Dcdc2a","Reln", #Stellate cells: cluster 4
                                        "Svep1","Ncam1" #Fibroblasts: clusters 2,3
                                        
))
#remove 1,2
p + theme(legend.title = element_text(size = 6), legend.text = element_text(size = 6)) + 
  theme(title = element_text(size = 10))+ theme(axis.text = element_text(size = 6)) + theme(axis.text.x = element_text(angle = 90)) 

#rename clusters that could be annotated 
Idents(merged_cl_subCl) <- "sub.cluster"
current.cluster.ids <- c("B_mem","B_plasma","Basophils","cDC1","cDC2", "Cholangiocytes","Hepatocytes","Kupffer", "LSECs","LVECs",
                         "Mac_C1q", "Mac_Ly6c","Metastasis","mixed",
                         "Mono_patrolling",
                         "Neutrophils","NKT","pDC", "RBC","Stellate_Stromal_Fib_1","Stellate_Stromal_Fib_2", "Stellate_Stromal_Fib_3",
                         "Stellate_Stromal_Fib_4","T_CD4","T_CD8")
new.cluster.ids <- c("B_mem","B_plasma","Basophils","cDC1","cDC2", "Cholangiocytes","Hepatocytes","Kupffer", "LSECs","LVECs",
                     "Mac_C1q", "Mac_Ly6c","Metastasis","mixed",
                     "Mono_patrolling",
                     "Neutrophils","NKT","pDC", "RBC","mixed2","Fibroblasts", "Fibroblasts",
                     "Stellate","T_CD4","T_CD8")
merged_cl_subCl$annotation <- plyr::mapvalues(x = merged_cl_subCl$sub.cluster, from = current.cluster.ids, to = new.cluster.ids)
DimPlot(merged_cl_subCl, reduction = "umap", label = TRUE, group.by = "annotation", label.size = 3)

###Kupffer cells
merged_cl <- merged_cl_subCl
Idents(merged_cl) <- "annotation"
merged_cl_subCl <- FindSubCluster(merged_cl, c("Kupffer"),graph.name = "RNA_snn", subcluster.name = "sub.cluster",resolution = 0.2, algorithm = 4)
DimPlot(merged_cl_subCl, reduction = "umap", label = TRUE, group.by = "sub.cluster")

Idents(merged_cl_subCl) <- "sub.cluster"
sub_celltype <- subset(merged_cl_subCl,idents = c("Kupffer_1","Kupffer_2","Kupffer_3","Kupffer_4","Kupffer_5"))
DimPlot(sub_celltype, reduction = "umap", label = TRUE, group.by = "sub.cluster", label.size = 3) 

#KC1 (ESAM- CD206 low) and KC2 (ESAM+ CD206 high) identity; CD206 = Mrc1, 
#https://www.sciencedirect.com/science/article/pii/S1074761321002119?via%3Dihub
#https://www.sciencedirect.com/science/article/pii/S1074761321003368?via%3Dihub 
#(KC1 = CD206(low), ESAM(negative) and KC2 = CD206(high), ESAM(positive) 
p <- DotPlot(sub_celltype, features = c("Clec4f","Vsig4","Timd4","Itgal","Cd5l", "Slc16a9","Slc40a1", #Kupffer cells (KC1): cluster 1,3,5
                                        "Esam","Mrc1", "Cd36","Cd63", "Cd81", "Lamp1", #Kupffer cells 2 (KC2) 
                                        "Ldb2","Rasip1","Ddx19b","Champ1","Tmem88","Cxcr4","Igfbp7","Clec4g" #Endothelial
)) 
#2 Kupffer_Endo, 4 mixed 
p + theme(legend.title = element_text(size = 6), legend.text = element_text(size = 6)) + 
  theme(title = element_text(size = 10))+ theme(axis.text = element_text(size = 6)) + theme(axis.text.x = element_text(angle = 90)) 

#rename clusters that could be annotated 
Idents(merged_cl_subCl) <- "sub.cluster"
current.cluster.ids <- c("B_mem","B_plasma","Basophils","cDC1","cDC2", "Cholangiocytes","Hepatocytes","Fibroblasts", "Kupffer_1",
                         "Kupffer_2","Kupffer_3","Kupffer_4","Kupffer_5", "LSECs","LVECs",
                         "Mac_C1q", "Mac_Ly6c","Metastasis","mixed","mixed2",
                         "Mono_patrolling",
                         "Neutrophils","NKT","pDC", "RBC", "Stellate", "T_CD4","T_CD8")
new.cluster.ids <- c("B_mem","B_plasma","Basophils","cDC1","cDC2", "Cholangiocytes","Hepatocytes","Fibroblasts", "Kupffer",
                     "Kupffer_Endo","Kupffer","mixed3","Kupffer", "LSECs","LVECs",
                     "Mac_C1q", "Mac_Ly6c","Metastasis","mixed","mixed2",
                     "Mono_patrolling",
                     "Neutrophils","NKT","pDC", "RBC", "Stellate", "T_CD4","T_CD8")
merged_cl_subCl$annotation <- plyr::mapvalues(x = merged_cl_subCl$sub.cluster, from = current.cluster.ids, to = new.cluster.ids)
DimPlot(merged_cl_subCl, reduction = "umap", label = TRUE, group.by = "annotation", label.size = 3)

#####assess mixed clusters for biases 
Idents(merged_cl_subCl) <- "annotation"
sub_celltype <- subset(merged_cl_subCl,idents = c("mixed","mixed2","mixed3"))
DimPlot(sub_celltype, reduction = "umap", label = TRUE, group.by = "sub.cluster", label.size = 3) 
DimPlot(sub_celltype, group.by = "seurat_clusters",split.by = "protocol")
#they are present in both conditions! no Bias 

Idents(sub_celltype) <- "seurat_clusters"
p <- DotPlot(sub_celltype, features = c(
  "Ighm", "Igkc", #B cells
  "Trac","Cd3d" ,#T cells
  "Klra8","Cma1" ,#NK cells
  "Clec4f","Vsig4", #Kupffer cells
  "Bcl11a", "Ccr9", #DCs
  "S100a4", "Itgax", #general Monocytes
  "Hba-a2", #Red blood cells
  "Csf3r", "S100a8",#Neutrophils
  "Siglech","Cox6a2", #pDCs
  "Fcer1a","Cpa3",#Basophils
  "Ttr", "Alb", "Cyp2e1", #Hepatocytes
  "Pecam1","Dll4", #LSECs
  "Lrat","Reln", #Stellate cells
  "Carmn","Nr1h5",  #Stromal cells
  "Svep1","Ncam1", #Fibroblasts
  "Krt19", "Epcam","Sox9", #Cholangiocytes
  "Sprr2a2","Pglyrp1","Gpx2" #Metastatic cells
)) 
p + theme(legend.title = element_text(size = 6), legend.text = element_text(size = 6)) + 
  theme(title = element_text(size = 10))+ theme(axis.text = element_text(size = 8)) +
  ggtitle("Marker genes for broad annotation") + theme(axis.text.x = element_text(angle = 90)) 
#they are MIXED 

#####Remove low quality cells (annotated as mixed) 
merged_cl <- merged_cl_subCl
Idents(merged_cl) <- "annotation"
merged_cl_sub <- subset(merged_cl, idents = c("B_mem","B_plasma","Basophils","cDC1","cDC2", "Cholangiocytes","Fibroblasts", "Hepatocytes",
                                              "Kupffer",
                                              "Kupffer_Endo", "LSECs","LVECs",
                                              "Mac_C1q", "Mac_Ly6c","Metastasis",
                                              "Mono_patrolling",
                                              "Neutrophils","NKT","pDC", "Stellate", "T_CD4","T_CD8"))
DimPlot(merged_cl_sub, reduction = "umap", label = TRUE, group.by = "annotation", label.size = 3)

#####Plot Annotated umap 
###colors 
#B cells (orange): B_mem: #EF975B; B_plasma: #F4741E    
#Granulocytes (yellow): Basophils:#56595B; Neutrophils: #AEAEAF  
#Endothelial (grey): LVECs: #F4CB1C; LSECs:#F4E740
#Stromal cells (brown): Fibroblasts: #8E5229; Stellate cells: #C1A08A
#Hepatocytes (redish): #F46042
#DCs (turquoise): cDC1: #91C6C4; cDC2: #315B5A; pDC: #7AEDD9  
#KC (purple): Kupffer: #E20FE8; Kupffer_Endo: #C491C6  
#Monocytes (green): M_C1q: #19E80F; M_patrolling: #B5EDB2; Mac_Ly6c: #566B44  
#T (blue): CD4: #2323E5 ; CD8: #9A9AE5; NKT: #2C91E5; NK:#95A0C6   
#Metastatic cells (red): #F90606
#Cholangiocytes (pink): #EDB2D4

p <- DimPlot(merged_cl_sub, label = TRUE,label.size = 5, group.by = "annotation", pt.size = 0.5, split.by = "protocol",
             cols = c("#EF975B","#F4741E","#56595B","#91C6C4",
                      "#315B5A","#EDB2D4",
                      "#8E5229","#F46042","#E20FE8",
                      "#C491C6","#F4E740", "#F4CB1C", "#19E80F",
                      "#566B44","#F90606","#B5EDB2","#AEAEAF",
                      "#2C91E5","#7AEDD9", "#C1A08A", "#2323E5","#9A9AE5")) + 
  theme(legend.title = element_text(size =10), legend.text = element_text(size = 10)) + 
  theme(title = element_text(size = 10))+ theme(axis.text = element_text(size = 10)) +
  ggtitle("Annotation") 
p + ggsave("./figures/8/Annotated_umap_sc_comp_SpS.pdf", width = 15, height = 10)
p + ggsave("./figures/8/Annotated_umap_sc_comp_SpS.svg", width = 15, height = 10)

#####extract only Sc and plot marker genes in heatmap 
Idents(merged_cl_sub) <- "protocol"
sc <- subset(merged_cl_sub, idents = "sc")

marker_genes <- c("Ttr", "Alb", "Cyp2e1", #Hepatocytes 
                  "Sprr2a2","Pglyrp1","Gpx2", #Metastatic cells 
                  "Krt19", "Epcam","Sox9", #Cholangiocytes  
                  
                  #Monocytes
                  "S100a4", "Itgax","Naaa","Adgre1", # #general Monocyte derived cells 
                  "Pglyrp1","Spn","Trem3", #Patrolling monocytes  
                  "Ly6c2","Lyz2", #Ly6c+ macrophages 
                  "C1qc","C1qb","C1qa", #C1q+ macrophages 
                  
                  #DCs
                  "Clec9a","Cd24a", #cDC1 
                  "Cd209a","Mgl2","Clec10a","Cd7","Tnfsf9" , #cDC2 
                  "Siglech","Cox6a2",#pDCs 
                  
                  #Endothelial
                  "Lyve1","Flt4","Stab2", #LSECs 
                  "Il33","Pdgfb","Cd9","Timp3", #LVECs
                  
                  #T and NK cells 
                  "Trac", "Trbc1", "Trbc2","Cd3d","Themis" ,#T cells
                  "Klra8","Cma1","Ncr1","Itga2","Itgam",#NK cells
                  "Itgax", #NKT cells
                  "Cd4", #CD4+ T cells 
                  "Cd8a","Sell", #Cd8 T cells 
                  "Cma1","Ncr1","Itga2","Itgam", #NK cells 
                  
                  #B cells
                  "Igkc", "Cd79b","Cd79a", #B cells general  
                  "Ighm","Ms4a1","Cd40","Igha" ,#memory B cells 
                  "Igha" ,"Cd27", "Sdc1" ,#Plasma cells
                  
                  #Granulocytes
                  "Csf3r", "S100a8", #Neutrophils 
                  "Fcer1a","Cpa3","Il6", #Basophils 
                  
                  #Stromal cells
                  "Lrat", "Rbp1","Nt5e","Itgb1","Dach1","Nrxn1","Dcdc2a","Reln", #Stellate cells
                  "Svep1","Ncam1", #Fibroblasts 
                  
                  #Kupffer cells
                  "Clec4f","Vsig4","Timd4","Itgal","Cd5l", "Slc16a9","Slc40a1"
)

# Re-level object@ident
sc@active.ident <- factor(x = sc@active.ident, 
                                     levels = c("Hepatocytes","Metastasis","Cholangiocytes","Mac_C1q" ,"Mono_patrolling","Mac_Ly6c",
                                                "cDC1","cDC2","pDC",
                                                "LSECs", "LVECs",
                                                "T_CD4", "T_CD8","NKT",
                                                "B_mem","B_plasma","Neutrophils",
                                                "Basophils","Fibroblasts","Stellate", "Kupffer","Kupffer_Endo"))
Idents(sc) <- "annotation"
p <- DoHeatmap(subset(sc, downsample = 100), features = marker_genes, assay = "RNA",slot = "data") +
  theme(legend.title = element_text(size = 15), legend.text = element_text(size = 15)) + theme(title = element_text(size = 25))+ 
  theme(axis.text = element_text(size = 7)) +scale_fill_gradientn(colours = colorRampPalette(rev(brewer.pal(n = 15, name ="RdBu")))(256)) +
  ggtitle("Heatmap marker genes used for annotation")
p + ggsave("./figures/8/Heatmap_anno.png", width = 15, height = 10)
p + ggsave("./figures/8/Heatmap_anno.svg", width = 15, height = 10)

########## Comparison of cellular quality sc vs. SpS ##########
#downsampled to 50K reads/cell, cells that have more are downsampled, cells that do not reach this threshold are removed 
#####annotate the zUMI output 
dge <- readRDS("/home/khandler/NAS/Kristina/Mets_sc1.dgecounts.rds")
Msc1 <- Annotation_mouse(dge)
Msc1_S <- CreateSeuratObject(Msc1, "sc1",min.cells = 3, min.features = 200)

dge <- readRDS("/home/khandler/NAS/Kristina/Mets_sc2.dgecounts.rds")
Msc2 <- Annotation_mouse(dge)
Msc2_S <- CreateSeuratObject(Msc2, "sc2",min.cells = 3, min.features = 200)

########## Investigation of quality measuremnts after downasmpling to 30K reads/cell ###################
#downsampled to 30K reads/cell, cells that have more reads are downsampled, cells that do not reach this threshold are removed 
#####annotate the zUMI output 
dge <- readRDS("/home/khandler/NAS/Kristina/Mets_sc1.dgecounts.rds")
Msc1 <- Annotation_mouse_downsampled(dge)
Msc1_S <- CreateSeuratObject(Msc1, "sc1",min.cells = 3, min.features = 200)

dge <- readRDS("/home/khandler/NAS/Kristina/Mets_sc2.dgecounts.rds")
Msc2 <- Annotation_mouse_downsampled(dge)
Msc2_S <- CreateSeuratObject(Msc2, "sc2",min.cells = 3, min.features = 200)

dge <- readRDS("/home/khandler/NAS/Kristina/SpSM1.dgecounts.rds")
S1 <- Annotation_mouse_downsampled(dge)
S1 <- CreateSeuratObject(S1, "M1",min.cells = 3, min.features = 200)

dge <- readRDS("/home/khandler/NAS/Kristina/SpSM3.dgecounts.rds")
S2 <- Annotation_mouse_downsampled(dge)
S2 <- CreateSeuratObject(S2, "M3",min.cells = 3, min.features = 200)

dge <- readRDS("/home/khandler/NAS/Kristina/SpSM4.dgecounts.rds")
S3 <- Annotation_mouse_downsampled(dge)
S3 <- CreateSeuratObject(S3, "M4",min.cells = 3, min.features = 200)

dge <- readRDS("/home/khandler/NAS/Kristina/SpSM5.dgecounts.rds")
S4 <- Annotation_mouse_downsampled(dge)
S4 <- CreateSeuratObject(S4, "M5",min.cells = 3, min.features = 200)

dge <- readRDS("/home/khandler/NAS/Kristina/SpS4M1.dgecounts.rds")
S5 <- Annotation_mouse_downsampled(dge)
S5 <- CreateSeuratObject(S5, "4M1",min.cells = 3, min.features = 200)

dge <- readRDS("/home/khandler/NAS/Kristina/SpS6M1.dgecounts.rds")
S7 <- Annotation_mouse_downsampled(dge)
S7 <- CreateSeuratObject(S7, "6M1",min.cells = 3, min.features = 200)

dge <- readRDS("/home/khandler/NAS/Kristina/SpS6M2.dgecounts.rds")
S8 <- Annotation_mouse_downsampled(dge)
S8 <- CreateSeuratObject(S8, "6M2",min.cells = 3, min.features = 200)

dge <- readRDS("/home/khandler/NAS/Kristina/SpS6M3.dgecounts.rds")
S9 <- Annotation_mouse_downsampled(dge)
S9 <- CreateSeuratObject(S9, "6M3",min.cells = 3, min.features = 200)

dge <- readRDS("/home/khandler/NAS/Kristina/SpSWT.dgecounts.rds")
S10 <- Annotation_mouse_downsampled(dge)
S10 <- CreateSeuratObject(S10, "WT",min.cells = 3, min.features = 200)

#####Ratio mito/cytopl genes 
mito.features <- grep(pattern = "^mt-", x = rownames(x = Msc1_S), value = TRUE)
percent.mito <- Matrix::colSums(x = GetAssayData(object = Msc1_S, slot = 'counts')[mito.features, ]) / Matrix::colSums(x = GetAssayData(object = Msc1_S, slot = 'counts'))
Msc1_S <- AddMetaData(object = Msc1_S, metadata = percent.mito, col.name = "percent.mito")

mito.features <- grep(pattern = "^mt-", x = rownames(x = Msc2_S), value = TRUE)
percent.mito <- Matrix::colSums(x = GetAssayData(object = Msc2_S, slot = 'counts')[mito.features, ]) / Matrix::colSums(x = GetAssayData(object = Msc2_S, slot = 'counts'))
Msc2_S <- AddMetaData(object = Msc2_S, metadata = percent.mito, col.name = "percent.mito")

mito.features <- grep(pattern = "^mt-", x = rownames(x = S1), value = TRUE)
percent.mito <- Matrix::colSums(x = GetAssayData(object = S1, slot = 'counts')[mito.features, ]) / Matrix::colSums(x = GetAssayData(object = S1, slot = 'counts'))
S1 <- AddMetaData(object = S1, metadata = percent.mito, col.name = "percent.mito")

mito.features <- grep(pattern = "^mt-", x = rownames(x = S2), value = TRUE)
percent.mito <- Matrix::colSums(x = GetAssayData(object = S2, slot = 'counts')[mito.features, ]) / Matrix::colSums(x = GetAssayData(object = S2, slot = 'counts'))
S2 <- AddMetaData(object = S2, metadata = percent.mito, col.name = "percent.mito")

mito.features <- grep(pattern = "^mt-", x = rownames(x = S3), value = TRUE)
percent.mito <- Matrix::colSums(x = GetAssayData(object = S3, slot = 'counts')[mito.features, ]) / Matrix::colSums(x = GetAssayData(object = S3, slot = 'counts'))
S3 <- AddMetaData(object = S3, metadata = percent.mito, col.name = "percent.mito")

mito.features <- grep(pattern = "^mt-", x = rownames(x = S4), value = TRUE)
percent.mito <- Matrix::colSums(x = GetAssayData(object = S4, slot = 'counts')[mito.features, ]) / Matrix::colSums(x = GetAssayData(object = S4, slot = 'counts'))
S4 <- AddMetaData(object = S4, metadata = percent.mito, col.name = "percent.mito")

mito.features <- grep(pattern = "^mt-", x = rownames(x = S5), value = TRUE)
percent.mito <- Matrix::colSums(x = GetAssayData(object = S5, slot = 'counts')[mito.features, ]) / Matrix::colSums(x = GetAssayData(object = S5, slot = 'counts'))
S5 <- AddMetaData(object = S5, metadata = percent.mito, col.name = "percent.mito")

mito.features <- grep(pattern = "^mt-", x = rownames(x = S7), value = TRUE)
percent.mito <- Matrix::colSums(x = GetAssayData(object = S7, slot = 'counts')[mito.features, ]) / Matrix::colSums(x = GetAssayData(object = S7, slot = 'counts'))
S7 <- AddMetaData(object = S7, metadata = percent.mito, col.name = "percent.mito")

mito.features <- grep(pattern = "^mt-", x = rownames(x = S8), value = TRUE)
percent.mito <- Matrix::colSums(x = GetAssayData(object = S8, slot = 'counts')[mito.features, ]) / Matrix::colSums(x = GetAssayData(object = S8, slot = 'counts'))
S8 <- AddMetaData(object = S8, metadata = percent.mito, col.name = "percent.mito")

mito.features <- grep(pattern = "^mt-", x = rownames(x = S9), value = TRUE)
percent.mito <- Matrix::colSums(x = GetAssayData(object = S9, slot = 'counts')[mito.features, ]) / Matrix::colSums(x = GetAssayData(object = S9, slot = 'counts'))
S9 <- AddMetaData(object = S9, metadata = percent.mito, col.name = "percent.mito")

mito.features <- grep(pattern = "^mt-", x = rownames(x = S10), value = TRUE)
percent.mito <- Matrix::colSums(x = GetAssayData(object = S10, slot = 'counts')[mito.features, ]) / Matrix::colSums(x = GetAssayData(object = S10, slot = 'counts'))
S10 <- AddMetaData(object = S10, metadata = percent.mito, col.name = "percent.mito")

sc1_m <- median(Msc1_S@meta.data$percent.mito)
sc2_m <- median(Msc2_S@meta.data$percent.mito)
S1_m <- median(S1@meta.data$percent.mito)
S2_m <- median(S2@meta.data$percent.mito)
S3_m <- median(S3@meta.data$percent.mito)
S4_m <- median(S4@meta.data$percent.mito)
S5_m <- median(S5@meta.data$percent.mito)
S7_m <- median(S7@meta.data$percent.mito)
S8_m <- median(S8@meta.data$percent.mito)
S9_m <- median(S9@meta.data$percent.mito)
S10_m <- median(S10@meta.data$percent.mito)

a <- c(sc1_m,sc2_m,S1_m,S2_m,S3_m,S4_m,S5_m,S7_m,S8_m,S9_m,S10_m)
b <- c("sc1","sc2","SpS1","SpS2","SpS3","SpS4","Sp5", "SpS7","SpS8","SpS9","SpS10")
c <- c("sc","sc","SpS","SpS","SpS","SpS","SpS","SpS", "SpS","SpS","SpS")

df <- data.frame(a,b,c)
colnames(df) <- c("percent.mito","Sample","Type")

###plot Proportions to compare
p <- ggplot(df,aes(x = Type,y = percent.mito, fill = Type)) +theme_classic() +
  geom_bar(stat="summary",width=0.75, fun = "mean")  +
  stat_summary(fun.data = mean_se,  
               geom = "errorbar", width = 0.5) +
  geom_jitter(position = position_jitter(seed = 1, width =0.4),size = 3) + 
  theme(axis.text = element_text(size = 30))  + scale_y_continuous(limits=c(0, 0.1))+ 
  ggtitle("Percent mito") + xlab("Protocol") + 
  ylab("Percent mito") + theme(axis.title= element_text(size = 25)) + 
  theme(plot.title = element_text(size = 25, face = "bold"))  +
  theme(legend.title = element_text(size = 30), legend.text = element_text(size = 30)) +
  theme(legend.title = element_text(size = 30), legend.text = element_text(size = 30)) +
  ggsignif::geom_signif(comparisons = list(c("sc", "SpS")),
                        textsize=7,
                        test = "wilcox.test",map_signif_level = c("***"=0.001,"**"=0.01,"*"=0.05))
p + ggsave("./figures/8/Percent_mito_sc_SpS.pdf",width = 12, height = 10)  
p + ggsave("./figures/8/Percent_mito_sc_SpS.svg",width = 12, height = 10)  

#####nFeature 
sc1_m <- median(Msc1_S@meta.data$nFeature_RNA)
sc2_m <- median(Msc2_S@meta.data$nFeature_RNA)
S1_m <- median(S1@meta.data$nFeature_RNA)
S2_m <- median(S2@meta.data$nFeature_RNA)
S3_m <- median(S3@meta.data$nFeature_RNA)
S4_m <- median(S4@meta.data$nFeature_RNA)
S5_m <- median(S5@meta.data$nFeature_RNA)
S7_m <- median(S7@meta.data$nFeature_RNA)
S8_m <- median(S8@meta.data$nFeature_RNA)
S9_m <- median(S9@meta.data$nFeature_RNA)
S10_m <- median(S10@meta.data$nFeature_RNA)

a <- c(sc1_m,sc2_m,S1_m,S2_m,S3_m,S4_m,S5_m,S7_m,S8_m,S9_m,S10_m)
b <- c("sc1","sc2","SpS1","SpS2","SpS3","SpS4","SpS5", "SpS7","SpS8","SpS9","SpS10")
c <- c("sc","sc","SpS","SpS","SpS","SpS","SpS","SpS", "SpS","SpS","SpS")

df <- data.frame(a,b,c)
colnames(df) <- c("nFeature","Sample","Type")

###plot Proportions to compare
p <- ggplot(df,aes(x = Type,y = nFeature, fill = Type)) +theme_classic() +
 geom_bar(stat="summary",width=0.75, fun = "mean")  +
  stat_summary(fun.data = mean_se,  
               geom = "errorbar", width = 0.5) + 
  geom_jitter(position = position_jitter(seed = 1, width =0.4),size = 3) + 
  theme(axis.text = element_text(size = 30)) + scale_y_continuous(limits=c(0, 6000))  + 
  ggtitle("nFeature RNA") + xlab("Protocol") + 
  ylab("nFeature RNA") + theme(axis.title= element_text(size = 25)) + 
  theme(plot.title = element_text(size = 25, face = "bold")) +
  theme(legend.title = element_text(size = 30), legend.text = element_text(size = 30)) +
  theme(legend.title = element_text(size = 30), legend.text = element_text(size = 30)) + 
  ggsignif::geom_signif(comparisons = list(c("sc", "SpS")),
                        textsize=7,
                        test = "wilcox.test",map_signif_level = c("***"=0.001,"**"=0.01,"*"=0.05))
p + ggsave("./figures/8/nFeature_sc_SpS.pdf",width = 12, height = 10)  
p + ggsave("./figures/8/nFeature_sc_SpS.svg",width = 12, height = 10)  

#####UMI counts 
sc1_m <- median(Msc1_S@meta.data$nCount_RNA)
sc2_m <- median(Msc2_S@meta.data$nCount_RNA)
S1_m <- median(S1@meta.data$nCount_RNA)
S2_m <- median(S2@meta.data$nCount_RNA)
S3_m <- median(S3@meta.data$nCount_RNA)
S4_m <- median(S4@meta.data$nCount_RNA)
S5_m <- median(S5@meta.data$nCount_RNA)
S7_m <- median(S7@meta.data$nCount_RNA)
S8_m <- median(S8@meta.data$nCount_RNA)
S9_m <- median(S9@meta.data$nCount_RNA)
S10_m <- median(S10@meta.data$nCount_RNA)

a <- c(sc1_m,sc2_m,S1_m,S2_m,S3_m,S4_m,S5_m,S7_m,S8_m,S9_m,S10_m)
b <- c("sc1","sc2","SpS1","SpS2","SpS3","SpS4","SpS5", "SpS7","SpS8","SpS9","SpS10")
c <- c("sc","sc","SpS","SpS","SpS","SpS","SpS","SpS", "SpS","SpS","SpS")

df <- data.frame(a,b,c)
colnames(df) <- c("nCounts","Sample","Type")

###plot Proportions to compare
p <- ggplot(df,aes(x = Type,y = nCounts, fill = Type)) +theme_classic() +
  geom_bar(stat="summary",width=0.75, fun = "mean")  +
  geom_jitter(position = position_jitter(seed = 1, width =0.4),size = 3) + 
  theme(axis.text = element_text(size = 30)) +
  stat_summary(fun.data = mean_se,  
               geom = "errorbar", width = 0.5) +
  ggtitle("nCount RNA") + xlab("species well") + 
  ylab("nCountRNA") + theme(axis.title= element_text(size = 25)) + 
  theme(plot.title = element_text(size = 25, face = "bold"))  + scale_colour_discrete("Species") +
  theme(legend.title = element_text(size = 30), legend.text = element_text(size = 30)) +
  theme(legend.title = element_text(size = 30), legend.text = element_text(size = 30))  + ylim(0,20000) + 
  ggsignif::geom_signif(comparisons = list(c("sc", "SpS")),
                        textsize=7,
                        test = "wilcox.test",map_signif_level = c("***"=0.001,"**"=0.01,"*"=0.05))
p + ggsave("./figures/8/nCount_sc_SpS.pdf",width = 12, height = 10)  
p + ggsave("./figures/8/nCount_sc_SpS.svg",width = 12, height = 10)  





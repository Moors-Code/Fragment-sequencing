########## Part 5: Crohn preliminary analysis ##########
#This part analyses preliminary Crohn's disease sphere-seq data  

########## Prepare environment ##########
###Setting the working directory 
setwd("/mnt/khandler/R_projects/Sphere-sequencing/Sphere-seq_analysis/")

###Load packages and functions 
source("./functions_and_packages/1.Packages.R")
source("./functions_and_packages/4.Functions_cells_per_sphere_cutoff.R")

###Load data 
CrohnP393 <- readRDS("./data_files_generated/SpS_Crohn_P387.Rda")
CrohnP387 <- readRDS("./data_files_generated/SpS_Crohn_P393.Rda")

########## merge and cluster data ##########
###add prefix to sphere for sample ID
CrohnP393$sphere <- paste(CrohnP393$sphere, "_1", sep="")
CrohnP387$sphere <- paste(CrohnP387$sphere, "_2", sep="")

###merge
merged <- merge(CrohnP393, CrohnP387)

###add percentage of mitochondrial features
mito.features <- grep(pattern = "^MT-", x = rownames(x = merged), value = TRUE)
percent.mito <- Matrix::colSums(x = GetAssayData(object = merged, slot = 'counts')[mito.features, ]) / Matrix::colSums(x = GetAssayData(object = merged, slot = 'counts'))
merged <- AddMetaData(object = merged, metadata = percent.mito, col.name = "percent.mito")
Idents(merged) <- "orig.ident"
FeatureScatter(merged, feature1 = "nCount_RNA", feature2 = "percent.mito")
VlnPlot(merged, features = "percent.mito") 
FeatureScatter(merged, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

###clustering 
merged <- SCTransform(merged, vars.to.regress = c("nCount_RNA","percent.mito"), verbose = FALSE)
merged <- RunPCA(object = merged, features = VariableFeatures(object = merged), npcs = 30, verbose = FALSE)
#check for batch effect  - no batch effect 
DimPlot(merged, reduction = "pca", group.by = "orig.ident")

#check for significant PCAs (10)
ElbowPlot(merged)

#cluster
merged <- FindNeighbors(object = merged, dims = 1:10, reduction = "pca")
merged_cl <- FindClusters(merged, resolution = 0.5, random.seed = 5, algorithm = 1, graph.name = "SCT_snn")
merged_cl <- RunUMAP(merged_cl, dims = 1:10, seed.use = 5, reduction = "pca")
DimPlot(merged_cl, reduction = "umap", label = TRUE, pt.size = .5, group.by = "seurat_clusters") 

########## Annotation based on marker gene from Martin et al 2019 and Smillie et al 2019 ##########
#first broadly annotate cells and then subcluster 
markers <- FindAllMarkers(object = merged_cl, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
View(markers %>% group_by(cluster) %>% top_n(n =10, wt = avg_log2FC))

#remove clusters with only MT genes in the first 5 top DEGs 
#cluster 7

DotPlot(merged_cl, features = c(
  "CD3D", "CD2", "CD7", # T cells cluster 
  "CCL5","CD8A",  #CD8 T: cluster 0,1
  "CD40LG","IL2","FOXP3","ICOS", #CD4 T : cluster 5
  "KLRC1","NKG7", #NK: cluster 13
  "TNFRSF17","MZB1", #Plasma cells: cluster 8,11,15,16,21,14
  "BANK1","CD79B","CD22","MS4A1", #B cells: cluster 6
  "FCN1","AIF1","S100A8","S100A9","IL1B","LYZ", #pro-inflammatory monocytes: cluster 9
  "C1QC","C1QA","TREM2", #Macrophages: cluster 3
  "S100A4","CD28", #Neutrophils 
  "HLA-DRA","HLA-DPA1", "HLA-DRB1","HLA-DQA1",#DCs: cluster  3
  "GZMB","LILRA4", #pDCs
  "TPSAB1","CMA1","KIT", #Mast cells: cluster 17
  "LYVE1","COL3A1","COL1A1","GPM6B","S100B", #Stroma/Glia: cluster 20
  "PLVAP","CLDN5","SEMA3G", #Endothelial: cluster 18
  "TFF1","TFF3","MUC2", #Goblet: cluster 18,12
  "SPOCK1","COL4A1","LGALS1", #Fibroblasts: cluster 18
  "EPCAM", #Epithelial: cluster 2,4,10
  "RGS5","COX4I2","IGFBP7", #Pericytes 
  "LGALS4","PIGR", #Transit amplyfier TA: cluster 2,4,10
  "AQP8","PLAC8","GUCA2A","CEACAM5" #Enterocytes: cluster 24,10
) )+   theme(axis.text.x = element_text(angle = 90)) 

current.cluster.ids <- c(0:21)
new.cluster.ids <- c("T","T","Stromal","Myeloid","Epithelial","T","B","remove",
                     "PC","Myeloid","Epithelial","PC","Epithelial", "T","PC","PC","PC","Mast","Stromal","remove","Stromal","PC")
merged_cl$annotation <- plyr::mapvalues(x = merged_cl$seurat_clusters, from = current.cluster.ids, to = new.cluster.ids)

#####Recluster for fine grade annotation 
###Myeloid 
Idents(merged_cl) <- "annotation"
merged_cl_subCl <- FindSubCluster(merged_cl, "Myeloid",graph.name = "SCT_snn", subcluster.name = "sub.cluster",resolution = 0.1, algorithm = 4)
DimPlot(merged_cl_subCl, reduction = "umap", label = TRUE, group.by = "sub.cluster")

Idents(merged_cl_subCl) <- "sub.cluster"
sub_celltype <- subset(merged_cl_subCl,idents = c( "Myeloid_1","Myeloid_2","Myeloid_3"))
DimPlot(sub_celltype, reduction = "umap", label = TRUE, group.by = "sub.cluster") 

markers <- FindAllMarkers(object = sub_celltype, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
View(markers %>% group_by(cluster) %>% top_n(n =10, wt = avg_log2FC))

DotPlot(sub_celltype, features = c(
  "FCN1","AIF1","S100A8","S100A9","IL1B","LYZ", #pro-inflammatory monocytes 
  "C1QC","C1QA","TREM2", #Macrophages 
  "HLA-DRA","HLA-DPA1", "HLA-DRB1","HLA-DQA1",#DCs 
  "GZMB","LILRA4" #pDCs
) )+   theme(axis.text.x = element_text(angle = 90)) 

#3 anti-inflammatory macrophages C1QC,C1QA,LYZ
#2 pro-inflammatory macrophages IL1B,S100A9,S100A8
#1 unclear --> remove 

Idents(merged_cl_subCl) <- "sub.cluster"
current.cluster.ids <- c("B","Epithelial", "Mast","Myeloid_1","Myeloid_2","Myeloid_3","PC" ,"remove","Stromal","T")
new.cluster.ids <- c("B","Epithelial", "Mast","remove","Pro_infl_Macs","Anti_infl_Macs","PC","remove","Stromal","T")
merged_cl_subCl$annotation <- plyr::mapvalues(x = merged_cl_subCl$sub.cluster, from = current.cluster.ids, to = new.cluster.ids)
DimPlot(merged_cl_subCl, reduction = "umap", label = TRUE, group.by = "annotation", label.size = 3)

###PCs check mito genes and remove cluster if high 
markers <- FindAllMarkers(object = merged_cl_subCl, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
View(markers %>% group_by(cluster) %>% top_n(n =10, wt = avg_log2FC))
#all good 

###Stromal 
Idents(merged_cl_subCl) <- "annotation"
merged_cl_subCl <- FindSubCluster(merged_cl_subCl, "Stromal",graph.name = "SCT_snn", subcluster.name = "sub.cluster",resolution = 0.1, algorithm = 4)
DimPlot(merged_cl_subCl, reduction = "umap", label = TRUE, group.by = "sub.cluster")

Idents(merged_cl_subCl) <- "sub.cluster"
sub_celltype <- subset(merged_cl_subCl,idents = c( "Stromal_1","Stromal_2","Stromal_3","Stromal_4","Stromal_5"))
DimPlot(sub_celltype, reduction = "umap", label = TRUE, group.by = "sub.cluster") 

markers <- FindAllMarkers(object = sub_celltype, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
View(markers %>% group_by(cluster) %>% top_n(n =10, wt = avg_log2FC))
#remove 2 --> high mito 

DotPlot(sub_celltype, features = c(
  "VWF","RAMP3","NPDC1","JAM2","PLVAP","CD36","ACKR1", #Endothelial cells 
  "LYVE1","TFF3","CCL21", #Lymphatics
  "RGS5","NDUFA4L2","ACTG2","MYH11", #Pericytes
  "CXCL14","LUM","ADAMDEC1","ABCA8", #SM
  "CTSK","MMP2","PTGDS", #Fibrobalsts
  "ALDH1A1","GPM6B","CRYAB","CLU","SPP1","MPZ", #Activated fibroblasts 
  "AQP8","PLAC8","GUCA2A","CEACAM5" ,"EPCAM",#Enterocytes 
  "S100A4" ,#Glial cells 
  "PLA2G2A","DEFA6","DEFA5","REG1A" #Paneth cells 2
) )+   theme(axis.text.x = element_text(angle = 90)) 

#remove 5 --> mixture of many cell types 
#stricture fibroblasts, express S100A4 --> cluster 4 
#https://journals.physiology.org/doi/full/10.1152/ajpgi.00351.2009?rfr_dat=cr_pub++0pubmed&url_ver=Z39.88-2003&rfr_id=ori%3Arid%3Acrossref.org
#3 Endothelial 
#1 Enterocytes

Idents(merged_cl_subCl) <- "sub.cluster"
current.cluster.ids <- c("Anti_infl_Macs","B","Epithelial", "Mast","PC","Pro_infl_Macs","remove","Stromal_1",
                         "Stromal_2","Stromal_3","Stromal_4","Stromal_5","T")
new.cluster.ids <- c("Anti_infl_Macs","B","Epithelial", "Mast","PC","Pro_infl_Macs","remove","Enterocytes",
                     "remove","Endothelial","Fib_stricture","remove","T")
merged_cl_subCl$annotation <- plyr::mapvalues(x = merged_cl_subCl$sub.cluster, from = current.cluster.ids, to = new.cluster.ids)
DimPlot(merged_cl_subCl, reduction = "umap", label = TRUE, group.by = "annotation", label.size = 3)

###Epithelial 
Idents(merged_cl_subCl) <- "annotation"
merged_cl_subCl <- FindSubCluster(merged_cl_subCl, "Epithelial",graph.name = "SCT_snn", subcluster.name = "sub.cluster",resolution = 0.1, algorithm = 4)
DimPlot(merged_cl_subCl, reduction = "umap", label = TRUE, group.by = "sub.cluster")

Idents(merged_cl_subCl) <- "sub.cluster"
sub_celltype <- subset(merged_cl_subCl,idents = c( "Epithelial_1","Epithelial_2","Epithelial_3","Epithelial_4"))
DimPlot(sub_celltype, reduction = "umap", label = TRUE, group.by = "sub.cluster") 

markers <- FindAllMarkers(object = sub_celltype, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
View(markers %>% group_by(cluster) %>% top_n(n =10, wt = avg_log2FC))

DotPlot(sub_celltype, features = c(
  "LYVE1","COL3A1","COL1A1","GPM6B","S100B", #Stroma/Glia cluster 
  "PLVAP","CLDN5","SEMA3G", #Endothelial 
  "TFF1","TFF3","MUC2", #Goblet 
  "SPOCK1","COL4A1","LGALS1", #Fibroblasts 
  "EPCAM", #Epithelial cluster 
  "RGS5","COX4I2","IGFBP7", #Pericytes 
  "LGALS4","PIGR", #Transit amplyfier TA 
  "AQP8","PLAC8","GUCA2A","CEACAM5", #Enterocytes 
  "PLA2G2A","DEFA6","DEFA5","REG1A" #Paneth cells 
) )+   theme(axis.text.x = element_text(angle = 90)) 

#3 Goblet 
#1 Enterocytes 
#4 Paneth 
#remove 2

Idents(merged_cl_subCl) <- "sub.cluster"
current.cluster.ids <- c("Anti_infl_Macs","B","Endothelial","Enterocytes", "Epithelial_1", "Epithelial_2", "Epithelial_3", "Epithelial_4",
                         "Mast","PC","Pro_infl_Macs","remove","T")
new.cluster.ids <- c("Anti_infl_Macs","B","Endothelial","Enterocytes", "Enterocytes", "remove", "Goblet", "Paneth",
                     "Mast","PC","Pro_infl_Macs","remove","T")
merged_cl_subCl$annotation <- plyr::mapvalues(x = merged_cl_subCl$sub.cluster, from = current.cluster.ids, to = new.cluster.ids)
DimPlot(merged_cl_subCl, reduction = "umap", label = TRUE, group.by = "annotation", label.size = 3)

###T 
Idents(merged_cl_subCl) <- "annotation"
merged_cl_subCl <- FindSubCluster(merged_cl_subCl, "T",graph.name = "SCT_snn", subcluster.name = "sub.cluster",resolution = 0.1, algorithm = 4)
DimPlot(merged_cl_subCl, reduction = "umap", label = TRUE, group.by = "sub.cluster")

Idents(merged_cl_subCl) <- "sub.cluster"
sub_celltype <- subset(merged_cl_subCl,idents = c( "T_1","T_2","T_3"))
DimPlot(sub_celltype, reduction = "umap", label = TRUE, group.by = "sub.cluster") 

markers <- FindAllMarkers(object = sub_celltype, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
View(markers %>% group_by(cluster) %>% top_n(n =10, wt = avg_log2FC))
#remove 2 

DotPlot(sub_celltype, features = c(
  "CD3D", "CD2", "CD7", # T cells 
  "CCL5","CD8A",  #CD8 T 
  "CD40LG","IL2","FOXP3","ICOS", #CD4 T 
  "KLRC1","NKG7" #NK
))+   theme(axis.text.x = element_text(angle = 90)) 

#1 CD8
#3 CD4 

Idents(merged_cl_subCl) <- "sub.cluster"
current.cluster.ids <- c("Anti_infl_Macs","B","Endothelial","Enterocytes","Paneth", "Fib_stricture","Goblet",
                         "Mast","PC","Pro_infl_Macs","remove","T_1","T_2","T_3")
new.cluster.ids <- c("Anti_infl_Macs","B","Endothelial","Enterocytes","Paneth", "Fib_stricture","Goblet",
                     "Mast","PC","Pro_infl_Macs","remove","T_CD8","remove","T_CD4")
merged_cl_subCl$annotation <- plyr::mapvalues(x = merged_cl_subCl$sub.cluster, from = current.cluster.ids, to = new.cluster.ids)
DimPlot(merged_cl_subCl, reduction = "umap", label = TRUE, group.by = "annotation", label.size = 3)

###remove the "remove" clusters 
Idents(merged_cl_subCl) <- "annotation"
merged_cl <- subset(merged_cl_subCl, idents = c("Anti_infl_Macs","B","Endothelial","Enterocytes","Paneth", "Fib_stricture","Goblet",
                                                "Mast","PC","Pro_infl_Macs","T_CD8","T_CD4")) 
DimPlot(merged_cl, reduction = "umap", label = TRUE, group.by = "annotation", label.size = 3)

p <- DimPlot(merged_cl, label = TRUE,label.size = 5, group.by = "annotation", pt.size = 0.5, 
             cols = c("#117508","#E8E813","#850fb7","#915533",
                      "#db2aea","#b70f3f",
                      "#4f3603","#f44552", "#E88C13","#2bed1b","#5DB9E2",
                      "#464EE8")) + 
  theme(legend.title = element_text(size = 22), legend.text = element_text(size = 22)) +
  theme(title = element_text(size = 25))+ theme(axis.text = element_text(size = 30)) +
  ggtitle("Annotation") 
p + ggsave("./figures/Annotated_umap.pdf", width = 15, height = 10)
p + ggsave("./figures/Annotated_umap.svg", width = 15, height = 10)

##Colours
#Immune cells:  CD8T #464EE8, CD4T #5DB9E2  
#B #E8E813 ,PC  #E88C13 , Mast.#4f3603  
#Anti macs:#117508   Pro macs: #2bed1b
#Fib_stricture: #db2aea   Endothelial: #850fb7
#Enterocytes:#915533    Paneth: #f44552    Goblet: #b70f3f

p <- DotPlot(merged_cl, features = c(
  "CCL5","CD8A",  #CD8 T 
  "CD40LG","IL2","FOXP3","ICOS", #CD4 T 
  "TNFRSF17","MZB1", #Plasma cells 
  "BANK1","CD79B","CD22","MS4A1", #B cells 
  "FCN1","AIF1","S100A8","S100A9","IL1B","LYZ", #pro-inflammatory monocytes 
  "C1QC","C1QA","TREM2", #anti-inflammatory Macrophages 
  "TPSAB1","CMA1","KIT", #Mast cells 
  "PLVAP","CLDN5","SEMA3G", #Endothelial 
  "TFF1","TFF3","MUC2", #Goblet 
  "SPOCK1","COL4A1","LGALS1", #Fibroblasts stricture with S100A4
  "S100A4",
  "AQP8","PLAC8","GUCA2A","CEACAM5","EPCAM", #Enterocytes 
  "PLA2G2A","DEFA6","DEFA5","REG1A" #Paneth cells 
  
) )+   theme(axis.text.x = element_text(angle = 90))  + theme(title = element_text(size = 15))+ 
  theme(axis.text = element_text(size = 15))
p + ggsave("./figures/Dotplot.pdf", width = 15, height = 10)
p + ggsave("./figures/Dotplot.svg", width = 15, height = 10)

########## Remove negative and doublets from Seurat object ##########
spheres <- as.data.frame(table(merged_cl$sphere))$Var1
spheres <- as.character(spheres)
spheres <- spheres[!spheres %in% c("Negative_1","Negative_2",
                                   "Doublet_1","Doublet_2")]
Idents(merged_cl) <- "sphere"
merged_cl <- subset(merged_cl, idents = spheres)

#save R object 
saveRDS(merged_cl, "./data_files_generated/CrohnP393_P387.Rda")

###apply 5 cells per sphere cutoff 
Sphere_cell_cutoff("./data_files_generated/CrohnP393_P387.Rda",5,
                   "./data_files_generated/CrohnP393_P387_5cells.Rda")

#read in R object 
merged_cl <- readRDS("./data_files_generated/CrohnP393_P387_5cells.Rda")

p <- DimPlot(merged_cl, label = TRUE,label.size = 5, group.by = "annotation", pt.size = 0.5, 
             cols = c("#117508","#E8E813","#850fb7","#915533",
                      "#db2aea","#b70f3f",
                      "#4f3603","#f44552", "#E88C13","#2bed1b","#5DB9E2",
                      "#464EE8")) + NoLegend() +
  theme(legend.title = element_text(size = 22), legend.text = element_text(size = 22)) +
  theme(title = element_text(size = 25))+ theme(axis.text = element_text(size = 30)) +
  ggtitle("Annotation 5 cells per sphere") 
p + ggsave("./figures/Annotated_umap_5cells_noLegend.pdf", width = 15, height = 10)
p + ggsave("./figures/Annotated_umap_5cells_noLegend.svg", width = 15, height = 10)

##Colours
#Immune cells:  CD8T #464EE8, CD4T #5DB9E2  
#B #E8E813 ,PC  #E88C13 , Mast.#4f3603  
#Anti macs:#117508   Pro macs: #2bed1b
#Fib_stricture: #db2aea   Endothelial: #850fb7
#Enterocytes:#915533    Paneth: #f44552    Goblet: #b70f3f

########## Barplot of cell type proportions per sphere ##########
df_sphere_cell_type <- table(merged_cl$sphere,merged_cl$annotation)

#include Total column
df_sphere_cell_type <- cbind(df_sphere_cell_type, Total = rowSums(df_sphere_cell_type))

#cell type breakdown per sphere in percentage
df_sphere_cell_type <- as.data.frame(df_sphere_cell_type)

df_sphere_cell_type_pct = lapply(df_sphere_cell_type[,], function(x) {
  x/df_sphere_cell_type$Total
})
df_sphere_cell_type_pct <- as.data.frame(df_sphere_cell_type_pct)

#add rownames from sphere BC information 
rownames(df_sphere_cell_type_pct) <- rownames(df_sphere_cell_type)

#remove total column 
df_sphere_cell_type_pct <- df_sphere_cell_type_pct[,-ncol(df_sphere_cell_type_pct)]
df_sphere_cell_type_pct$sphere <- rownames(df_sphere_cell_type_pct)

#subset each celltype and add column name for cell type, then merge again 
ct1 <- df_sphere_cell_type_pct[,colnames(df_sphere_cell_type_pct) %in% c("Anti_infl_Macs","sphere")]
colnames(ct1) <- c("proportion","sphere")
ct1$annotation <- "Anti_infl_Macs"
ct2 <- df_sphere_cell_type_pct[,colnames(df_sphere_cell_type_pct) %in% c("B","sphere")]
colnames(ct2) <- c("proportion","sphere")
ct2$annotation <- "B"
ct3 <- df_sphere_cell_type_pct[,colnames(df_sphere_cell_type_pct) %in% c("Endothelial","sphere")]
colnames(ct3) <- c("proportion","sphere")
ct3$annotation <- "Endothelial"
ct4 <- df_sphere_cell_type_pct[,colnames(df_sphere_cell_type_pct) %in% c("Enterocytes","sphere")]
colnames(ct4) <- c("proportion","sphere")
ct4$annotation <- "Enterocytes"
ct5 <- df_sphere_cell_type_pct[,colnames(df_sphere_cell_type_pct) %in% c("Paneth","sphere")]
colnames(ct5) <- c("proportion","sphere")
ct5$annotation <- "Paneth"
ct6 <- df_sphere_cell_type_pct[,colnames(df_sphere_cell_type_pct) %in% c("Fib_stricture","sphere")]
colnames(ct6) <- c("proportion","sphere")
ct6$annotation <- "Fib_stricture"
ct7 <- df_sphere_cell_type_pct[,colnames(df_sphere_cell_type_pct) %in% c("Goblet","sphere")]
colnames(ct7) <- c("proportion","sphere")
ct7$annotation <- "Goblet"
ct8 <- df_sphere_cell_type_pct[,colnames(df_sphere_cell_type_pct) %in% c("Mast","sphere")]
colnames(ct8) <- c("proportion","sphere")
ct8$annotation <- "Mast"
ct9 <- df_sphere_cell_type_pct[,colnames(df_sphere_cell_type_pct) %in% c("PC","sphere")]
colnames(ct9) <- c("proportion","sphere")
ct9$annotation <- "PC"
ct10 <- df_sphere_cell_type_pct[,colnames(df_sphere_cell_type_pct) %in% c("Pro_infl_Macs","sphere")]
colnames(ct10) <- c("proportion","sphere")
ct10$annotation <- "Pro_infl_Macs"
ct11 <- df_sphere_cell_type_pct[,colnames(df_sphere_cell_type_pct) %in% c("T_CD4","sphere")]
colnames(ct11) <- c("proportion","sphere")
ct11$annotation <- "T_CD4"
ct12 <- df_sphere_cell_type_pct[,colnames(df_sphere_cell_type_pct) %in% c("T_CD8","sphere")]
colnames(ct12) <- c("proportion","sphere")
ct12$annotation <- "T_CD8"

#merge dataframes 
ct <- rbind(ct1,ct2)
ct <- rbind(ct,ct3)
ct <- rbind(ct,ct4)
ct <- rbind(ct,ct5)
ct <- rbind(ct,ct6)
ct <- rbind(ct,ct7)
ct <- rbind(ct,ct8)
ct <- rbind(ct,ct9)
ct <- rbind(ct,ct10)
ct <- rbind(ct,ct11)
ct <- rbind(ct,ct12)

#define order of spheres based on precence of CD4 T cells
df_sort <- table(merged_cl$sphere,merged_cl$annotation)
df_sort <- cbind(df_sort, Total = rowSums(df_sort))
df_sort <- as.data.frame(df_sort)
df_sort_pct = lapply(df_sort[,], function(x) {
  x/df_sort$Total
})
df_sort_pct <- as.data.frame(df_sort_pct)
rownames(df_sort_pct) <- rownames(df_sort)

#order decreasing base on CD4 T cells
df_sort_pct <- df_sort_pct[order(df_sort_pct$T_CD4, decreasing = TRUE), ]

ct$sphere <- factor(ct$sphere, levels = rownames(df_sort_pct))

p <- ggplot(ct, aes(fill=annotation, y=proportion, x=sphere)) + theme_classic() +
  geom_bar(position="stack", stat="identity" ) + 
  ggtitle("Cell type per sphere")+  theme(axis.text = element_text(size = 5)) + 
  theme(axis.title= element_text(size = 25))  +  theme(plot.title = element_text(size = 25, face = "bold")) + 
  theme(legend.title = element_text(size = 30), legend.text = element_text(size = 30))  + 
  xlab("Sphere") + ylab("Cell type proportion") + 
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_manual(values=  c("#117508","#E8E813","#850fb7","#915533",
                               "#db2aea","#b70f3f",
                               "#4f3603","#f44552", "#E88C13","#2bed1b","#5DB9E2",
                               "#464EE8"))
p + ggsave("./figures/Proportion.pdf", width = 15, height = 10)
p + ggsave("./figures/Proportion.svg", width = 15, height = 10)

########## Add inflamed and non-ninflamed Phenotype ##########
#P393 _1 non-inflamed = Bar97-Bar192, inflamed = Bar1-Bar96, Bar193-288
#P387 _2 non-inflamed = Bar145-Bar192, inflamed = Bar1-144, 193-288
spheres_infl1 <- paste0("Bar",c(1:96,193:288),"_1")
spheres_noninfl1 <- paste0("Bar",c(97:192),"_1")
spheres_infl2 <- paste0("Bar",c(1:144,193:288),"_2")
spheres_noninfl2 <- paste0("Bar",c(145:192),"_2")

spheres_infl <- c(spheres_infl1, spheres_infl2)
spheres_noninfl <- c(spheres_noninfl1, spheres_noninfl2)

merged_cl$condition <- NA
merged_cl@meta.data <- merged_cl@meta.data %>%
  mutate(condition = case_when(
    sphere %in% spheres_infl ~ "inflamed",
    sphere %in% spheres_noninfl ~ "nonInflamed",
    TRUE ~ NA_character_))

p <- DimPlot(merged_cl, label = TRUE,label.size = 5, group.by = "annotation", pt.size = 0.5, split.by = "condition",
             cols = c("#117508","#E8E813","#850fb7","#915533",
                      "#db2aea","#b70f3f",
                      "#4f3603","#f44552", "#E88C13","#2bed1b","#5DB9E2",
                      "#464EE8")) + 
  theme(legend.title = element_text(size = 22), legend.text = element_text(size = 22)) +
  theme(title = element_text(size = 25))+ theme(axis.text = element_text(size = 30)) +
  ggtitle("Annotation 5 cells per sphere") 
p + ggsave("./figures/umap_infl_non_infl.pdf", width = 15, height = 10)
p + ggsave("./figures/umap_infl_non_infl.svg", width = 15, height = 10)

#save R object 
saveRDS(merged_cl, "./data_files_generated/CrohnP393_P387_5cell_annotated_infl_noninfl.Rda")


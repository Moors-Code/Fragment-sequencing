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

########## merge and cluster ##########
###add prefix to sphere for sample ID
CrohnP393$sphere <- paste(CrohnP393$sphere, "_1", sep="")
CrohnP387$sphere <- paste(CrohnP387$sphere, "_2", sep="")

###merge
merged <- merge(CrohnP393, CrohnP387)

###Remove negative and doublets from Seurat object 
spheres <- as.data.frame(table(merged$sphere))$Var1
spheres <- as.character(spheres)
spheres <- spheres[!spheres %in% c("Negative_1","Negative_2",
                                   "Doublet_1","Doublet_2")]
Idents(merged) <- "sphere"
merged <- subset(merged, idents = spheres)

###Remove low quality cells
#add percentage of mitochondrial features
mito.features <- grep(pattern = "^MT-", x = rownames(x = merged), value = TRUE)
percent.mito <- Matrix::colSums(x = GetAssayData(object = merged, slot = 'counts')[mito.features, ]) / Matrix::colSums(x = GetAssayData(object = merged, slot = 'counts'))
merged <- AddMetaData(object = merged, metadata = percent.mito, col.name = "percent.mito")
Idents(merged) <- "orig.ident"
FeatureScatter(merged, feature1 = "nCount_RNA", feature2 = "percent.mito")
VlnPlot(merged, features = "percent.mito") 
FeatureScatter(merged, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

#remove low quality cells 
merged <- subset(merged, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mito < 0.25 )

#save R object 
saveRDS(merged, "./data_files_generated/CrohnP393_P387.Rda")

###apply 5 cells per sphere cutoff 
Sphere_cell_cutoff("./data_files_generated/CrohnP393_P387.Rda",5,
                   "./data_files_generated/CrohnP393_P387_5cells.Rda")

###clustering 
merged <- readRDS("./data_files_generated/CrohnP393_P387_5cells.Rda")

merged <- SCTransform(merged, vars.to.regress = c("nCount_RNA","percent.mito"), verbose = FALSE)
merged <- RunPCA(object = merged, features = VariableFeatures(object = merged), npcs = 30, verbose = FALSE)
#check for batch effect  - no batch effect 
DimPlot(merged, reduction = "pca", group.by = "orig.ident")

#check for significant PCAs
ElbowPlot(merged)

#cluster
merged <- FindNeighbors(object = merged, dims = 1:10, reduction = "pca")
merged_cl <- FindClusters(merged, resolution = 0.5, random.seed = 5, algorithm = 1, graph.name = "SCT_snn")
merged_cl <- RunUMAP(merged_cl, dims = 1:10, seed.use = 5, reduction = "pca")
DimPlot(merged_cl, reduction = "umap", label = TRUE, pt.size = .5, group.by = "seurat_clusters") + ggsave("./figures/5/umap_clustering.pdf")

########## Annotation based on marker gene from Martin et al 2019, Cell ##########
DotPlot(merged_cl, features = c(
  "CD3D", "CD2", "CD7", # T cells cluster 0,1,2,12
  "TNFRSF17","MZB1", #Plasma cells cluster 6,7,8,9,10,15
  "BANK1","CD79B","CD22","MS4A1", #B cells cluster 5
  "HLA-DRB1","HLA-DQA1", "LYZ", #MNPs Mononuclear phagocyte system (monocytes, macrophages ) cluster 3,4
  "GZMB","LILRA4", #pDCs 
  "TPSAB1","CMA1","KIT", #Mast cells cluser 14
  "LYVE1","COL3A1","COL1A1","GPM6B","S100B", #Stroma/Glia cluster 11
  "EPCAM" #Epithelial cluster 13
) )+   theme(axis.text.x = element_text(angle = 90)) + ggsave("./figures/5/marker_clustering_dotplot.pdf")

current.cluster.ids <- c(0:15)
new.cluster.ids <- c("T","T","T","MNPs","MNPs","B","PC","PC",
                     "PC","PC","PC","Stromal","T", "Epithelial","Mast","PC")
merged_cl$annotation <- plyr::mapvalues(x = merged_cl$seurat_clusters, from = current.cluster.ids, to = new.cluster.ids)

DimPlot(merged_cl, group.by = "annotation",label = TRUE, cols = c("#2323E5","#19E80F","#EF975B","#F4741E",
                                                                  "#8E5229","#D2E51C","#E20FE8")) + ggsave("./figures/5/annotated_umap.pdf")

##Colours
#Immune cells:  T #2323E5,  B#EF975B ,PC  #F4741E , Mast.#E20FE8   MNPs #19E80F
#Stromal:  #8E5229
#Epithelial:  yellow #D2E51C

##### Barplot of cell type proportions #####
df_sphere_cell_type <- table(merged_cl$sphere,merged_cl$annotation)

#include Total column
df_sphere_cell_type <- cbind(df_sphere_cell_type, Total = rowSums(df_sphere_cell_type))

#cell type breakdown per sphere in percentage, generates at data frame with percentages  
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
ct1 <- df_sphere_cell_type_pct[,colnames(df_sphere_cell_type_pct) %in% c("B","sphere")]
colnames(ct1) <- c("proportion","sphere")
ct1$annotation <- "B"
ct2 <- df_sphere_cell_type_pct[,colnames(df_sphere_cell_type_pct) %in% c("T","sphere")]
colnames(ct2) <- c("proportion","sphere")
ct2$annotation <- "T"
ct3 <- df_sphere_cell_type_pct[,colnames(df_sphere_cell_type_pct) %in% c("MNPs","sphere")]
colnames(ct3) <- c("proportion","sphere")
ct3$annotation <- "MNPs"
ct4 <- df_sphere_cell_type_pct[,colnames(df_sphere_cell_type_pct) %in% c("PC","sphere")]
colnames(ct4) <- c("proportion","sphere")
ct4$annotation <- "PC"
ct5 <- df_sphere_cell_type_pct[,colnames(df_sphere_cell_type_pct) %in% c("Epithelial","sphere")]
colnames(ct5) <- c("proportion","sphere")
ct5$annotation <- "Epithelial"
ct6 <- df_sphere_cell_type_pct[,colnames(df_sphere_cell_type_pct) %in% c("Stromal","sphere")]
colnames(ct6) <- c("proportion","sphere")
ct6$annotation <- "Stromal"
ct7 <- df_sphere_cell_type_pct[,colnames(df_sphere_cell_type_pct) %in% c("Mast","sphere")]
colnames(ct7) <- c("proportion","sphere")
ct7$annotation <- "Mast"

#merge dataframes 
ct <- rbind(ct1,ct2)
ct <- rbind(ct,ct3)
ct <- rbind(ct,ct4)
ct <- rbind(ct,ct5)
ct <- rbind(ct,ct6)
ct <- rbind(ct,ct7)

#define order of spheres based on presense of MNPs
df_zone <- table(merged_cl$sphere,merged_cl$annotation)

#include Total column
df_zone <- cbind(df_zone, Total = rowSums(df_zone))

#cell type breakdown per sphere in percentage, generates at data frame with percentages  
df_zone <- as.data.frame(df_zone)

df_zone_pct = lapply(df_zone[,], function(x) {
  x/df_zone$Total
})
df_zone_pct <- as.data.frame(df_zone_pct)

#add rownames from sphere BC information 
rownames(df_zone_pct) <- rownames(df_zone)

#order decreasing base on T zone 
df_zone_pct <- df_zone_pct[order(df_zone_pct$T, decreasing = TRUE), ]

ct$sphere <- factor(ct$sphere, levels = rownames(df_zone_pct))

##Colours
#Immune cells:  T #2323E5,  B#EF975B ,PC  #F4741E , Mast.#E20FE8   MNPs #19E80F
#Stromal:  #8E5229
#Epithelial:  yellow #D2E51C

p <- ggplot(ct, aes(fill=annotation, y=proportion, x=sphere)) + theme_classic() +
  geom_bar(position="stack", stat="identity" ) + 
  ggtitle("Cell type per sphere")+  theme(axis.text = element_text(size = 5)) + 
  theme(axis.title= element_text(size = 25))  +  theme(plot.title = element_text(size = 25, face = "bold")) + 
  theme(legend.title = element_text(size = 30), legend.text = element_text(size = 30))  + 
  xlab("Sphere") + ylab("Cell type proportion") + 
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_manual(values=  c("#EF975B","#D2E51C","#E20FE8","#19E80F","#F4741E","#8E5229","#2323E5"))
p + ggsave("./figures/5/Proportion_per_Crohn_legend.pdf", width = 15, height = 10)
p + ggsave("./figures/5/Proportion_per_Crohn_legend.svg", width = 15, height = 10)


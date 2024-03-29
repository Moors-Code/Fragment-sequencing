########## Part 4: Spleen preliminary analysis ##########
#This part analyses preliminary spleen fragment-seq data  

########## Prepare environment ##########
###Setting the working directory 
setwd("/mnt/khandler/R_projects/Fragment-sequencing/Fragment-seq_analysis/")

###Load packages and functions 
source("./functions_and_packages/1.Packages.R")
source("./functions_and_packages/4.Functions_cells_per_fragment_cutoff.R")

###Load data 
spleenS1 <- readRDS("./data_files_generated/SpS_Spleen_Sample1.Rda")
spleenS2 <- readRDS("./data_files_generated/SpS_Spleen_Sample2.Rda")

########## merge samples and cluster ##########
###add prefix to fragment for sample ID
spleenS1$fragment <- paste(spleenS1$fragment, "_1", sep="")
spleenS2$fragment <- paste(spleenS2$fragment, "_2", sep="")

###merge
merged <- merge(spleenS1, spleenS2)

###Remove negative and doublets from Seurat object 
fragments <- as.data.frame(table(merged$fragment))$Var1
fragments <- as.character(fragments)
fragments <- fragments[!fragments %in% c("Negative_1","Negative_2",
                                   "Doublet_1","Doublet_2")]
Idents(merged) <- "fragment"
merged <- subset(merged, idents = fragments)

###Remove low quality cells
#add percentage of mitochondrial features
mito.features <- grep(pattern = "^mt-", x = rownames(x = merged), value = TRUE)
percent.mito <- Matrix::colSums(x = GetAssayData(object = merged, slot = 'counts')[mito.features, ]) / Matrix::colSums(x = GetAssayData(object = merged, slot = 'counts'))
merged <- AddMetaData(object = merged, metadata = percent.mito, col.name = "percent.mito")
Idents(merged) <- "orig.ident"
FeatureScatter(merged, feature1 = "nCount_RNA", feature2 = "percent.mito")
VlnPlot(merged, features = "percent.mito") 
FeatureScatter(merged, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

#remove low quality cells
merged <- subset(merged, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mito < 0.1 )

#save R object 
saveRDS(merged,"./data_files_generated/Spleen_merged.Rda")

###apply 5 cells per fragment cutoff 
Fragment_cell_cutoff("./data_files_generated/Spleen_merged.Rda",5,
                   "./data_files_generated/Spleen_merged_5cells.Rda")

###clustering 
merged <- readRDS("./data_files_generated/Spleen_merged_5cells.Rda")

merged <- SCTransform(merged, vars.to.regress = c("nCount_RNA","percent.mito"), verbose = FALSE)
merged <- RunPCA(object = merged, features = VariableFeatures(object = merged), npcs = 30, verbose = FALSE)

#check for batch effect  - no batch effect 
DimPlot(merged, reduction = "pca", group.by = "orig.ident")

#check for significant PCAs
ElbowPlot(merged)

#cluster 
merged <- FindNeighbors(object = merged, dims = 1:10, reduction = "pca")
merged_cl <- FindClusters(merged, resolution = 0.5, random.seed = 5, algorithm = 4, graph.name = "SCT_snn")
merged_cl <- RunUMAP(merged_cl, dims = 1:10, seed.use = 5, reduction = "pca")
DimPlot(merged_cl, reduction = "umap", label = TRUE, pt.size = .5, group.by = "seurat_clusters") + ggsave("./figures/4/umap_clustering.pdf")

########## Annotation based on marker genes from Medaglia et al 2017 Science ##########
DotPlot(merged_cl, features = c(
  "Cd79b", "Ly6d", # B cluster 1,4,5,8
  "Cr2", #B Cd21 high  cluster 3
  "Sdc1", #Plasma cells PC 
  "Trbc2","Tcf7", "Igfbp4", #CD4 T cluster 7
  "Cd8a", #CD8 T  cluster 2,10,17
  "Gzma","Ccl5", #NK cluster 9
  "S100a4", #Mono Ly6c low cluster 6
  "Ccr2","Lyz2", #Mono Ly6c high cluster 12
  "Cd5l","Vcam1", #Mac cluster 13
  "Bst2", #pDC cluster 16
  "Fscn1", #DC cluster none
  "S100a8", #Neut cluster 15
  "Hba-a2" #RBC cluster 11
) )+   theme(axis.text.x = element_text(angle = 90)) + ggsave("./figures/4/marker_clustering_dotplot.pdf")

current.cluster.ids <- c(1:17)
new.cluster.ids <- c("B","CD8_T","B_CD21_high","B","B","Mono_Ly6c_low","CD4_T","B","NK",
                     "CD8_T","RBC", "Mono_Ly6c_high","Macrophages","?","Neutrophils","pDC","CD8_T")
merged_cl$annotation <- plyr::mapvalues(x = merged_cl$seurat_clusters, from = current.cluster.ids, to = new.cluster.ids)

#remove RBCs 
Idents(merged_cl) <- "annotation"
merged_cl <- subset(merged_cl, idents = c("B","B_CD21_high", "CD8_T","CD4_T","Mono_Ly6c_low","NK","Mono_Ly6c_high","Macrophages","Neutrophils","pDC"))

#T cell areas: red colours:  CD4 T , CD8 T , NK 
#B cell areas: Blue colours:  B, B CD21 high, PC 
#Marginal zones MZ: purple colours: Mono Ly6c low, mono ly6c high, Mac, pDC, Neut

current.cluster.ids <- c("B","B_CD21_high", "CD8_T","CD4_T","Mono_Ly6c_low","NK","Mono_Ly6c_high","Macrophages","Neutrophils","pDC")
new.cluster.ids <- c("B","B", "T","T","MZ","T","MZ","MZ","MZ","MZ")
merged_cl$zone <- plyr::mapvalues(x = merged_cl$annotation, from = current.cluster.ids, to = new.cluster.ids)

#remove RBCs
Idents(merged_cl) <- "zone"
merged_cl <- subset(merged_cl, idents = c("B","T","MZ"))

#add zone into meta data based on cell type presence 
DimPlot(merged_cl, group.by = "zone",label = TRUE) + ggsave("./figures/4/zone_umap.pdf")
p <- DimPlot(merged_cl, group.by = "annotation",label = TRUE, cols = c("#1C3EE5","#EA32EA","#1CA7E5", "#D2E51C",  "#B50926",
                                                                       "#F9AFC0","#E5A21C",
                                                                       "#EFCD8B","#F2EB08","#EFEF8B" )) 
p + ggsave("./figures/4/annotated_umap.pdf")
p + ggsave("./figures/4/annotated_umap.svg")

###Colours
#T cell areas: red colors:  CD4 T #B50926, CD8 T#EA32EA , NK #F9AFC0
#B cell areas: Blue colors:  B #1C3EE5, B CD21 high #1CA7E5
#Marginal zones MZ: purple colors: Mono Ly6c low #D2E51C, mono ly6c high #E5A21C, Mac #EFCD8B, pDC #EFEF8B, Neut #F2EB08

########## Barplot of cell type proportions ##########
df_fragment_cell_type <- table(merged_cl$fragment,merged_cl$annotation)

#include Total column
df_fragment_cell_type <- cbind(df_fragment_cell_type, Total = rowSums(df_fragment_cell_type))

#cell type breakdown per fragment in percentage, generates at data frame with percentages  
df_fragment_cell_type <- as.data.frame(df_fragment_cell_type)

df_fragment_cell_type_pct = lapply(df_fragment_cell_type[,], function(x) {
  x/df_fragment_cell_type$Total
})
df_fragment_cell_type_pct <- as.data.frame(df_fragment_cell_type_pct)

#add rownames from fragment BC information 
rownames(df_fragment_cell_type_pct) <- rownames(df_fragment_cell_type)

#remove total column 
df_fragment_cell_type_pct <- df_fragment_cell_type_pct[,-ncol(df_fragment_cell_type_pct)]
df_fragment_cell_type_pct$fragment <- rownames(df_fragment_cell_type_pct)
#subset each celltype and add column name for cell type, then merge again 
ct1 <- df_fragment_cell_type_pct[,colnames(df_fragment_cell_type_pct) %in% c("B","fragment")]
colnames(ct1) <- c("proportion","fragment")
ct1$annotation <- "B"
ct2 <- df_fragment_cell_type_pct[,colnames(df_fragment_cell_type_pct) %in% c("CD8_T","fragment")]
colnames(ct2) <- c("proportion","fragment")
ct2$annotation <- "CD8_T"
ct3 <- df_fragment_cell_type_pct[,colnames(df_fragment_cell_type_pct) %in% c("B_CD21_high","fragment")]
colnames(ct3) <- c("proportion","fragment")
ct3$annotation <- "B_CD21_high"
ct4 <- df_fragment_cell_type_pct[,colnames(df_fragment_cell_type_pct) %in% c("Mono_Ly6c_low","fragment")]
colnames(ct4) <- c("proportion","fragment")
ct4$annotation <- "Mono_Ly6c_low"
ct5 <- df_fragment_cell_type_pct[,colnames(df_fragment_cell_type_pct) %in% c("CD4_T","fragment")]
colnames(ct5) <- c("proportion","fragment")
ct5$annotation <- "CD4_T"
ct6 <- df_fragment_cell_type_pct[,colnames(df_fragment_cell_type_pct) %in% c("NK","fragment")]
colnames(ct6) <- c("proportion","fragment")
ct6$annotation <- "NK"
ct7 <- df_fragment_cell_type_pct[,colnames(df_fragment_cell_type_pct) %in% c("Mono_Ly6c_high","fragment")]
colnames(ct7) <- c("proportion","fragment")
ct7$annotation <- "Mono_Ly6c_high"
ct8 <- df_fragment_cell_type_pct[,colnames(df_fragment_cell_type_pct) %in% c("Macrophages","fragment")]
colnames(ct8) <- c("proportion","fragment")
ct8$annotation <- "Macrophages"
ct9 <- df_fragment_cell_type_pct[,colnames(df_fragment_cell_type_pct) %in% c("Neutrophils","fragment")]
colnames(ct9) <- c("proportion","fragment")
ct9$annotation <- "Neutrophils"
ct10 <- df_fragment_cell_type_pct[,colnames(df_fragment_cell_type_pct) %in% c("pDC","fragment")]
colnames(ct10) <- c("proportion","fragment")
ct10$annotation <- "pDC"

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

#define order of fragments based on T cell zone percentage high to low 
df_zone <- table(merged_cl$fragment,merged_cl$zone)

#include Total column
df_zone <- cbind(df_zone, Total = rowSums(df_zone))

#cell type breakdown per fragment in percentage, generates at data frame with percentages  
df_zone <- as.data.frame(df_zone)

df_zone_pct = lapply(df_zone[,], function(x) {
  x/df_zone$Total
})
df_zone_pct <- as.data.frame(df_zone_pct)

#add rownames from fragment BC information 
rownames(df_zone_pct) <- rownames(df_zone)

#order decreasing base on T zone 
df_zone_pct <- df_zone_pct[order(df_zone_pct$B, decreasing = TRUE), ]

ct$fragment <- factor(ct$fragment, levels = rownames(df_zone_pct))

###Colours
#T cell areas: red colors:  CD4 T #B50926, CD8 T#EA32EA , NK #F9AFC0
#B cell areas: Blue colors:  B #1C3EE5, B CD21 high #1CA7E5
#Marginal zones MZ: purple colors: Mono Ly6c low #D2E51C, mono ly6c high #E5A21C, Mac #EFCD8B, pDC #EFEF8B, Neut #F2EB08

#reorder cell types for plotting 
ct$annotation <- factor(ct$annotation,levels = c("B","B_CD21_high","CD4_T","CD8_T","NK", "Macrophages","Mono_Ly6c_high",
                                                 "Mono_Ly6c_low","Neutrophils","pDC"))
p <- ggplot(ct, aes(fill=annotation, y=proportion, x=fragment)) + theme_classic() +
  geom_bar(position="stack", stat="identity" ) + 
  ggtitle("Cell type per fragment")+  theme(axis.text = element_text(size = 5)) + 
  theme(axis.title= element_text(size = 25))  +  theme(plot.title = element_text(size = 25, face = "bold")) + 
  theme(legend.title = element_text(size = 30), legend.text = element_text(size = 30))  + 
  xlab("Fragment") + ylab("Cell type proportion") + 
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_manual(values=  c("#1C3EE5","#1CA7E5", "#B50926","#EA32EA", "#F9AFC0", "#EFCD8B","#E5A21C","#D2E51C","#F2EB08",
                               "#EFEF8B"))
p + ggsave("./figures/4/Proportion_per_fragment_legend.pdf", width = 15, height = 10)
p + ggsave("./figures/4/Proportion_per_fragment_legend.svg", width = 15, height = 10)

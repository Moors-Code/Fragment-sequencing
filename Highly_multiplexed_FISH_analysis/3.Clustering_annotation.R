########## Part 3: Clustering and annotation ##########
#This part does clustering and annotation on highly multiplexed FISH data 

########## Prepare environment ##########
###Setting the working directory 
setwd("/mnt/khandler/R_projects/Sphere-sequencing/Highly_multiplexed_FISH_analysis/")

###Load packages and functions 
source("./functions_and_packages/1.Packages.R")
source("./functions_and_packages/2.Functions_ImageJ_connection.R")

###Load data 
#mets_samples = with visible metastasis, noMets_samples = no visible metastasis
mets_samples <- readRDS(file = "./data_files_generated/mets_QC_imageJcoordinate.rds")
noMets_samples <- readRDS(file = "./data_files_generated/noMets_QC_imageJcoordinate.rds")

mets_samples$samples <- "mets"
noMets_samples$samples <- "no_mets"

########## Merging and clustering ##########
###merge samples 
merged <- merge(mets_samples, noMets_samples)
###cluster cells 
merged <- FindVariableFeatures(merged, selction.method = "vst", nfeatures = 100)
all.genes <- rownames(merged)
merged <- ScaleData(merged,features = all.genes)
merged <- RunPCA(object = merged, features = VariableFeatures(object = merged),approx=FALSE)
ElbowPlot(merged)
merged <- FindNeighbors(object = merged, dims = 1:10, reduction = "pca")
merged <- FindClusters(merged, resolution = 0.2, random.seed = 2, algorithm = 1, graph.name = "originalexp_snn")
merged <- RunUMAP(merged, dims = 1:10, seed.use = 5, reduction = "pca")
DimPlot(merged, reduction = "umap", label = TRUE)

###remove all the little random cluster, these are probably low quality cells only expression a view molecules 
Idents(merged) <- "seurat_clusters"
merged <- subset(merged, idents = c(0,1,2,3,4,5))
DimPlot(merged, reduction = "umap", label = TRUE) 

###cluster again 
merged <- FindNeighbors(object = merged, dims = 1:10, reduction = "pca")
merged <- FindClusters(merged, resolution = 1, random.seed = 2, algorithm = 1, graph.name = "originalexp_snn") 
merged <- RunUMAP(merged, dims = 1:10, seed.use = 5, reduction = "pca")
DimPlot(merged, reduction = "umap", label = TRUE) 

########## Annotation ##########
###first broad annotation, for clusters that are not clear do DGE and compare top genes with Sphere-seq analysis 
marker_genes <- c("Pglyrp1","Gpx2", #Metastasis 1,4
                  "Cyp2e1","Cyp1a2" , #Hepatocytes CV 2,3,14
                  "Fgb","Cyp2f2",  #Hepatocytes PV 0,7,18
                  "Clec4f","Vsig4", #Kupffer 8
                  "Lyz2","C1qc",#Monocytes 6
                  "Csf3r", #Neutrophils (also 6, but difficult to say with only one marker, annotate as Monocytes)
                  "Galnt15","Acer2","Pecam1",#LECs 5,9,15
                  "App","Spp1",#Cholangiocytes 
                  "Plvap","Lrat", #Stellate 10
                  "Il2rb","Ccl5","Cd3d", #T 11
                  "Cald1","Fn1", #Fibroblasts
                  "Ighm","Itga4" #B 11
)
DotPlot(merged,features = marker_genes) + theme(axis.text.x = element_text(angle = 90)) 

#non conclusive 12,13,16,17,19,20,21,22,23,24
#recluster T and B cluster 11

###check top DEGs from clusters that are non-conclusive and compare top genes with sphere-seq data 
dge_genes <- FindAllMarkers(merged, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
View(dge_genes %>% group_by(cluster) %>% top_n(n =10, wt = avg_log2FC))

#load sphere-seq data 
liverSpS5C <- readRDS("/mnt/khandler/R_projects/Sphere-sequencing/Sphere-seq_analysis/data_files_generated/LiverMerged_afterBC_anno.Rda")
Idents(liverSpS5C) <- "annotation.broad"

dge_genes <- FindAllMarkers(merged, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
View(dge_genes %>% group_by(cluster) %>% top_n(n =10, wt = avg_log2FC))

#Clusters 12 --> non conclusive --> remove
DotPlot(liverSpS5C, features =c("App","Pglyrp1","Dhrs3","Nrp1","Spp1"))

#Clusters 13 --> Fibroblasts
DotPlot(liverSpS5C, features =c("Fn1","Rrbp1","Cald1","App","Spp1"))

#Clusters 16 --> Kupffer/LECs doublets --> remove 
DotPlot(liverSpS5C, features =c("Cd5l","Clec4f","Vsig4","C1qc","Cyp2e1"))

#Clusters 17 --> non conclusive --> remove
DotPlot(liverSpS5C, features =c("Pglyrp1","Ltbr","App","Sema4a","Psen1"))

#Clusters 19 --> non conclusive --> remove
DotPlot(liverSpS5C, features =c("Il2rb","Sema4d","Gprc5c","Fn1","Pglyrp1"))

#Clusters 20 --> non conclusive --> remove
DotPlot(liverSpS5C, features =c("Ltbr","Tgfbi","Rrbp1","App","Dab2"))

#Clusters 21 --> non conclusive --> remove
DotPlot(liverSpS5C, features =c("Lhx6","Spp1","Pglyrp1","App","Fn1","Jup"))

#cluster 22,23,24 only express one gene, remove from the analysis

#Clusters 22
DotPlot(liverSpS5C, features =c("Pglyrp1","Gpx2","App","Spp1","Mecom","Jup"))


current.cluster.ids <- c(0:24)
new.cluster.ids <- c("Hepatocytes_PV","Metastasis","Hepatocytes_CV","Hepatocytes_CV","Metastasis","LECs","Monocytes",
                     "Hepatocytes_PV","Kupffer","LECs","Stellate","T_B","Fibroblasts",
                     "remove","Hepatocytes_CV","LECs","remove","remove",
                     "Hepatocytes_PV","remove","remove","remove","remove","remove",
                     "remove")
merged$annotation <- plyr::mapvalues(x = merged$seurat_clusters, from = current.cluster.ids, to = new.cluster.ids)
DimPlot(merged, reduction = "umap", label = TRUE, group.by = "annotation", label.size = 3)

###subcluster B_T cell cluster 
Idents(merged) <- "annotation"
merged <- FindSubCluster(merged, "T_B",graph.name = "originalexp_snn", subcluster.name = "sub.cluster",resolution = 0.2, algorithm = 4)
DimPlot(merged, reduction = "umap", label = TRUE, group.by = "sub.cluster")

Idents(merged) <- "sub.cluster"
sub_celltype <- subset(merged,idents = c( "T_B_1","T_B_2","T_B_3","T_B_4"))
DimPlot(sub_celltype, reduction = "umap", label = TRUE, group.by = "sub.cluster")

Idents(sub_celltype) <- "sub.cluster"
marker_genes <- c("Il2rb","Ccl5","Cd3d", #T 1,3
                  "Ighm","Itga4" #B 2,4
) 
DotPlot(sub_celltype,features = marker_genes) + theme(axis.text.x = element_text(angle = 90)) 

current.cluster.ids <- c("Fibroblasts","Hepatocytes_CV","Hepatocytes_PV","Kupffer","LECs","Metastasis","Monocytes","remove","Stellate",
                         "T_B_1","T_B_2","T_B_3","T_B_4")
new.cluster.ids <- c("Fibroblasts","Hepatocytes_CV","Hepatocytes_PV","Kupffer","LECs","Metastasis","Monocytes","remove","Stellate",
                     "T","B","T","B")
merged$annotation <- plyr::mapvalues(x = merged$sub.cluster, from = current.cluster.ids, to = new.cluster.ids)
DimPlot(merged, reduction = "umap", label = TRUE, group.by = "annotation", label.size = 3)

#remove Remove cluster 
merged <- subset(merged, idents = c("Fibroblasts","Hepatocytes_CV","Hepatocytes_PV","Kupffer","LECs","Metastasis","Monocytes","Stellate",
                                    "T","B"))

###
#DotPlot of marker 
marker_genes <- c("Pglyrp1","Gpx2", #Metastasis
                  "Cyp2e1","Cyp1a2" , #Hepatocytes CV 
                  "Fgb","Cyp2f2",  #Hepatocytes PV
                  "Clec4f","Vsig4", #Kupffer
                  "Lyz2","Csf3r","C1qc",#Monocytes, Neutrophils Csf3r top  
                  "Galnt15","Acer2","Pecam1",#LECs
                  "App","Spp1",#Cholangiocytes
                  "Tgfbi","Plvap","Lrat", #Stellate
                  "Il2rb","Ccl5","Cd3d", #T
                  "Cald1","Fn1", #Fibroblasts
                  "Ighm","Itga4" #B 
)

Idents(merged) <- "annotation"
#define order of cell types 
p <- DotPlot(merged,features = marker_genes) + theme(axis.text.x = element_text(angle = 90)) + 
  ggtitle("Marker genes in annotation of Resolve data") + 
  theme(legend.title = element_text(size = 22), legend.text = element_text(size = 22)) + 
  theme(axis.title= element_text(size = 25)) + theme(axis.text = element_text(size = 15))+ 
  theme(plot.title = element_text(size = 25, face = "bold")) 
p+ ggsave("./figures/3/marker_dotPlot.svg",width = 12, height = 10)
p+ ggsave("./figures/3/marker_dotPlot.pdf",width = 12, height = 10)


#Annotated umap plots 
p <- DimPlot(merged, reduction = "umap", label = TRUE, group.by = "annotation", label.size = 3, cols = 
               c("#EF975B","#8E5229","#F46042","#56040C",
                 "#E20FE8","#F4E740","#F90606","#19E80F" ,"#2323E5"), pt.size = 0.5) + 
  theme(legend.title = element_text(size = 22), legend.text = element_text(size = 22)) + 
  theme(title = element_text(size = 25))+ theme(axis.text = element_text(size = 30)) 
p + ggsave("./figures/3/annotated_umap.pdf",width = 15, height = 10)
p + ggsave("./figures/3/annotated_umap.svg",width = 15, height = 10)

#colours: 
#Hepatocytes (red): CV #F46042 PV #56040C
#KC: #E20FE8   
#LECs:  #F4E740
#Metastasis (pink) #F90606
#Stellate #C1A08A
#Fibroblasts #8E5229 
#Monocytes #19E80F 
#T #2323E5 
#B #EF975B

###change identity to have a vein meta data column and then combine sm and nm to distal and mets to proximal in Mets_distance
current.cluster.ids <- c("cv","cv_nm","cv_sm","mets","pv","pv_nm","pv_sm")
new.cluster.ids <- c("CV","CV","CV","Mets","PV","PV","PV")
merged$vein <- plyr::mapvalues(x = merged$spatial_feature, from = current.cluster.ids, to = new.cluster.ids)

current.cluster.ids <- c("cv","cv_nm","cv_sm","mets","pv","pv_nm","pv_sm")
new.cluster.ids <- c("distal","distal","distal","proximal","distal","distal","distal")
merged$Mets_distance <- plyr::mapvalues(x = merged$spatial_feature, from = current.cluster.ids, to = new.cluster.ids)

########## check CV and PV heps in veins you draw on ImageJ ##########
Idents(merged) <- "annotation"
hep <- subset(merged, idents = c("Hepatocytes_CV","Hepatocytes_PV"))
p <- DimPlot(hep,group.by = "annotation", split.by = "vein", cols = c("#56040C","#F46042"), pt.size = 0.5) + theme(legend.title = element_text(size = 22), legend.text = element_text(size = 22)) + 
  theme(title = element_text(size = 25))+ theme(axis.text = element_text(size = 30)) 
p + ggsave("./figures/3/split_plot_pv_cv_hep.pdf",width = 15, height = 10)
p + ggsave("./figures/3/split_plot_pv_cv_hep.svg",width = 15, height = 10)

#change to sample ids 
current.cluster.ids <- c("Slide1_A1-1" ,"Slide1_A2-1"   ,"Slide1_A2-2"     ,"Slide1_B1-1" ,"Slide1_B1-2"  ,"Slide1_B2-1" )
new.cluster.ids <- c("metsA1","no_metsA2","no_metsA2","metsB1","metsB1","no_metsB2")
merged$sampleID <- plyr::mapvalues(x = merged$Slide, from = current.cluster.ids, to = new.cluster.ids)

########## save R object ##########
saveRDS(merged,file = "./data_files_generated/Resolve_seurat_anno.rds")

########## Get cell annotation per slide and celltype ##########
setwd("./data_files_generated/annotation_file_for_ImageJ")
Idents(merged) <- "Slide"

cell_types_mets <- c("Metastasis","Hepatocytes_CV","Kupffer","Monocytes","Hepatocytes_PV","LECs",
                     "Stellate","T","Fibroblasts","B")
cell_types_noMets <- c("Hepatocytes_CV","Kupffer","Monocytes","Hepatocytes_PV","LECs",
                       "Stellate","T","Fibroblasts","B")

A1_1 <- subset(merged, idents = "Slide1_A1-1")
Idents(A1_1) <- "annotation"
for (i in cell_types_mets) {
  text_file_generation_anno_resolve(A1_1,i,"^Slide1_A1-1_Cell","A1_1")
}

A2_1 <- subset(merged, idents = "Slide1_A2-1")
Idents(A2_1) <- "annotation"
for (i in cell_types_noMets) {
  text_file_generation_anno_resolve(A2_1,i,"^Slide1_A2-1_Cell","A2_1")
}

A2_2 <- subset(merged, idents = "Slide1_A2-2")
Idents(A2_2) <- "annotation"
for (i in cell_types_noMets) {
  text_file_generation_anno_resolve(A2_2,i,"^Slide1_A2-2_Cell","A2_2")
}

B1_1 <- subset(merged, idents = "Slide1_B1-1")
Idents(B1_1) <- "annotation"
for (i in cell_types_mets) {
  text_file_generation_anno_resolve(B1_1,i,"^Slide1_B1-1_Cell","B1_1")
}

B1_2 <- subset(merged, idents = "Slide1_B1-2")
Idents(B1_2) <- "annotation"
for (i in cell_types_mets) {
  text_file_generation_anno_resolve(B1_2,i,"^Slide1_B1-2_Cell","B1_2")
}

B2_1 <- subset(merged, idents = "Slide1_B2-1")
Idents(B2_1) <- "annotation"
for (i in cell_types_noMets) {
  text_file_generation_anno_resolve(B2_1,i,"^Slide1_B2-1_Cell","B2_1")
}

##per slide all cell types 
Idents(merged) <- "Slide"
merged_slides <- c("^Slide1_A1-1_Cell","^Slide1_A2-1_Cell","^Slide1_A2-2_Cell","^Slide1_B1-1_Cell",
                   "^Slide1_B1-2_Cell","^Slide1_B2-1_Cell")

slides_names <- as.data.frame(table(merged$Slide))$Var1

text_file_generation_per_slide_resolve(merged,"Slide1_A1-1","^Slide1_A1-1_Cell")
text_file_generation_per_slide_resolve(merged,"Slide1_A2-1","^Slide1_A2-1_Cell")
text_file_generation_per_slide_resolve(merged,"Slide1_A2-2","^Slide1_A2-2_Cell")
text_file_generation_per_slide_resolve(merged,"Slide1_B1-1","^Slide1_B1-1_Cell")
text_file_generation_per_slide_resolve(merged,"Slide1_B1-2","^Slide1_B1-2_Cell")
text_file_generation_per_slide_resolve(merged,"Slide1_B2-1","^Slide1_B2-1_Cell")

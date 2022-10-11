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
merged <- FindClusters(merged, resolution = 2, random.seed = 2, algorithm = 1, graph.name = "originalexp_snn") 
merged <- RunUMAP(merged, dims = 1:10, seed.use = 5, reduction = "pca")
DimPlot(merged, reduction = "umap", label = TRUE) 

########## Annotation ##########
###check top DEGs from clusters and compare with sphere-seq data 
dge_genes <- FindAllMarkers(merged, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
View(dge_genes %>% group_by(cluster) %>% top_n(n =10, wt = avg_log2FC))

#load sphere-seq data 
liverSpS5C <- readRDS("/mnt/khandler/R_projects/Sphere-sequencing/Sphere-seq_analysis/data_files_generated/LiverMerged_afterBC_anno.Rda")
Idents(liverSpS5C) <- "annotation.broad"

#Clusters 0,1,20,21,27 --> Metastasis 
DotPlot(liverSpS5C, features =c("Pglyrp1","Gpx2","App","Spp1","Mecom","Jup"))

#Clusters 2,9,11,13,17 --> Hepatocytes 
DotPlot(liverSpS5C, features =c("Cyp2e1","Cyp1a2"))

#Clusters 2,24 -->  Kupffer 
DotPlot(liverSpS5C, features =c("Clec4f","Cd5l","Vsig4","C1qc"))

#Clusters 4 -->  Monocytes/Neutrophils --> mono marker top 
DotPlot(liverSpS5C, features =c("Lyz2","Csf3r","C1qc","Tgfbi","Cx3cr1"))

#Clusters 5,6,8,18 --> Hepatocytes 
DotPlot(liverSpS5C, features =c("Fgb","Cyp2f2","Cfh","Ghr"))

#Clusters 7-->  LECs
DotPlot(liverSpS5C, features =c("Galnt15","Acer2","Pecam1","Dll4","Cd36"))

#Clusters 10,12 -->   LEC
DotPlot(liverSpS5C, features =c("Plxnc1","Plpp1","Pecam1","Dab2"))

#Clusters 14 -->  Cholangiocytes 
DotPlot(liverSpS5C, features =c("App", "Pglyrp1","Spp1","Fn1","Dhrs3"))

#Clusters 15,23 --> Stellate  
DotPlot(liverSpS5C, features =c("Tgfbi","Plvap","Lrat","Cyp2e1","Vcam1","Cyp1a2"))

#Clusters 16 -->  NK/T/B --> T marker top 
DotPlot(liverSpS5C, features =c("Il2rb","Ighm","Ccl5","Sema4d","Itgb1"))

#Clusters 19 -->  Fibroblasts
DotPlot(liverSpS5C, features =c("Cald1","Fn1","Tgfbi","Spp1"))

#Clusters 22 -->   LEC
DotPlot(liverSpS5C, features =c("Cd34","Flt1","Dll4"))

#Clusters 25 -->   B/T/NK --> B marker top 
DotPlot(liverSpS5C, features =c("Ighm","Itga4","Il2rb"))

#Clusters 26 -->  could be anything 
DotPlot(liverSpS5C, features =c("App","Ltbr","Gprc5c","Psen1","Sema4a","Pglyrp1"))

#Clusters 29 -->  Hep/Cholangiocytes --> more Hep 
DotPlot(liverSpS5C, features =c("Spp1","Cald1", "Fgb","Acly","Gprc5c", "Dhrs3"))

#Clusters 30 -->   hep 
DotPlot(liverSpS5C, features =c("Cd5l","Cyp2f2","Fgb","Cfh"))

#Clusters 31 -->   Fibroblasts 
DotPlot(liverSpS5C, features =c("Fn1","Rrbp1","App"))

#Clusters 32 -->  Stellate/LECs --> more stellate 
DotPlot(liverSpS5C, features =c("Tgfbi","Ltbr","Nrp1","Dab2","Flt1","App","Rrbp1"))

#Clusters 33 -->  can be anything 
DotPlot(liverSpS5C, features =c("Sema4d","Ighm","Gprc5c","Fn1","Csf3r"))

#Clusters 34 -->   neutrophils top 
DotPlot(liverSpS5C, features =c("Csf3r","Jup","Dhrs3","Cald1","C1qb"))

#Clusters 35 -->   Metastasis  
DotPlot(liverSpS5C, features =c("Lhx6","Spp1","Pglyrp1","App","Fn1","Gpx2"))

#Clusters 36 -->   Stellate
DotPlot(liverSpS5C, features =c("Plvap","Tgfbi","Nrp1","App"))

#Clusters 37 -->  Cholangiocytes 
DotPlot(liverSpS5C, features =c("App","Spp1","Fn1","Pglyrp1"))

#Clusters 38,39,40 only express one gene each, difficult to put to a cell type --> remove 

current.cluster.ids <- c(0:40)
new.cluster.ids <- c("Metastasis","Metastasis","Hepatocytes_CV","Kupffer","Monocytes","Hepatocytes_PV","Hepatocytes_PV",
                     "LECs","Hepatocytes_PV","Hepatocytes_CV","LECs","Hepatocytes_CV","LECs",
                     "Hepatocytes_CV","Cholangiocytes","Stellate","T","Hepatocytes_CV",
                     "Hepatocytes_PV","Fibroblasts","Metastasis","Metastasis","LECs","Stellate",
                     "Kupffer","B","anything","Metastasis","Metastasis","Hepatocytes_PV","Hepatocytes_PV",
                     "Fibroblasts","Stellate","anything","Neutrophils","Metastasis","Stellate","Cholangiocytes","?","?","?")
merged$annotation <- plyr::mapvalues(x = merged$seurat_clusters, from = current.cluster.ids, to = new.cluster.ids)
DimPlot(merged, reduction = "umap", label = TRUE, group.by = "annotation", label.size = 3)

#remove ? and anything
Idents(merged) <- "annotation"
merged <- subset(merged, idents = c("Metastasis","Hepatocytes_CV","Kupffer","Monocytes","Hepatocytes_PV","LECs","Cholangiocytes",
                                    "Stellate","T","Fibroblasts","B","Neutrophils"))


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
               c("#F90606","#F46042","#E20FE8","#19E80F",
                 "#56040C","#F4E740","#EDB2D4","#C1A08A" ,"#2323E5",
                 "#8E5229","#EF975B","#AEAEAF"), pt.size = 0.5) + theme(legend.title = element_text(size = 22), legend.text = element_text(size = 22)) + 
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
#Cholangiocytes #EDB2D4
#Monocytes #19E80F 
#T #2323E5 
#B #EF975B
#Neutrophils #AEAEAF

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

cell_types_mets <- c("Metastasis","Hepatocytes_CV","Kupffer","Monocytes","Hepatocytes_PV","LECs","Cholangiocytes",
                     "Stellate","T","Fibroblasts","B","Neutrophils")
cell_types_noMets <- c("Hepatocytes_CV","Kupffer","Monocytes","Hepatocytes_PV","LECs","Cholangiocytes",
                       "Stellate","T","Fibroblasts","B","Neutrophils")

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

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
#merge samples 
merged <- merge(mets_samples, noMets_samples)

merged <- FindVariableFeatures(merged, selction.method = "vst", nfeatures = 100)
all.genes <- rownames(merged)
merged <- ScaleData(merged,features = all.genes)
merged <- RunPCA(object = merged, features = VariableFeatures(object = merged),approx=FALSE)
ElbowPlot(merged)

#plot PCA space to check for batch effect 
DimPlot(merged, reduction = "pca", group.by = "Slide") #no batch effect 

merged <- FindNeighbors(object = merged, dims = 1:10, reduction = "pca")
merged <- FindClusters(merged, resolution = 0.2, random.seed = 2, algorithm = 1, graph.name = "originalexp_snn")
merged <- RunUMAP(merged, dims = 1:10, seed.use = 5, reduction = "pca")
DimPlot(merged, reduction = "umap", label = TRUE)

#remove all the little random cluster 
Idents(merged) <- "seurat_clusters"
merged <- subset(merged, idents = c(0,1,2,3,4,5))
DimPlot(merged, reduction = "umap", label = TRUE) 

#cluster again 
merged <- FindNeighbors(object = merged, dims = 1:10, reduction = "pca")
merged <- FindClusters(merged, resolution = 2, random.seed = 2, algorithm = 1, graph.name = "originalexp_snn") 
merged <- RunUMAP(merged, dims = 1:10, seed.use = 5, reduction = "pca")
DimPlot(merged, reduction = "umap", label = TRUE) 

######### Annotation #########
#check DEGs from each cluster 
dge_genes <- FindAllMarkers(merged, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
View(dge_genes %>% group_by(cluster) %>% top_n(n =10, wt = avg_log2FC))

#compare top DGE with annotation of Sphere-seq data
liverSpS <- readRDS(file = "/mnt/khandler/R_projects/Sphere-sequencing/Sphere-seq_analysis/data_files_generated/LiverMerged_afterBC_anno_BS_5cells.Rda")

Idents(liverSpS) <- "annotation.broad"

#Clusters 0,1,20,21,27 --> Metastatic cells 
DotPlot(liverSpS, features =c("Pglyrp1","Gpx2","App","Spp1","Mecom","Jup"))

#Clusters 2,9,11,13,17 --> Hepatocytes 
DotPlot(liverSpS, features =c("Cyp2e1","Cyp1a2"))

#Clusters 2,24 -->  Kupffer cells
DotPlot(liverSpS, features =c("Clec4f","Cd5l","Vsig4","C1qc"))

#Clusters 4 -->  Monocytes/Neutrophils --> monocytes marker top 
DotPlot(liverSpS, features =c("Lyz2","Csf3r","C1qc","Tgfbi","Cx3cr1"))

#Clusters 5,6,8,18 --> Hepatocytes 
DotPlot(liverSpS, features =c("Fgb","Cyp2f2","Cfh","Ghr"))

#Clusters 7-->  LECs
DotPlot(liverSpS, features =c("Galnt15","Acer2","Pecam1","Dll4","Cd36"))

#Clusters 10,12 -->   LECs
DotPlot(liverSpS, features =c("Plxnc1","Plpp1","Pecam1","Dab2"))

#Clusters 14 -->  Cholangiocytes 
DotPlot(liverSpS, features =c("App", "Pglyrp1","Spp1","Fn1","Dhrs3"))

#Clusters 15,23 --> Stellate cells
DotPlot(liverSpS, features =c("Tgfbi","Plvap","Lrat","Cyp2e1","Vcam1","Cyp1a2"))

#Clusters 16 -->  NK/T/B --> T/NK marker top 
DotPlot(liverSpS, features =c("Il2rb","Ighm","Ccl5","Sema4d","Itgb1"))

#Clusters 19 -->  Fibroblasts
DotPlot(liverSpS, features =c("Cald1","Fn1","Tgfbi","Spp1"))

#Clusters 22 -->   LECs
DotPlot(liverSpS, features =c("Cd34","Flt1","Dll4"))

#Clusters 25 -->   B/T/NK --> B cell marker top 
DotPlot(liverSpS, features =c("Ighm","Itga4","Il2rb"))

#Clusters 16 -->  NK/T/B --> T/NK cell marker top 
DotPlot(liverSpS, features =c("Il2rb","Ighm","Ccl5","Sema4d","Itgb1"))

#Clusters 19 -->  Fibroblasts
DotPlot(liverSpS, features =c("Cald1","Fn1","Tgfbi","Spp1"))

#Clusters 22 -->   LECs
DotPlot(liverSpS, features =c("Cd34","Flt1","Dll4"))

#Clusters 25 -->   B/T/NK --> B cell marker top 
DotPlot(liverSpS, features =c("Ighm","Itga4","Il2rb"))

#Clusters 26 -->  could be anything 
DotPlot(liverSpS, features =c("App","Ltbr","Gprc5c","Psen1","Sema4a","Pglyrp1"))

#Clusters 29 -->  Hepatocytes/Cholangiocytes --> more Hepatocytes 
DotPlot(liverSpS, features =c("Spp1","Cald1", "Fgb","Acly","Gprc5c", "Dhrs3"))

#Clusters 30 -->   Hepatocytes 
DotPlot(liverSpS, features =c("Cd5l","Cyp2f2","Fgb","Cfh"))

#Clusters 31 -->   Fibroblasts 
DotPlot(liverSpS, features =c("Fn1","Rrbp1","App"))

#Clusters 32 -->  Stellate cells/LECs --> more Stellate cells 
DotPlot(liverSpS, features =c("Tgfbi","Ltbr","Nrp1","Dab2","Flt1","App","Rrbp1"))

#Clusters 33 -->  could be anything 
DotPlot(liverSpS, features =c("Sema4d","Ighm","Gprc5c","Fn1","Csf3r"))

#Clusters 34 -->   Neutrophils top 
DotPlot(liverSpS, features =c("Csf3r","Jup","Dhrs3","Cald1","C1qb"))

#Clusters 35 -->   Metastatic cells  
DotPlot(liverSpS, features =c("Lhx6","Spp1","Pglyrp1","App","Fn1","Gpx2"))

#Clusters 36 -->   Stellate cells
DotPlot(liverSpS, features =c("Plvap","Tgfbi","Nrp1","App"))

#Clusters 37 -->  Cholangiocytes 
DotPlot(liverSpS, features =c("App","Spp1","Fn1","Pglyrp1"))

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


########## plot Marker genes in DotPlot ##########
Idents(merged) <- "annotation"
marker_genes <- c("Pglyrp1","Gpx2","App","Spp1","Mecom","Jup", #Metastasis
                  "Cyp2e1","Cyp1a2" , #Hepatocytes CV 
                  "Fgb","Cyp2f2","Cfh",  #Hepatocytes PV
                  "Clec4f","Vsig4", #Kupffer cells
                  "Lyz2","Csf3r","C1qc","Tgfbi","Cx3cr1", #Monocytes, Neutrophils Csf3r top  
                  "Galnt15","Acer2","Pecam1","Dll4","Cd36", #LECs
                  "App", "Pglyrp1","Spp1","Fn1","Dhrs3", #Cholangiocytes
                  "Tgfbi","Plvap","Lrat","Cyp2e1","Vcam1","Cyp1a2", #Stellate cells
                  "Il2rb","Ighm","Ccl5","Sema4d","Itgb1", #T cells
                  "Cald1","Fn1","Tgfbi","Spp1", #Fibroblasts
                  "Ighm","Itga4","Il2rb" #B cells
)

marker_genes <- c("Pglyrp1","Gpx2", #Metastasis
                  "Cyp2e1","Cyp1a2" , #Hepatocytes CV 
                  "Fgb","Cyp2f2",  #Hepatocytes PV
                  "Clec4f","Vsig4", #Kupffer cells
                  "Lyz2","Csf3r","C1qc",#Monocytes, Neutrophils Csf3r top  
                  "Galnt15","Acer2","Pecam1",#LECs
                  "App","Spp1",#Cholangiocytes
                  "Tgfbi","Plvap","Lrat", #Stellate cells
                  "Il2rb","Ccl5","Cd3d", #T cells
                  "Cald1","Fn1", #Fibroblasts
                  "Ighm","Itga4" #B cells
)

Idents(merged) <- "annotation"

p <- DotPlot(merged,features = marker_genes) + theme(axis.text.x = element_text(angle = 90)) + 
  ggtitle("Marker genes in annotation of Resolve data") + 
  theme(legend.title = element_text(size = 22), legend.text = element_text(size = 22)) + 
  theme(axis.title= element_text(size = 25)) + theme(axis.text = element_text(size = 15))+ 
  theme(plot.title = element_text(size = 25, face = "bold")) 
p + ggsave("./figures/3/marker_dotPlot.pdf",width = 12, height = 10)
p + ggsave("./figures/3/marker_dotPlot.svg",width = 12, height = 10)

########## plot annotated UMAP plots ##########
#colors: 
#Hepatocytes: CV #F46042 PV #56040C
#Kupffer cells: #E20FE8   
#LECs:  #F4E740
#Metastatic cells: #F90606
#Stellate cells: #C1A08A
#Fibroblasts: #8E5229 
#Cholangiocytes: #EDB2D4
#Monocytes: #19E80F 
#T cells: #2323E5 
#B cells: #EF975B
#Neutrophils: #AEAEAF

p <- DimPlot(merged, reduction = "umap", label = TRUE, group.by = "annotation", label.size = 3, cols = 
               c("#F90606","#F46042","#E20FE8","#19E80F",
                 "#56040C","#F4E740","#EDB2D4","#C1A08A" ,"#2323E5",
                 "#8E5229","#EF975B","#AEAEAF"), pt.size = 0.5) + theme(legend.title = element_text(size = 22), legend.text = element_text(size = 22)) + 
  theme(title = element_text(size = 25))+ theme(axis.text = element_text(size = 30)) 
p + ggsave("./figures/3/annotated_umap.pdf",width = 15, height = 10)
p + ggsave("./figures/3/annotated_umap.svg",width = 15, height = 10)

######### add vein and Mets_distance identity ##########
#combine the surrounding PV/CV and CV/PV to far and Mets to close 
current.cluster.ids <- c("cv","cv_nm","cv_sm","mets","pv","pv_nm","pv_sm")
new.cluster.ids <- c("CV","CV","CV","Mets","PV","PV","PV")
merged$vein <- plyr::mapvalues(x = merged$spatial_feature, from = current.cluster.ids, to = new.cluster.ids)

current.cluster.ids <- c("cv","cv_nm","cv_sm","mets","pv","pv_nm","pv_sm")
new.cluster.ids <- c("far","far","far","close","far","far","far")
merged$Mets_distance <- plyr::mapvalues(x = merged$spatial_feature, from = current.cluster.ids, to = new.cluster.ids)

#change sample IDs depending on visible metastasis, mets..=visible metastasis; no_mets...=no visible metastasis 
current.cluster.ids <- c("Slide1_A1-1" ,"Slide1_A2-1"   ,"Slide1_A2-2"     ,"Slide1_B1-1" ,"Slide1_B1-2"  ,"Slide1_B2-1" )
new.cluster.ids <- c("metsA1","no_metsA2","no_metsA2","metsB1","metsB1","no_metsB2")
merged$sampleID <- plyr::mapvalues(x = merged$Slide, from = current.cluster.ids, to = new.cluster.ids)

###plot split UMAP plot 
p <- DimPlot(merged, reduction = "umap", label = TRUE, group.by = "annotation",split.by = "Mets_distance", label.size = 3, cols = 
               c("#F90606","#F46042","#E20FE8","#19E80F",
                 "#56040C","#F4E740","#EDB2D4","#C1A08A" ,"#2323E5",
                 "#8E5229","#EF975B","#AEAEAF"),pt.size = 0.5) + theme(legend.title = element_text(size = 22), legend.text = element_text(size = 22)) + 
  theme(title = element_text(size = 25))+ theme(axis.text = element_text(size = 30)) 
p + ggsave("./figures/3/annotated_umap_split_mets.pdf",width = 15, height = 10)
p + ggsave("./figures/3/annotated_umap_split_mets.svg",width = 15, height = 10)

######### check CV/PV manual drawing based on Hepatocyte landmark genes ##########
Idents(merged) <- "annotation"
hep <- subset(merged, idents = c("Hepatocytes_CV","Hepatocytes_PV"))
p <- DimPlot(hep,group.by = "annotation", split.by = "vein", cols = c("#56040C","#F46042"), pt.size = 0.5) + theme(legend.title = element_text(size = 22), legend.text = element_text(size = 22)) + 
  theme(title = element_text(size = 25))+ theme(axis.text = element_text(size = 30)) 
p + ggsave("./figures/3/split_plot_pv_cv_hep.pdf",width = 15, height = 10)
p + ggsave("./figures/3/split_plot_pv_cv_hep.svg",width = 15, height = 10)

########## save R object ########## 
saveRDS(merged,file = "./data_files_generated/Resolve_seurat_anno.rds")

########## Get annotation of cell types per slide and for each cell ##########
#this is used for the overlay with the Dapi image in ImageJ 
setwd("./data_files_generated/annotation_file_for_ImageJ")
Idents(merged) <- "Slide"

cell_types_mets <- c("Metastasis","Hepatocytes_CV","Kupffer","Monocytes","Hepatocytes_PV","LECs","Cholangiocytes",
                     "Stellate","NK_T","Fibroblasts","B","Neutrophils")
cell_types_noMets <- c("Hepatocytes_CV","Kupffer","Monocytes","Hepatocytes_PV","LECs","Cholangiocytes",
                       "Stellate","NK_T","Fibroblasts","B","Neutrophils")

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




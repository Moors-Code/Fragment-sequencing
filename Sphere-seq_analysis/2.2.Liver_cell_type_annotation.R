########## Part 2.2: Cell type annotation of liver samples ##########
#This part annotates clustered Seurat object using cell type markers from https://www.livercellatlas.org and based on top DEGs between clusters 
#Annotation is done broadly first and then broad cell type clusters are subclustered to annotate cell type subtypes 

########## Prepare environment ##########
###Setting the working directory 
setwd("/mnt/khandler/R_projects/Sphere-sequencing/Sphere-seq_analysis/")

###Load packages and functions 
source("./functions_and_packages/1.Packages.R")

###load R object 
merged_cl <- readRDS(file = "./data_files_generated/LiverMerged_afterBC.Rda")

########## Broad annotation ##########
Idents(merged_cl) <- "seurat_clusters"
p <- DotPlot(merged_cl, features = c(
  "Ighm", "Igkc", #B cells: clusters 6,14,25
  "Trac","Cd3d" ,#T cells: clusters 8,11
  "Klra8","Cma1" ,#NK cells: clusters 8
  "Clec4f","Vsig4", #Kupffer cells: clusters 1,2,3,9,22
  "Bcl11a", "Ccr9", #DCs: cluster 21
  "S100a4", "Itgax", #general Monocytes: clusters 10,12,15
  "Hba-a2", #Red blood cells: cluter 19
  "Csf3r", "S100a8",#Neutrophils: cluster 17
  "Siglech","Cox6a2", #pDCs: cluster 21
  "Fcer1a","Cpa3",#Basophils: cluster 17
  "Ttr", "Alb", "Cyp2e1", #Hepatocytes: cluster 18
  "Pecam1","Dll4", #LSECs: clusters 0,4,5,7,13,23,24
  "Lrat","Reln", #Stellate cells: cluster 20
  "Carmn","Nr1h5",  #Stromal cells: cluster 20
  "Svep1","Ncam1", #Fibroblasts: cluster 20
  "Krt19", "Epcam","Sox9", #Cholangiocytes: cluster  20
  "Sprr2a2","Pglyrp1","Gpx2" #Metastatic cells: cluster 16
)) 
p + theme(legend.title = element_text(size = 6), legend.text = element_text(size = 6)) + 
  theme(title = element_text(size = 10))+ theme(axis.text = element_text(size = 8)) +
  ggtitle("Marker genes for broad annotation") + theme(axis.text.x = element_text(angle = 90)) 

#Add idents of broad annotation 
Idents(merged_cl) <- "seurat_clusters"
current.cluster.ids <- c(0:25)
new.cluster.ids <- c("Endothelial","Kupffer","Kupffer","Kupffer","Endothelial","Endothelial","B","Endothelial","T_NK",
                     "Kupffer","Monocytes_DC","T_NK","Monocytes_DC","Endothelial","B","Monocytes_DC","Mets_Chol","Baso_Neut",
                     "Hepatocytes","RBC","Fib_Stellate","Monocytes_DC","Kupffer","Endothelial","Endothelial","B")
merged_cl$annotation <- plyr::mapvalues(x = merged_cl$seurat_clusters, from = current.cluster.ids, to = new.cluster.ids)

########## Fine grade annotation ##########
###Monocytes/DCs 
Idents(merged_cl) <- "annotation"
merged_cl_subCl <- FindSubCluster(merged_cl, "Monocytes_DC",graph.name = "RNA_snn", subcluster.name = "sub.cluster",resolution = 0.4, algorithm = 4)
DimPlot(merged_cl_subCl, reduction = "umap", label = TRUE, group.by = "sub.cluster")

Idents(merged_cl_subCl) <- "sub.cluster"
sub_celltype <- subset(merged_cl_subCl,idents = c( "Monocytes_DC_1","Monocytes_DC_2","Monocytes_DC_3","Monocytes_DC_4","Monocytes_DC_5",
                                                   "Monocytes_DC_6","Monocytes_DC_7"))
DimPlot(sub_celltype, reduction = "umap", label = TRUE, group.by = "sub.cluster")

p <- DotPlot(sub_celltype, features = c("S100a4", "Itgax", "Crip1","Naaa","Adgre1","Axl","Mafb","Cx3cr1","C5ar1","Chil3","Sell","Gm9733", # general monocyte-derived cells
                                        "Pglyrp1","Spn","Trem3", #Patrolling monocytes  
                                        "Clec9a","Cd24a","Xcr1", #cDC1: cluster 5
                                        "Cd209a","Mgl2","Clec10a","Cd7","Tnfsf9" , #cDC2: cluster 4
                                        "Siglech","Cox6a2",#pDCs: cluster 6
                                        "Clec4f","Vsig4" #Kupffer cells: cluster 7
)) 
#Monocytes: cluster 1,2,3
p + theme(legend.title = element_text(size = 6), legend.text = element_text(size = 6)) + 
  theme(title = element_text(size = 10))+ theme(axis.text = element_text(size = 6)) +
  ggtitle("Marker genes Monocytes/DC") + theme(axis.text.x = element_text(angle = 90)) 

Idents(merged_cl_subCl) <- "sub.cluster"
current.cluster.ids <- c("B","Baso_Neut","Endothelial","Fib_Stellate","Hepatocytes","Kupffer","Mets_Chol","Monocytes_DC_1",
                         "Monocytes_DC_2","Monocytes_DC_3","Monocytes_DC_4","Monocytes_DC_5","Monocytes_DC_6","Monocytes_DC_7",
                         "RBC","T_NK")
new.cluster.ids <- c("B","Baso_Neut","Endothelial","Fib_Stellate","Hepatocytes","Kupffer","Mets_Chol","Monocytes",
                     "Monocytes","Monocytes","cDC2","cDC1","pDC","Kupffer",
                     "RBC","T_NK")
merged_cl_subCl$annotation <- plyr::mapvalues(x = merged_cl_subCl$sub.cluster, from = current.cluster.ids, to = new.cluster.ids)
DimPlot(merged_cl_subCl, reduction = "umap", label = TRUE, group.by = "annotation", label.size = 3)

###Monocytes
Idents(merged_cl_subCl) <- "annotation"
merged_cl_subCl <- FindSubCluster(merged_cl_subCl, "Monocytes",graph.name = "RNA_snn", subcluster.name = "sub.cluster",resolution = 0.2, algorithm = 4)
DimPlot(merged_cl_subCl, reduction = "umap", label = TRUE, group.by = "sub.cluster")

Idents(merged_cl_subCl) <- "sub.cluster"
sub_celltype <- subset(merged_cl_subCl,idents = c( "Monocytes_1","Monocytes_2","Monocytes_3"))
DimPlot(sub_celltype, reduction = "umap", label = TRUE, group.by = "sub.cluster") 

p <- DotPlot(sub_celltype, features = c("Pglyrp1","Spn","Trem3","Ly6c1", #Patrolling monocytes: cluster 2
                                        "Tgfb2","Cbr2", #Peritoneal Macrophages
                                        "Gpnmb","Trem2","Lgals1","Ftl1", #Bile-duct LAMs (lipid associated macrophages)
                                        "Clec4f","Vsig4", #Kupffer cells 
                                        "Thbs1","Clec4l","Cd209a", #Thbs1+ macrophages Qi 2022 et al, Nature Communication (Thbs1 activates M1-like TAMS ) 
                                        "Vcan","Anpep", #Vcan+ macrophages Qi 2022 et al, Nature Communication
                                        "Ly6c2","Lyz2", #Ly6c+ macrophages: cluster 1 
                                        "C1qc","C1qb","C1qa" #C1q+ macrophages: cluster 3 
)) 
p + theme(legend.title = element_text(size = 6), legend.text = element_text(size = 6)) + 
  theme(title = element_text(size = 10))+ theme(axis.text = element_text(size = 6)) +
  ggtitle("Marker genes Monocytes") + theme(axis.text.x = element_text(angle = 90)) + 
  ggsave("./figures/2.2/DotPlot_monocytes_anno.pdf")

dge_genes <- FindAllMarkers(sub_celltype, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
View(dge_genes %>% group_by(cluster) %>% top_n(n =10, wt = avg_log2FC))


Idents(merged_cl_subCl) <- "sub.cluster"
current.cluster.ids <- c("B","Baso_Neut","cDC1","cDC2", "Endothelial","Fib_Stellate","Hepatocytes","Kupffer","Mets_Chol","Monocytes_1",
                         "Monocytes_2","Monocytes_3","pDC",
                         "RBC","T_NK")
new.cluster.ids <- c("B","Baso_Neut","cDC1","cDC2", "Endothelial","Fib_Stellate","Hepatocytes","Kupffer","Mets_Chol","Mac_Ly6c",
                     "Mono_patrolling","Mac_C1q","pDC",
                     "RBC","T_NK")
merged_cl_subCl$annotation <- plyr::mapvalues(x = merged_cl_subCl$sub.cluster, from = current.cluster.ids, to = new.cluster.ids)
DimPlot(merged_cl_subCl, reduction = "umap", label = TRUE, group.by = "annotation", label.size = 3)


###Metastatic cells/Cholangiocytes
merged_cl <- merged_cl_subCl
Idents(merged_cl) <- "annotation"
merged_cl_subCl <- FindSubCluster(merged_cl, "Mets_Chol",graph.name = "RNA_snn", subcluster.name = "sub.cluster",resolution = 0.1, algorithm = 4)
DimPlot(merged_cl_subCl, reduction = "umap", label = TRUE, group.by = "sub.cluster")

Idents(merged_cl_subCl) <- "sub.cluster"
sub_celltype <- subset(merged_cl_subCl,idents = c( "Mets_Chol_1","Mets_Chol_2","Mets_Chol_3"))
DimPlot(sub_celltype, reduction = "umap", label = TRUE, group.by = "sub.cluster") 

p <- DotPlot(sub_celltype, features = c(  "Krt19", "Epcam","Sox9","Spp1", #Cholangiocytes: cluster  2
                                          "Sprr2a2","Pglyrp1","Gpx2" #Metastatic cells: clusters 1,3
))
p + theme(legend.title = element_text(size = 6), legend.text = element_text(size = 6)) + 
  theme(title = element_text(size = 10))+ theme(axis.text = element_text(size = 8)) + theme(axis.text.x = element_text(angle = 90)) 

#rename clusters that could be annotated 
Idents(merged_cl) <- "annotation"
current.cluster.ids <- c("B","Baso_Neut","cDC1","cDC2", "Endothelial","Fib_Stellate","Hepatocytes","Kupffer","Mets_Chol_1",
                         "Mets_Chol_2","Mets_Chol_3", "Mac_C1q","Mono_patrolling","Mac_Ly6c",
                         "pDC","RBC","T_NK")
new.cluster.ids <- c("B","Baso_Neut","cDC1","cDC2", "Endothelial","Fib_Stellate","Hepatocytes","Kupffer","Metastasis",
                     "Cholangiocytes","Metastasis", "Mac_C1q","Mono_patrolling","Mac_Ly6c",
                     "pDC","RBC","T_NK")
merged_cl_subCl$annotation <- plyr::mapvalues(x = merged_cl_subCl$sub.cluster, from = current.cluster.ids, to = new.cluster.ids)
DimPlot(merged_cl_subCl, reduction = "umap", label = TRUE, group.by = "annotation", label.size = 3) 

###Endothelial cells 
merged_cl <- merged_cl_subCl
Idents(merged_cl) <- "annotation"
merged_cl_subCl <- FindSubCluster(merged_cl, "Endothelial",graph.name = "RNA_snn", subcluster.name = "sub.cluster",resolution = 0.2, algorithm = 4)
DimPlot(merged_cl_subCl, reduction = "umap", label = TRUE, group.by = "sub.cluster")

Idents(merged_cl_subCl) <- "sub.cluster"
sub_celltype <- subset(merged_cl_subCl, idents = c( "Endothelial_1","Endothelial_2","Endothelial_3","Endothelial_4",
                                                    "Endothelial_5","Endothelial_6","Endothelial_7"))
DimPlot(sub_celltype, reduction = "umap", label = TRUE, group.by = "sub.cluster", label.size = 3) 

p <- DotPlot(sub_celltype, features = c("Lyve1","Flt4","Stab2", #LSECs: clusters 2,1,4,6
                                        "Il33","Pdgfb","Cd9","Timp3", #LVECs: cluster 3
                                        "Rspo3","Wnt9b", #Central vein ECs: cluster 3
                                        "Gja5", #Protain vein ECs: cluster 3
                                        "Ccl21a","Mmrn1","Thy1" #LECs (lymphatic ECs): cluster 3
)) 
#remove: clusters 5,7
p + theme(legend.title = element_text(size = 6), legend.text = element_text(size = 6)) + 
  theme(title = element_text(size = 10))+ theme(axis.text = element_text(size = 8)) + theme(axis.text.x = element_text(angle = 90)) 

#rename clusters that could be annotated 
Idents(merged_cl_subCl) <- "sub.cluster"
current.cluster.ids <- c("B","Baso_Neut","cDC1","cDC2", "Endothelial_1", "Endothelial_2", "Endothelial_3", "Endothelial_4",
                         "Endothelial_5", "Endothelial_6", "Endothelial_7","Fib_Stellate","Hepatocytes","Kupffer","Metastasis",
                         "Mac_C1q", "Mono_patrolling","Cholangiocytes", "Mac_Ly6c","pDC","RBC","T_NK")
new.cluster.ids <- c("B","Baso_Neut","cDC1","cDC2", "LSECs", "LSECs", "ECs", "LSECs",
                     "remove", "LSECs", "remove","Fib_Stellate","Hepatocytes","Kupffer","Metastasis",
                     "Mac_C1q", "Mono_patrolling","Cholangiocytes", "Mac_Ly6c","pDC","RBC","T_NK")
merged_cl_subCl$annotation <- plyr::mapvalues(x = merged_cl_subCl$sub.cluster, from = current.cluster.ids, to = new.cluster.ids)
DimPlot(merged_cl_subCl, reduction = "umap", label = TRUE, group.by = "annotation", label.size = 3)

###Endothelial cells that are not LSECs (ECs)
merged_cl <- merged_cl_subCl
Idents(merged_cl) <- "annotation"
merged_cl_subCl <- FindSubCluster(merged_cl, "ECs",graph.name = "RNA_snn", subcluster.name = "sub.cluster",resolution = 0.4, algorithm = 4)
DimPlot(merged_cl_subCl, reduction = "umap", label = TRUE, group.by = "sub.cluster")

Idents(merged_cl_subCl) <- "sub.cluster"
sub_celltype <- subset(merged_cl_subCl, idents = c( "ECs_1","ECs_2","ECs_3","ECs_4","ECs_5","ECs_6"))
DimPlot(sub_celltype, reduction = "umap", label = TRUE, group.by = "sub.cluster", label.size = 3) 

p <- DotPlot(sub_celltype, features = c("Lyve1","Flt4","Stab2", #LSECs: clusters 1,3,4
                                        "Rspo3","Wnt9b", #Central vein ECs clusters 
                                        "Gja5", #Protain vein ECs 
                                        "Ccl21a","Mmrn1","Thy1", #LECs (lymphatic ECs) 
                                        "Il33","Pdgfb","Cd9","Timp3" #LVECs: clusters 5,6,2
)) 
p + theme(legend.title = element_text(size = 6), legend.text = element_text(size = 6)) + 
  theme(title = element_text(size = 10))+ theme(axis.text = element_text(size = 8)) +theme(axis.text.x = element_text(angle = 90)) 

#rename clusters that could be annotated 
Idents(merged_cl_subCl) <- "sub.cluster"
current.cluster.ids <- c("B","Baso_Neut","cDC1","cDC2", "ECs_1", "ECs_2","ECs_3","ECs_4","ECs_5","ECs_6","Fib_Stellate",
                         "Hepatocytes","Kupffer","Metastasis","LSECs",
                         "Mac_C1q", "Mono_patrolling","Cholangiocytes", "Mac_Ly6c","pDC","RBC","T_NK")
new.cluster.ids <- c("B","Baso_Neut","cDC1","cDC2", "LSECs", "LVECs","LSECs","LVECs","LVECs","LVECs","Fib_Stellate",
                     "Hepatocytes","Kupffer","Metastasis","LSECs",
                     "Mac_C1q", "Mono_patrolling","Cholangiocytes", "Mac_Ly6c","pDC","RBC","T_NK")
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
p <- DotPlot(sub_celltype, features = c( "Trac", "Trbc1", "Trbc2","Cd3d","Themis" ,#T cells: clusters 1,3,4,5,6
                                         "Cma1","Ncr1","Itga2","Itgam"#NK cells: cluster 2
))
p + theme(legend.title = element_text(size = 6), legend.text = element_text(size = 6)) + 
  theme(title = element_text(size = 10))+ theme(axis.text = element_text(size = 8)) + theme(axis.text.x = element_text(angle = 90)) 

#rename clusters that could be annotated 
Idents(merged_cl_subCl) <- "sub.cluster"
current.cluster.ids <- c("B","Baso_Neut","cDC1","cDC2","Fib_Stellate",
                         "Hepatocytes","Kupffer","LSECs","LVECs", "Metastasis",
                         "Mac_C1q", "Mono_patrolling","Cholangiocytes", "Mac_Ly6c","pDC","RBC","remove","T_NK_1","T_NK_2","T_NK_3",
                         "T_NK_4","T_NK_5","T_NK_6")
new.cluster.ids <- c("B","Baso_Neut","cDC1","cDC2", "Fib_Stellate",
                     "Hepatocytes","Kupffer","LSECs","LVECs", "Metastasis",
                     "Mac_C1q", "Mono_patrolling","Cholangiocytes", "Mac_Ly6c","pDC","RBC","remove","T","NK","T",
                     "T","T","T")
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
current.cluster.ids <- c("B_1","B_2","B_3","Baso_Neut","cDC1","cDC2","Fib_Stellate",
                         "Hepatocytes","Kupffer","LSECs","LVECs","Metastasis",
                         "Mac_C1q", "Mono_patrolling","Cholangiocytes", "Mac_Ly6c","NK",
                         "pDC","RBC","remove","T")
new.cluster.ids <- c("B_mem","B_mem","B_plasma","Baso_Neut","cDC1","cDC2","Fib_Stellate",
                     "Hepatocytes","Kupffer","LSECs","LVECs", "Metastasis",
                     "Mac_C1q", "Mono_patrolling","Cholangiocytes", "Mac_Ly6c","NK",
                     "pDC","RBC","remove","T")
merged_cl_subCl$annotation <- plyr::mapvalues(x = merged_cl_subCl$sub.cluster, from = current.cluster.ids, to = new.cluster.ids)
DimPlot(merged_cl_subCl, reduction = "umap", label = TRUE, group.by = "annotation", label.size = 3)

###Fibroblasts and Stellate cells
merged_cl <- merged_cl_subCl
Idents(merged_cl) <- "annotation"
merged_cl_subCl <- FindSubCluster(merged_cl, c("Fib_Stellate"),graph.name = "RNA_snn", subcluster.name = "sub.cluster",resolution = 0.4, algorithm = 4)
DimPlot(merged_cl_subCl, reduction = "umap", label = TRUE, group.by = "sub.cluster")

Idents(merged_cl_subCl) <- "sub.cluster"
sub_celltype <- subset(merged_cl_subCl,idents = c("Fib_Stellate_1","Fib_Stellate_2","Fib_Stellate_3","Fib_Stellate_4"))
DimPlot(sub_celltype, reduction = "umap", label = TRUE, group.by = "sub.cluster", label.size = 3) 

p <- DotPlot(sub_celltype, features = c("Lrat", "Rbp1","Nt5e","Itgb1","Dach1","Nrxn1","Dcdc2a","Reln", #Stellate cells: cluster 4
                                        "Svep1","Ncam1" #Fibroblasts: clusters 1,2,3
                                        
))
p + theme(legend.title = element_text(size = 6), legend.text = element_text(size = 6)) + 
  theme(title = element_text(size = 10))+ theme(axis.text = element_text(size = 6)) + theme(axis.text.x = element_text(angle = 90)) 

#rename clusters that could be annotated 
Idents(merged_cl_subCl) <- "sub.cluster"
current.cluster.ids <- c("B_mem","B_plasma","Baso_Neut","cDC1","cDC2","Fib_Stellate_1","Fib_Stellate_2","Fib_Stellate_3","Fib_Stellate_4",
                         "Hepatocytes","Kupffer","LSECs","LVECs", "Metastasis",
                         "Mac_C1q", "Mono_patrolling","Cholangiocytes", "Mac_Ly6c","NK",
                         "pDC","RBC","remove","T")
new.cluster.ids <- c("B_mem","B_plasma","Baso_Neut","cDC1","cDC2","Fibroblasts","Fibroblasts","Fibroblasts","Stellate",
                     "Hepatocytes","Kupffer","LSECs","LVECs", "Metastasis",
                     "Mac_C1q", "Mono_patrolling","Cholangiocytes", "Mac_Ly6c","NK",
                     "pDC","RBC","remove","T")
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
p <- DotPlot(sub_celltype, features = c("Clec4f","Vsig4","Timd4","Itgal","Cd5l", "Slc16a9","Slc40a1", #Kupffer cells (KC1): cluster 1,3,4
                                        "Esam","Mrc1", "Cd36","Cd63", "Cd81", "Lamp1", #Kupffer cells 2 (KC2) 
                                        "Ldb2","Rasip1","Ddx19b","Champ1","Tmem88","Cxcr4","Igfbp7","Clex4g" #Endothelial
)) 
#Kupffer_Endo: cluster 2 (no so clear if these are KC2 or KC/LECs doublets therefore annotate as Kupffer_Endo doublets), remove cluster 5 
p + theme(legend.title = element_text(size = 6), legend.text = element_text(size = 6)) + 
  theme(title = element_text(size = 10))+ theme(axis.text = element_text(size = 6)) + theme(axis.text.x = element_text(angle = 90)) 

#rename clusters that could be annotated 
Idents(merged_cl_subCl) <- "sub.cluster"
current.cluster.ids <- c("B_mem","B_plasma","Baso_Neut","cDC1","cDC2","Fibroblasts",
                         "Hepatocytes","Kupffer_1","Kupffer_2","Kupffer_3","Kupffer_4","Kupffer_5","LSECs","LVECs", "Metastasis",
                         "Mac_C1q", "Mono_patrolling","Cholangiocytes", "Mac_Ly6c","NK",
                         "pDC","RBC","remove","Stellate", "T")
new.cluster.ids <- c("B_mem","B_plasma","Baso_Neut","cDC1","cDC2","Fibroblasts",
                     "Hepatocytes","Kupffer","Kupffer_Endo","Kupffer","Kupffer","remove","LSECs","LVECs", "Metastasis",
                     "Mac_C1q", "Mono_patrolling","Cholangiocytes", "Mac_Ly6c","NK",
                     "pDC","RBC","remove","Stellate", "T")
merged_cl_subCl$annotation <- plyr::mapvalues(x = merged_cl_subCl$sub.cluster, from = current.cluster.ids, to = new.cluster.ids)
DimPlot(merged_cl_subCl, reduction = "umap", label = TRUE, group.by = "annotation", label.size = 3)

###Basophils/Neutrophils
merged_cl <- merged_cl_subCl
Idents(merged_cl) <- "annotation"
merged_cl_subCl <- FindSubCluster(merged_cl, c("Baso_Neut"),graph.name = "RNA_snn", subcluster.name = "sub.cluster",resolution = 0.05, algorithm = 4)
DimPlot(merged_cl_subCl, reduction = "umap", label = TRUE, group.by = "sub.cluster")

Idents(merged_cl_subCl) <- "sub.cluster"
sub_celltype <- subset(merged_cl_subCl,idents = c("Baso_Neut_1","Baso_Neut_2"))
DimPlot(sub_celltype, reduction = "umap", label = TRUE, group.by = "sub.cluster", label.size = 3) 

p <- DotPlot(sub_celltype, features = c("Csf3r", "S100a8",#Neutrophils: cluster 1
                                        "Fcer1a","Cpa3","Il6" #Basophils: cluster 2
)) 

p + theme(legend.title = element_text(size = 6), legend.text = element_text(size = 6)) + 
  theme(title = element_text(size = 10))+ theme(axis.text = element_text(size = 6)) + theme(axis.text.x = element_text(angle = 90)) 

Idents(merged_cl_subCl) <- "sub.cluster"
current.cluster.ids <- c("B_mem","B_plasma","Baso_Neut_1","Baso_Neut_2","cDC1","cDC2","Fibroblasts",
                         "Hepatocytes","Kupffer","Kupffer_Endo","LSECs","LVECs", "Metastasis",
                         "Mac_C1q", "Mono_patrolling","Cholangiocytes", "Mac_Ly6c","NK",
                         "pDC","RBC","remove","Stellate", "T")
new.cluster.ids <- c("B_mem","B_plasma","Neutrophils","Basophils","cDC1","cDC2","Fibroblasts",
                     "Hepatocytes","Kupffer","Kupffer_Endo","LSECs","LVECs", "Metastasis",
                     "Mac_C1q", "Mono_patrolling","Cholangiocytes", "Mac_Ly6c","NK",
                     "pDC","RBC","remove","Stellate", "T")
merged_cl_subCl$annotation <- plyr::mapvalues(x = merged_cl_subCl$sub.cluster, from = current.cluster.ids, to = new.cluster.ids)
DimPlot(merged_cl_subCl, reduction = "umap", label = TRUE, group.by = "annotation", label.size = 3)

##T cells 
merged_cl <- merged_cl_subCl
Idents(merged_cl) <- "annotation"
merged_cl_subCl <- FindSubCluster(merged_cl, c("T"),graph.name = "RNA_snn", subcluster.name = "sub.cluster",resolution = 0.4, algorithm = 4)
DimPlot(merged_cl_subCl, reduction = "umap", label = TRUE, group.by = "sub.cluster")

Idents(merged_cl_subCl) <- "sub.cluster"
sub_celltype <- subset(merged_cl_subCl,idents = c("T_1","T_2","T_3","T_4","T_5"))
DimPlot(sub_celltype, reduction = "umap", label = TRUE, group.by = "sub.cluster", label.size = 3) 

#Cd59b = Itga2, Cd11b = Itgam, Cd11c = Itgax, Tcrb = Trb, Cd11a = Itgal, Cd278 = Icos, Cd73 = Nt5e, Cd62l = Sell , Cd43 = Spn, Cd49d = Klra32
p <- DotPlot(sub_celltype, features = c( "Trac", "Trbc1", "Trbc2","Cd3d","Themis" ,#T cells general 
                                         "Klra8","Cma1","Ncr1","Itga2","Itgam","Klra32",#NK cells
                                         "Itgax", #NKT cells: 4
                                         "Gzmc","Klrb1b", #ILC1s 
                                         "Cd5","Trb","Cd4", #CD4+ T cells: clusters 2,1,5
                                         "Itgal","Cd90", #Th1s
                                         "Icos","Nt5e","Itgb8", #Th17
                                         "Cd8a","Sell","Cd8b", #CD+8 T cells: cluster 3
                                         "Ly6c1","Spn" #TEMs
))

p + theme(legend.title = element_text(size = 6), legend.text = element_text(size = 6)) + 
  theme(title = element_text(size = 10))+ theme(axis.text = element_text(size = 6)) + theme(axis.text.x = element_text(angle = 90)) 

#rename clusters that could be annotated 
Idents(merged_cl_subCl) <- "sub.cluster"
current.cluster.ids <- c("B_mem","B_plasma","Neutrophils","Basophils","cDC1","cDC2","Fibroblasts",
                         "Hepatocytes","Kupffer","Kupffer_Endo","LSECs","LVECs", "Metastasis",
                         "Mac_C1q", "Mono_patrolling","Cholangiocytes", "Mac_Ly6c","NK",
                         "pDC","RBC","remove","Stellate", "T_1","T_2","T_3","T_4","T_5")
new.cluster.ids <- c("B_mem","B_plasma","Neutrophils","Basophils","cDC1","cDC2","Fibroblasts",
                     "Hepatocytes","Kupffer","Kupffer_Endo","LSECs","LVECs", "Metastasis",
                     "Mac_C1q", "Mono_patrolling","Cholangiocytes", "Mac_Ly6c","NK",
                     "pDC","RBC","remove","Stellate", "T_CD4","T_CD4","T_CD8","NKT","T_CD4")
merged_cl_subCl$annotation <- plyr::mapvalues(x = merged_cl_subCl$sub.cluster, from = current.cluster.ids, to = new.cluster.ids)
DimPlot(merged_cl_subCl, reduction = "umap", label = TRUE, group.by = "annotation", label.size = 3)


########## Remove RBCs and low quality cells (annotated as remove) ##########
merged_cl <- merged_cl_subCl
Idents(merged_cl) <- "annotation"
merged_cl_sub <- subset(merged_cl, idents = c("B_mem","B_plasma","Neutrophils","Basophils","cDC1","cDC2","Fibroblasts",
                                              "Hepatocytes","Kupffer","Kupffer_Endo","LSECs","LVECs", "Metastasis",
                                              "Mac_C1q", "Mono_patrolling","Cholangiocytes", "Mac_Ly6c","NK","NKT",
                                              "pDC","Stellate", "T_CD4","T_CD8"))
DimPlot(merged_cl_sub, reduction = "umap", label = TRUE, group.by = "annotation", label.size = 3)

########## Plot Annotated umap ##########
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

p <- DimPlot(merged_cl_sub, label = TRUE,label.size = 5, group.by = "annotation", pt.size = 0.5, 
             cols = c("#EF975B","#F4741E","#56595B","#91C6C4",
                      "#315B5A","#EDB2D4",
                      "#8E5229","#F46042","#E20FE8",
                      "#C491C6","#F4E740", "#F4CB1C", "#19E80F",
                      "#566B44","#F90606","#B5EDB2","#AEAEAF",
                      "#95A0C6","#2C91E5",  "#7AEDD9","#C1A08A", "#2323E5","#9A9AE5")) + 
  theme(legend.title = element_text(size = 22), legend.text = element_text(size = 22)) + 
  theme(title = element_text(size = 25))+ theme(axis.text = element_text(size = 30)) +
  ggtitle("Annotation") 
p + ggsave("./figures/2.2/Annotated_umap.pdf", width = 15, height = 10)
p + ggsave("./figures/2.2/Annotated_umap.svg", width = 15, height = 10)


########## Plot marker genes of annotated clusters in Heatmap ##########
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
merged_cl_sub@active.ident <- factor(x = merged_cl_sub@active.ident, 
                                     levels = c("Hepatocytes","Metastasis","Cholangiocytes","Mac_C1q" ,"Mono_patrolling","Mac_Ly6c",
                                                "cDC1","cDC2","pDC",
                                                "LSECs", "LVECs",
                                                "T_CD4", "T_CD8","NKT",
                                                "NK","B_mem","B_plasma","Neutrophils",
                                                "Basophils","Fibroblasts","Stellate", "Kupffer","Kupffer_Endo"))
p <- DoHeatmap(subset(merged_cl_sub, downsample = 100), features = marker_genes, assay = "RNA",slot = "data") +
  theme(legend.title = element_text(size = 15), legend.text = element_text(size = 15)) + theme(title = element_text(size = 25))+ 
  theme(axis.text = element_text(size = 7)) +scale_fill_gradientn(colours = colorRampPalette(rev(brewer.pal(n = 15, name ="RdBu")))(256)) +
  ggtitle("Heatmap marker genes used for annotation")
p + ggsave("./figures/2.2/Heatmap_anno.png", width = 15, height = 10)
p + ggsave("./figures/2.2/Heatmap_anno.svg", width = 15, height = 10)

########## Add broad annotation ##########
Idents(merged_cl_sub) <- "annotation"
current.cluster.ids <- c("B_mem","B_plasma","Neutrophils","Basophils","cDC1","cDC2","Fibroblasts",
                         "Hepatocytes","Kupffer","Kupffer_Endo","LSECs","LVECs", "Metastasis",
                         "Mac_C1q", "Mono_patrolling","Cholangiocytes", "Mac_Ly6c","NK","NKT",
                         "pDC","Stellate", "T_CD4","T_CD8")
new.cluster.ids <- c("B","B","Neutrophils","Basophils","DCs","DCs","Fibroblasts",
                     "Hepatocytes","Kupffer","Kupffer_Endo","LECs","LECs", "Metastasis",
                     "Monocytes", "Monocytes","Cholangiocytes", "Monocytes","NK","T",
                     "DCs","Stellate", "T","T")
merged_cl_sub$annotation.broad <- plyr::mapvalues(x = merged_cl_sub$annotation, from = current.cluster.ids, to = new.cluster.ids)

p <- DimPlot(merged_cl_sub, label = TRUE,label.size = 5, group.by = "annotation.broad", pt.size = 0.5, 
             cols = c("#EF975B","#56595B","#EDB2D4","#91C6C4",
                      "#8E5229","#F46042","#E20FE8",
                      "#C491C6","#F4E740","#F90606","#19E80F",
                      "#AEAEAF" ,"#2C91E5", "#C1A08A","#2323E5")) + 
  theme(legend.title = element_text(size = 22), legend.text = element_text(size = 22)) +
  theme(title = element_text(size = 25))+ theme(axis.text = element_text(size = 30)) +
  ggtitle("Broad annotation") 
p + ggsave("./figures/2.2/Annotated_broad_umap.pdf", width = 15, height = 10)
p + ggsave("./figures/2.2/Annotated_broad_umap.svg", width = 15, height = 10)


########## Monocytes subtypes feature plot analysis ##########
#subset monocytes 
Idents(merged_cl_sub) <- "annotation.broad"
mono <- subset(merged_cl_sub, idents = "Monocytes")
Idents(mono) <- "annotation"

p <- FeaturePlot(mono, features = "Ly6c2", pt.size = 3) + theme(title = element_text(size = 25))+ theme(axis.text = element_text(size = 30)) 
p + ggsave("./figures/2.2/Feature_plot_Ly6c2.pdf", width = 12, height = 10)
p + ggsave("./figures/2.2/Feature_plot_Ly6c2.svg", width = 12, height = 10)

p <- FeaturePlot(mono, features = "C1qc", pt.size = 3) + theme(title = element_text(size = 25))+ theme(axis.text = element_text(size = 30)) 
p + ggsave("./figures/2.2/Feature_plot_C1qc.pdf", width = 12, height = 10)
p + ggsave("./figures/2.2/Feature_plot_C1qc.svg", width = 12, height = 10)

p <- FeaturePlot(mono, features = "Cx3cr1", pt.size = 3) + theme(title = element_text(size = 25))+ theme(axis.text = element_text(size = 30)) 
p + ggsave("./figures/2.2/Feature_plot_Cx3cr1.pdf", width = 12, height = 10)
p + ggsave("./figures/2.2/Feature_plot_Cx3cr1.svg", width = 12, height = 10)

p <- FeaturePlot(mono, features = "Spn", pt.size = 3) + theme(title = element_text(size = 25))+ theme(axis.text = element_text(size = 30)) 
p + ggsave("./figures/2.2/Feature_plot_Spn.pdf", width = 12, height = 10)
p + ggsave("./figures/2.2/Feature_plot_Spn.svg", width = 12, height = 10)

p <- FeaturePlot(mono, features = "Lyz2", pt.size = 3) + theme(title = element_text(size = 25))+ theme(axis.text = element_text(size = 30)) 
p + ggsave("./figures/2.2/Feature_plot_Lyz2.pdf", width = 12, height = 10)
p + ggsave("./figures/2.2/Feature_plot_Lyz2.svg", width = 12, height = 10)

p <- FeaturePlot(mono, features = "Tgfbi", pt.size = 3) + theme(title = element_text(size = 25))+ theme(axis.text = element_text(size = 30)) 
p + ggsave("./figures/2.2/Feature_plot_Tgfbi.pdf", width = 12, height = 10)
p + ggsave("./figures/2.2/Feature_plot_Tgfbi.svg", width = 12, height = 10)

########## Save annotated Seurat object ##########
saveRDS(merged_cl_sub,file = "./data_files_generated/LiverMerged_afterBC_anno.Rda")



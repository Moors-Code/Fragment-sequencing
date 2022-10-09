########## Part 1: Clustering, annotation and number of Features comparison with sphere-seq data ##########
#This part normalizes Visium data, merges them, annotates them based on the feature area and compares nFeatures with sphere-seq data 

########## Prepare environment ##########
###Setting the working directory 
setwd("/mnt/khandler/R_projects/Sphere-sequencing/Visium_analysis/")

###Load packages and functions 
source("./functions_and_packages/1.Packages.R")

###load R object 
V1 <- Load10X_Spatial("/home/khandler/NAS/Kristina/20220218_VisiumOutputs/Mets6M1/outs/")
V2 <- Load10X_Spatial("/home/khandler/NAS/Kristina/20220218_VisiumOutputs/Mets4M1/outs/")


V1$sample <- "V1"
V2$sample <- "V2"

########### nFeatures of spots ##########
SpatialFeaturePlot(V1, features = c("nFeature_Spatial")) + ggsave("./figures/1/V1_nFeature_Spatial.pdf")
SpatialFeaturePlot(V2, features = c("nFeature_Spatial")) + ggsave("./figures/1/V2_nFeature_Spatial.pdf")

########## Clustering and annotation of areas separately ##########
###Sample 1
V1 <- SCTransform(V1, assay = "Spatial", verbose = FALSE)
V1 <- RunPCA(V1, assay = "SCT", verbose = FALSE)
V1 <- FindNeighbors(V1, reduction = "pca", dims = 1:10)
V1 <- FindClusters(V1, verbose = FALSE, resolution = 0.2)
V1 <- RunUMAP(V1, reduction = "pca", dims = 1:10)
SpatialDimPlot(V1, label = TRUE, label.size = 3) 

#annotate clusters 
DotPlot(V1, features = c("Cyp2e1","Cyp1a2", "Alb", "Cyp2f2","Gpx2","Pglyrp1")) + 
  theme(legend.title = element_text(size = 22), legend.text = element_text(size = 22)) + 
  theme(title = element_text(size = 25))+ theme(axis.text = element_text(size = 30)) +
  theme(axis.text.x = element_text(angle = 90)) 

current.cluster.ids <- c(0:4)
new.cluster.ids <- c("Mets","CV","Mets","PV","Mixed")
V1$zone <- plyr::mapvalues(x = V1$seurat_clusters, from = current.cluster.ids, to = new.cluster.ids)
Idents(V1) <- "zone"
VlnPlot(V1, features = "nFeature_Spatial") 

V1 <- subset(V1, idents = c("Mets","CV","PV"))

current.cluster.ids <- c("Mets","CV","PV")
new.cluster.ids <- c("proximal","distal","distal")
V1$Mets_distance <- plyr::mapvalues(x = V1$zone, from = current.cluster.ids, to = new.cluster.ids)

Idents(V1) <- "Mets_distance"
V1_proximal <- subset(V1, idents = "proximal")
V1_distal <- subset(V1, idents = "distal")

###Sample 2 
V2 <- SCTransform(V2, assay = "Spatial", verbose = FALSE)
V2 <- RunPCA(V2, assay = "SCT", verbose = FALSE)
V2 <- FindNeighbors(V2, reduction = "pca", dims = 1:10)
V2 <- FindClusters(V2, verbose = FALSE, resolution = 0.2)
V2 <- RunUMAP(V2, reduction = "pca", dims = 1:10)
SpatialDimPlot(V2, label = TRUE, label.size = 3) 

#annotate clusters 
DotPlot(V2, features = c("Cyp2e1","Cyp1a2", "Alb", "Cyp2f2","Gpx2","Pglyrp1")) + 
  theme(legend.title = element_text(size = 22), legend.text = element_text(size = 22)) + 
  theme(title = element_text(size = 25))+ theme(axis.text = element_text(size = 30)) +
  theme(axis.text.x = element_text(angle = 90)) 

current.cluster.ids <- c(0:4)
new.cluster.ids <- c("Mets","Mixed","Mixed","CV","PV")
V2$zone <- plyr::mapvalues(x = V2$seurat_clusters, from = current.cluster.ids, to = new.cluster.ids)
Idents(V2) <- "zone"
VlnPlot(V2, features = "nFeature_Spatial") 

V2 <- subset(V2, idents = c("Mets","CV","PV"))

current.cluster.ids <- c("Mets","CV","PV")
new.cluster.ids <- c("proximal","distal","distal")
V2$Mets_distance <- plyr::mapvalues(x = V2$zone, from = current.cluster.ids, to = new.cluster.ids)

Idents(V2) <- "Mets_distance"
V2_proximal <- subset(V2, idents = "proximal")
V2_distal <- subset(V2, idents = "distal")

########## Plot number of Features per sample and feature area and compare with sphere-seq data ##########
###Visium sample 1
V1_nFeature_RNA_proximal <- median(V1_proximal@meta.data$nFeature_Spatial)
V1_nFeature_RNA_distal <- median(V1_distal@meta.data$nFeature_Spatial)

###Visium sample 2 
V2_nFeature_RNA_proximal <- median(V2_proximal@meta.data$nFeature_Spatial)
V2_nFeature_RNA_distal <- median(V2_distal@meta.data$nFeature_Spatial)

###get number of features of sphere-seq data
sphereSeq_mets <- readRDS("/mnt/khandler/R_projects/Sphere-sequencing/Sphere-seq_analysis/data_files_generated/LiverMerged_afterBC_anno_BS_5cells_zC_lobules_mets_distance.Rda")

Idents(sphereSeq_mets) <- "orig.ident"
SpS1 <- subset(sphereSeq_mets, idents = "6M1")
SpS2 <- subset(sphereSeq_mets, idents = "4M1")

Idents(SpS1) <- "Mets_distance"
proximal_SpS1 <- subset(SpS1, idents = "proximal") 
distal_SpS1 <- subset(SpS1, idents = "distal") 
median_nFeature_RNA_proximal_SpS1 <- median(proximal_SpS1@meta.data$nFeature_RNA)
median_nFeature_RNA_distal_SpS1 <- median(distal_SpS1@meta.data$nFeature_RNA)

Idents(SpS2) <- "Mets_distance"
proximal_SpS2 <- subset(SpS2, idents = "proximal") 
distal_SpS2 <- subset(SpS2, idents = "distal") 
median_nFeature_RNA_proximal_SpS2 <- median(proximal_SpS2@meta.data$nFeature_RNA)
median_nFeature_RNA_distal_SpS2 <- median(distal_SpS2@meta.data$nFeature_RNA)

###Plot nFeature per area and sample in barplot 
sample <- c("S1_SpS", "S1_Visium","S2_SpS","S2_Visium")
proximal_values <- c(median_nFeature_RNA_proximal_SpS1,V1_nFeature_RNA_proximal ,median_nFeature_RNA_proximal_SpS2,V2_nFeature_RNA_proximal)
distal_values <- c(median_nFeature_RNA_distal_SpS1, V1_nFeature_RNA_distal,median_nFeature_RNA_distal_SpS2,V2_nFeature_RNA_distal)

df_proximal <- data.frame(sample,proximal_values)
colnames(df_proximal) <- c("Sample","nFeature")
df_proximal$nFeature <- -(df_proximal$nFeature)
df_distal <- data.frame(sample,distal_values)
colnames(df_distal) <- c("Sample","nFeature")

df <- rbind(df_proximal,df_distal)

p <- ggplot(df, aes(x=nFeature, y=Sample,fill = nFeature)) + theme_classic() +
  geom_bar(stat = "identity",position = "identity",aes(fill = ifelse(nFeature<0,"proximal","distal"))) + 
  theme(axis.text = element_text(size = 15)) +
  theme(axis.title= element_text(size = 25)) + 
  ggtitle("nFeature of Visium and SpS experiments") + xlab("nFeatures") + 
  ylab("Sample")  +
  theme(plot.title = element_text(size = 15, face = "bold")) + theme(legend.text = element_text(size = 22),
                                                                     legend.title= element_text(size = 22)) 
p + ggsave("./figures/1/nFeatue_different_samples.pdf", width = 12, height = 7)
p + ggsave("./figures/1/nFeatue_different_samples.svg", width = 12, height = 7)

########## save R object ##########
saveRDS(V1,"./data_files_generated/Visium_sample1.rds")
saveRDS(V2,"./data_files_generated/Visium_sample2.rds")

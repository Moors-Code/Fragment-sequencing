########## Part 3: Spot deconvolution ##########
#This part does deconvolution of spots by integrating data with scRNAseq data of sphere-seq 

########## Prepare environment ##########
###Setting the working directory 
setwd("/mnt/khandler/R_projects/Sphere-sequencing/Visium_analysis/")

###Load packages and functions 
source("./functions_and_packages/1.Packages.R")

###load R object 
merged <- readRDS("./data_files_generated/Visium_merged_batch_corrected.rds")

########## apply another cutoff of nFeatures per spot ##########
#otherwise deconvolution is not working properly, in total spots are left with >300 nFeatures 
merged = merged[, merged$nFeature_Spatial > 200 ]

########## Prepare reference for deconvolution ##########
reference <- readRDS(file = "/mnt/khandler/R_projects/Sphere-sequencing/Sphere-seq_analysis/data_files_generated/LiverMerged_afterBC_anno.Rda")
#select 200 cells per subclass, fist set subclass as active.ident
Idents(reference) <- reference$annotation
reference <- subset(reference, cells = WhichCells(reference, downsample = 200))
#First run SCTransform and PCA
reference <- SCTransform(reference, ncells = 3000, verbose = FALSE, method = "poisson") %>%
  RunPCA(verbose = FALSE) %>%
  RunUMAP(dims = 1:30)

#the annotation is stored in the 'subclass' column of object metadata
DimPlot(reference, label = TRUE)

###Select genes for deconvolution 
reference@active.assay = "RNA"

markers_sc <- FindAllMarkers(reference, only.pos = TRUE, logfc.threshold = 0.1,
                             test.use = "wilcox", min.pct = 0.05, min.diff.pct = 0.1, max.cells.per.ident = 200,
                             return.thresh = 0.05, assay = "RNA")

#Filter for genes that are also present in the ST data
markers_sc <- markers_sc[markers_sc$gene %in% rownames(merged), ]

#Select top 20 genes per cluster, select top by first p-value, then absolute diff in pct, then quota of pct.
markers_sc$pct.diff <- markers_sc$pct.1 - markers_sc$pct.2
markers_sc$log.pct.diff <- log2((markers_sc$pct.1 * 99 + 1)/(markers_sc$pct.2 * 99 +
                                                               1))
markers_sc %>%
  group_by(cluster) %>%
  top_n(-100, p_val) %>%
  top_n(50, pct.diff) %>%
  top_n(20, log.pct.diff) -> top20
m_feats <- unique(as.character(top20$gene))

########## Deconvolution ##########
eset_SC <- ExpressionSet(assayData = as.matrix(reference@assays$RNA@counts[m_feats,]), 
                         phenoData = AnnotatedDataFrame(reference@meta.data))
eset_ST <- ExpressionSet(assayData = as.matrix(merged@assays$Spatial@counts[m_feats,]), 
                         phenoData = AnnotatedDataFrame(merged@meta.data))

###Deconvolve
deconvolution_crc <- SCDC::SCDC_prop(bulk.eset = eset_ST, sc.eset = eset_SC, ct.varname = "annotation",
                                     ct.sub = as.character(unique(eset_SC$annotation)))
head(deconvolution_crc$prop.est.mvw)

deconvolution_crc_df <- as.data.frame(deconvolution_crc$prop.est.mvw)
write.csv(deconvolution_crc_df, "./figures/3/deconvoluted_spots.csv")

#load data frame 
deconvolution_crc_df <- read.csv("./figures/3/deconvoluted_spots.csv")

########## annotate spots with cell type if at elast 75% are that cell type, others are mixed ##########
#extract cell-barcodes for individual cell types to then match it with Seurat object 
#assign a spot to a cell type if 75% of transcripts are assigned to one cell type 
cDC1 <- deconvolution_crc_df[deconvolution_crc_df$cDC1 >= 0.75,] #none
hepato <- deconvolution_crc_df[deconvolution_crc_df$Hepatocytes >= 0.75,] 
kc_endo <- deconvolution_crc_df[deconvolution_crc_df$Kupffer_Endo >= 0.75,] #none
Mono_pat <- deconvolution_crc_df[deconvolution_crc_df$Mono_patrolling >= 0.75,] #none
kc <- deconvolution_crc_df[deconvolution_crc_df$Kupffer >= 0.75,] #none
Baso <- deconvolution_crc_df[deconvolution_crc_df$Basophils >= 0.75,] #none
Mac_ly6c <- deconvolution_crc_df[deconvolution_crc_df$Mac_Ly6c >= 0.75,] 
NK <- deconvolution_crc_df[deconvolution_crc_df$NK >= 0.75,] #none
Neut <- deconvolution_crc_df[deconvolution_crc_df$Neutrophils >= 0.75,]
LSECs <- deconvolution_crc_df[deconvolution_crc_df$LSECs >= 0.75,]
cDC2 <- deconvolution_crc_df[deconvolution_crc_df$cDC2 >= 0.75,] #none
T_CD4 <- deconvolution_crc_df[deconvolution_crc_df$T_CD4 >= 0.75,] #none
pDC <- deconvolution_crc_df[deconvolution_crc_df$pDC >= 0.75,] #none
Mets <- deconvolution_crc_df[deconvolution_crc_df$Metastasis >= 0.75,] #non3
NKT <- deconvolution_crc_df[deconvolution_crc_df$NKT >= 0.75,] 
B_pl <- deconvolution_crc_df[deconvolution_crc_df$B_plasma >= 0.75,] 
T_CD8 <- deconvolution_crc_df[deconvolution_crc_df$T_CD8 >= 0.75,] #none
B_mem <- deconvolution_crc_df[deconvolution_crc_df$B_mem >= 0.75,] #none
Stellate <- deconvolution_crc_df[deconvolution_crc_df$Stellate >= 0.75,] 
LVECs <- deconvolution_crc_df[deconvolution_crc_df$LVECs >= 0.75,] #none
Chol <- deconvolution_crc_df[deconvolution_crc_df$Cholangiocytes >= 0.75,] 
Mac_C1q <- deconvolution_crc_df[deconvolution_crc_df$Mac_C1q >= 0.75,] #none
Fib <- deconvolution_crc_df[deconvolution_crc_df$Fibroblasts >= 0.75,] 

hepato <- hepato$X
Mac_ly6c <- Mac_ly6c$X
Neut <- Neut$X
LSECs <- LSECs$X
NKT <- NKT$X
B_pl <- B_pl$X
Stellate <- Stellate$X
Chol <- Chol$X
Fib <- Fib$X

all <- deconvolution_crc_df$X
mixed <- setdiff(all,c(hepato,Mac_ly6c,Neut,LSECs,NKT,B_pl,Stellate,Chol,Fib))

#calculate percentages for plotting 
percent.hepato <- length(hepato) *100 / length(all)
percent.mac_ly6c <- length(Mac_ly6c) *100 / length(all)
percent.neut <- length(Neut) *100 / length(all)
percent.lsecs <- length(LSECs) *100 / length(all)
percent.nkt <- length(NKT) *100 / length(all)
percent.Bpl <- length(B_pl) *100 / length(all)
percent.Stel <- length(Stellate) *100 / length(all)
percent.Chol <- length(Chol) *100 / length(all)
percent.Fib <- length(Fib) *100 / length(all)
percent.mixed <- length(mixed) *100 / length(all)
percent.others <- sum(percent.mac_ly6c,percent.neut,percent.lsecs,percent.nkt,percent.Bpl,percent.Stel,
                      percent.Fib)

########## Plot pie chart of percentage of cell types ##########
cell_types <- c("Hepatocytes","others",
                "Cholangiocytes","Mixed")
percentage <- c(percent.hepato,percent.others,percent.Chol,percent.mixed)
df <- data.frame(cell_types,percentage)

p <- ggplot(df,aes(x="",y=percentage,fill=cell_types)) + geom_bar(stat="identity",width=1,color = "white") + 
  coord_polar("y",start=0) + scale_fill_manual(values = c("#EDB2D4","#F46042","#2323E5",
                                                          "#56595B")) + 
  theme_void() 
p + ggsave("./figures/3/piechart.pdf",width = 12, height = 10)
p + ggsave("./figures/3/piechart.svg",width = 12, height = 10)




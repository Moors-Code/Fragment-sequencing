########## Part 9: Liver fragment_seq_vs_MC ##########
#This part compares cell type proportions between fragment-seq and MC data of CRC injected liver samples 

########## Prepare environment ##########
###Setting the working directory 
setwd("/mnt/khandler/R_projects/Fragment-sequencing/Fragment-seq_analysis/")

###Load packages and functions 
source("./functions_and_packages/1.Packages.R")

########## plot proportions in one plot ##########
#load R objects
MC <- readRDS("/mnt/khandler/R_projects/Fragment-sequencing/Highly_multiplexed_FISH_analysis/data_files_generated/Resolve_seurat_anno.rds")
SpS <- readRDS("./data_files_generated/LiverMerged_afterBC_anno.Rda")

#remove WT sample 
Idents(SpS) <- "orig.ident"
SpS <- subset(SpS, idents = c("M1","M3","M4","4M1","4M2","M5","6M1","6M2","6M3"))

#change Hepatocytes PV CV to only Hepatocytes in MC data 
current.cluster.ids <- c("B","Fibroblasts","Hepatocytes_CV","Hepatocytes_PV","Kupffer","LECs","Metastasis",
                         "Monocytes","Stellate","T")
new.cluster.ids <- c("B","Fibroblasts","Hepatocytes","Hepatocytes","Kupffer","LECs","Metastasis",
                     "Monocytes","Stellate","T")
MC$annotation <- plyr::mapvalues(x = MC$annotation, from = current.cluster.ids, to = new.cluster.ids)

#combine NK and T to T in SpS
current.cluster.ids <- c("B","Basophils","Cholangiocytes","DCs","Fibroblasts","Hepatocytes","Kupffer",
                         "Kupffer_Endo","LECs","Metastasis","Monocytes","Neutrophils","NK","Stellate","T")
new.cluster.ids <- c("B","Basophils","Cholangiocytes","DCs","Fibroblasts","Hepatocytes","Kupffer",
                     "Kupffer_Endo","LECs","Metastasis","Monocytes","Neutrophils","T","Stellate","T")
SpS$annotation.broad <- plyr::mapvalues(x = SpS$annotation.broad, from = current.cluster.ids, to = new.cluster.ids)

MC_tab<-  table(MC@meta.data[["annotation"]])
ann_tab <- t(MC_tab)
#add cell types that are missing 
Basophils <- c(0)
Cholangiocytes <- c(0)
DCs <- c(0)
Kupffer_Endo <- c(0)
Neutrophils <- c(0)

ann_tab <- cbind(ann_tab,Basophils)
ann_tab <- cbind(ann_tab,Cholangiocytes)
ann_tab <- cbind(ann_tab,DCs)
ann_tab <- cbind(ann_tab,Kupffer_Endo)
ann_tab <- cbind(ann_tab,Neutrophils)

ann_tab <- cbind(ann_tab, Total = rowSums(ann_tab))
ann_tab <- as.data.frame(ann_tab)
ann_tab_pct = lapply(ann_tab[,], function(x) {
  x/ann_tab$Total})
ann_tab_pct <- as.data.frame(ann_tab_pct)
colnames(ann_tab_pct) <- colnames(ann_tab)
ann_tab_pct$Total <- NULL

SpS_tab<-  table(SpS@meta.data[["annotation.broad"]])
SpS_tab <- t(SpS_tab)
SpS_tab <- cbind(SpS_tab, Total = rowSums(SpS_tab))
SpS_tab <- as.data.frame(SpS_tab)
SpS_tab_pct = lapply(SpS_tab[,], function(x) {
  x/SpS_tab$Total})
SpS_tab_pct <- as.data.frame(SpS_tab_pct)
colnames(SpS_tab_pct) <- colnames(SpS_tab)
SpS_tab_pct$Total <- NULL

merged_df <- full_join(SpS_tab_pct, ann_tab_pct)
rownames(merged_df) <- c("SpS","MC")

merged_df <- t(merged_df)

SpS <- merged_df[,"SpS"]
SpS <- as.data.frame(SpS)
SpS$sample <- "SpS"
colnames(SpS) <- c("proportion","sample")
SpS$cell_types <- rownames(SpS)

mc <- merged_df[,"MC"]
mc <- as.data.frame(mc)
mc$sample <- "MC"
colnames(mc) <- c("proportion","sample")
mc$cell_types <- rownames(mc)
#mc$proportion <- -mc$proportion

both <- rbind(mc,SpS)

#plot zoomed-in barplot 
p <- ggplot(both, aes(x=proportion, y=cell_types, groups = sample)) + theme_classic() +
  geom_bar(stat = "identity",position = "dodge",aes(colour = sample), fill =c("#EF975B","#56595B","#EDB2D4","#91C6C4",
                                                            "#8E5229","#F46042","#E20FE8",
                                                            "#C491C6","#F4E740","#F90606","#19E80F",
                                                            "#AEAEAF" , "#C1A08A","#2323E5",
                                                            "#EF975B","#56595B","#EDB2D4","#91C6C4",
                                                            "#8E5229","#F46042","#E20FE8",
                                                            "#C491C6","#F4E740","#F90606","#19E80F",
                                                            "#AEAEAF" , "#C1A08A","#2323E5") ) + 
  ###colors 
  #B cells (orange): B: #EF975B 
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
scale_color_manual(values = c("black","red")) +

theme(axis.text = element_text(size = 10))  + theme(axis.text.x = element_text(angle = 90)) +
  theme(axis.title= element_text(size = 10)) + 
  ggtitle("Cell type prop MC vs SpS") + xlab("Cell type proportions") + 
  ylab("Cell type")  + facet_zoom(xlim = c(0,0.08)) +
  theme(plot.title = element_text(size = 15, face = "bold")) + theme(legend.text = element_text(size = 22),
                                                                     legend.title= element_text(size = 22))  
p + ggsave("./figures/9/MC_sps_prop.pdf", width = 12, height = 10)
p + ggsave("./figures/9/MC_sps_prop.svg", width = 12, height = 10)


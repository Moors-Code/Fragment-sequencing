########## Part 2.4: Cells per sphere cutoff - Liver samples ##########
#This part applies a cutoff of at least 5 cells per sphere, all other spheres are not used in further analysis 

########## Prepare environment ##########
###Setting the working directory 
setwd("/mnt/khandler/R_projects/Sphere-sequencing/Sphere-seq_analysis/")

###Load packages and functions 
source("./functions_and_packages/1.Packages.R")
source("./functions_and_packages/4.Functions_cells_per_sphere_cutoff.R")

########## 5 cells per sphere cutoff ##########
Sphere_cell_cutoff("./data_files_generated/LiverMerged_afterBC_anno_BS.Rda",5,"./data_files_generated/LiverMerged_afterBC_anno_BS_5cells.Rda")

###load R object 
liverSpS5C <- readRDS("./data_files_generated/LiverMerged_afterBC_anno_BS_5cells.Rda")

########## Plot UMAP of spheres with at least 5 cells ##########
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

#fine grade annotation
p <- DimPlot(liverSpS5C, label = TRUE,label.size = 5, group.by = "annotation", pt.size = 0.5, 
             cols = c("#EF975B","#F4741E","#56595B","#91C6C4",
                      "#315B5A","#EDB2D4",
                      "#8E5229","#F46042","#E20FE8",
                      "#C491C6","#F4E740", "#F4CB1C", "#19E80F",
                      "#566B44","#F90606","#B5EDB2","#AEAEAF",
                      "#95A0C6","#2C91E5",  "#7AEDD9","#C1A08A", "#2323E5","#9A9AE5")) + 
  theme(legend.title = element_text(size = 22), legend.text = element_text(size = 22)) + theme(title = element_text(size = 25))+ 
  theme(axis.text = element_text(size = 30)) +
  ggtitle("Annotation - spheres ≥ 5 cells") 
p + ggsave("./figures/2.4/Annotated_umap5Cells.pdf", width = 15, height = 10)
p + ggsave("./figures/2.4/Annotated_umap5Cells.svg", width = 15, height = 10)

#broad annotation
p <- DimPlot(liverSpS5C, label = TRUE,label.size = 5, group.by = "annotation.broad", pt.size = 0.5, 
             cols = c("#EF975B","#56595B","#EDB2D4","#91C6C4",
                      "#8E5229","#F46042","#E20FE8",
                      "#C491C6","#F4E740","#F90606","#19E80F",
                      "#AEAEAF" ,"#2C91E5", "#C1A08A","#2323E5")) + 
  theme(legend.title = element_text(size = 22), legend.text = element_text(size = 22)) + theme(title = element_text(size = 25))+ 
  theme(axis.text = element_text(size = 30)) +
  ggtitle("Annotation - spheres ≥ 5 cells") 
p + ggsave("./figures/2.4/Annotated_umap5Cells_broad.pdf", width = 15, height = 10)
p + ggsave("./figures/2.4/Annotated_umap5Cells_broad.svg", width = 15, height = 10)

########## Plot barplot of cell type proportions per sphere ##########
###generate dataframe of cell type proportions per sphere 
df_sphere_cell_type <- table(liverSpS5C$sphere,liverSpS5C$annotation)
#include Total column
df_sphere_cell_type <- cbind(df_sphere_cell_type, Total = rowSums(df_sphere_cell_type))
df_sphere_cell_type <- as.data.frame(df_sphere_cell_type)
#calculate cell type proportions per sphere 
df_sphere_cell_type_pct = lapply(df_sphere_cell_type[,], function(x) {
  x/df_sphere_cell_type$Total
})
df_sphere_cell_type_pct <- as.data.frame(df_sphere_cell_type_pct)
rownames(df_sphere_cell_type_pct) <- rownames(df_sphere_cell_type)
#remove total column 
df_sphere_cell_type_pct <- df_sphere_cell_type_pct[,-ncol(df_sphere_cell_type_pct)]
df_sphere_cell_type_pct$sphere <- rownames(df_sphere_cell_type_pct)

###subset each celltype and add column name for cell type, then merge again 
ct1 <- df_sphere_cell_type_pct[,colnames(df_sphere_cell_type_pct) %in% c("B_mem","sphere")]
colnames(ct1) <- c("proportion","sphere")
ct1$annotation <- "B_mem"
ct2 <- df_sphere_cell_type_pct[,colnames(df_sphere_cell_type_pct) %in% c("B_plasma","sphere")]
colnames(ct2) <- c("proportion","sphere")
ct2$annotation <- "B_plasma"
ct3 <- df_sphere_cell_type_pct[,colnames(df_sphere_cell_type_pct) %in% c("Basophils","sphere")]
colnames(ct3) <- c("proportion","sphere")
ct3$annotation <- "Basophils"
ct4 <- df_sphere_cell_type_pct[,colnames(df_sphere_cell_type_pct) %in% c("cDC1","sphere")]
colnames(ct4) <- c("proportion","sphere")
ct4$annotation <- "cDC1"
ct5 <- df_sphere_cell_type_pct[,colnames(df_sphere_cell_type_pct) %in% c("cDC2","sphere")]
colnames(ct5) <- c("proportion","sphere")
ct5$annotation <- "cDC2"
ct6 <- df_sphere_cell_type_pct[,colnames(df_sphere_cell_type_pct) %in% c("Cholangiocytes","sphere")]
colnames(ct6) <- c("proportion","sphere")
ct6$annotation <- "Cholangiocytes"
ct7 <- df_sphere_cell_type_pct[,colnames(df_sphere_cell_type_pct) %in% c("Fibroblasts","sphere")]
colnames(ct7) <- c("proportion","sphere")
ct7$annotation <- "Fibroblasts"
ct8 <- df_sphere_cell_type_pct[,colnames(df_sphere_cell_type_pct) %in% c("Hepatocytes","sphere")]
colnames(ct8) <- c("proportion","sphere")
ct8$annotation <- "Hepatocytes"
ct9 <- df_sphere_cell_type_pct[,colnames(df_sphere_cell_type_pct) %in% c("Kupffer","sphere")]
colnames(ct9) <- c("proportion","sphere")
ct9$annotation <- "Kupffer"
ct10 <- df_sphere_cell_type_pct[,colnames(df_sphere_cell_type_pct) %in% c("Kupffer_Endo","sphere")]
colnames(ct10) <- c("proportion","sphere")
ct10$annotation <- "Kupffer_Endo"
ct11 <- df_sphere_cell_type_pct[,colnames(df_sphere_cell_type_pct) %in% c("LSECs","sphere")]
colnames(ct11) <- c("proportion","sphere")
ct11$annotation <- "LSECs"
ct12 <- df_sphere_cell_type_pct[,colnames(df_sphere_cell_type_pct) %in% c("LVECs","sphere")]
colnames(ct12) <- c("proportion","sphere")
ct12$annotation <- "LVECs"
ct13 <- df_sphere_cell_type_pct[,colnames(df_sphere_cell_type_pct) %in% c("Mac_C1q","sphere")]
colnames(ct13) <- c("proportion","sphere")
ct13$annotation <- "Mac_C1q"
ct14 <- df_sphere_cell_type_pct[,colnames(df_sphere_cell_type_pct) %in% c("Mac_Ly6c","sphere")]
colnames(ct14) <- c("proportion","sphere")
ct14$annotation <- "Mac_Ly6c"
ct15 <- df_sphere_cell_type_pct[,colnames(df_sphere_cell_type_pct) %in% c("Metastasis","sphere")]
colnames(ct15) <- c("proportion","sphere")
ct15$annotation <- "Metastasis"
ct16 <- df_sphere_cell_type_pct[,colnames(df_sphere_cell_type_pct) %in% c("Mono_patrolling","sphere")]
colnames(ct16) <- c("proportion","sphere")
ct16$annotation <- "Mono_patrolling"
ct17 <- df_sphere_cell_type_pct[,colnames(df_sphere_cell_type_pct) %in% c("Neutrophils","sphere")]
colnames(ct17) <- c("proportion","sphere")
ct17$annotation <- "Neutrophils"
ct18 <- df_sphere_cell_type_pct[,colnames(df_sphere_cell_type_pct) %in% c("NK","sphere")]
colnames(ct18) <- c("proportion","sphere")
ct18$annotation <- "NK"
ct19 <- df_sphere_cell_type_pct[,colnames(df_sphere_cell_type_pct) %in% c("NKT","sphere")]
colnames(ct19) <- c("proportion","sphere")
ct19$annotation <- "NKT"
ct20 <- df_sphere_cell_type_pct[,colnames(df_sphere_cell_type_pct) %in% c("pDC","sphere")]
colnames(ct20) <- c("proportion","sphere")
ct20$annotation <- "pDC"
ct21 <- df_sphere_cell_type_pct[,colnames(df_sphere_cell_type_pct) %in% c("Stellate","sphere")]
colnames(ct21) <- c("proportion","sphere")
ct21$annotation <- "Stellate"
ct22 <- df_sphere_cell_type_pct[,colnames(df_sphere_cell_type_pct) %in% c("T_CD4","sphere")]
colnames(ct22) <- c("proportion","sphere")
ct22$annotation <- "T_CD4"
ct23 <- df_sphere_cell_type_pct[,colnames(df_sphere_cell_type_pct) %in% c("T_CD8","sphere")]
colnames(ct23) <- c("proportion","sphere")
ct23$annotation <- "T_CD8"

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
ct <- rbind(ct,ct13)
ct <- rbind(ct,ct14)
ct <- rbind(ct,ct15)
ct <- rbind(ct,ct16)
ct <- rbind(ct,ct17)
ct <- rbind(ct,ct18)
ct <- rbind(ct,ct19)
ct <- rbind(ct,ct20)
ct <- rbind(ct,ct21)
ct <- rbind(ct,ct22)
ct <- rbind(ct,ct23)

###plot with colors to match cell types in UMAP 
p <- ggplot(ct, aes(fill=annotation, y=proportion, x=sphere)) + theme_classic() +
  geom_bar(position="stack", stat="identity" ) + 
  ggtitle("Cell type per sphere")+  theme(axis.text = element_text(size = 5)) + 
  theme(axis.title= element_text(size = 25))  +  theme(plot.title = element_text(size = 25, face = "bold")) + 
  theme(legend.title = element_text(size = 30), legend.text = element_text(size = 30))  + 
  xlab("Sphere") + ylab("Cell type proportion") + 
  theme(axis.text.x = element_text(angle = 90)) + theme(legend.position="none") + 
  scale_fill_manual(values=  c("#EF975B","#F4741E","#56595B","#91C6C4","#315B5A","#EDB2D4","#8E5229","#F46042",
                               "#E20FE8","#C491C6","#F4E740","#F4CB1C","#19E80F","#566B44","#F90606",
                               "#B5EDB2","#AEAEAF","#95A0C6",
                               "#2C91E5","#7AEDD9","#C1A08A","#2323E5","#9A9AE5"))
p + ggsave("./figures/2.4/Porportion_per_sphere.pdf", width = 15, height = 10)
p + ggsave("./figures/2.4/Porportion_per_sphere.svg", width = 15, height = 10)






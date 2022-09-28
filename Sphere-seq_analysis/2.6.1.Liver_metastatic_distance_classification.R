########## Part 2.6.1: Metastatic distance classification of liver samples ##########
#This part groups spheres into close and far distances to metastatic sites
#close = spheres with metastatic cells, far = spheres without metastatic cells 

########## Prepare environment ##########
###Setting the working directory 
setwd("/mnt/khandler/R_projects/Sphere-sequencing/Sphere-seq_analysis/")

###Load packages and functions 
source("./functions_and_packages/1.Packages.R")

###load R object
liverSpS5C <- readRDS("./data_files_generated/LiverMerged_afterBC_anno_BS_5cells_zC_lobules.Rda")

#Injectes samples  
Idents(liverSpS5C) <- "orig.ident"
Injected <- subset(liverSpS5C, idents = c("M1","M3","M4","4M1","4M2","M5","6M1","6M2","6M3"))

########## Identify samples with metastatic cells ##########
#only samples with at least 20 metastatic cells are considered in downstream analysis 
Idents(Injected) <- "annotation"
mets <- subset(Injected, idents = "Metastasis")
cell_counts_mets <- as.data.frame(table(mets$orig.ident))
p <- ggplot(cell_counts_mets, aes( y=Freq, x=Var1)) + theme_classic() +
  geom_bar( stat="identity", fill = "#F90606") + 
  xlab("Sample") + ylab("Cell count") +ggtitle("Amount of metastastatic cells in different samples") + 
  theme(axis.text = element_text(size = 30)) + 
  theme(axis.title= element_text(size = 25)) + 
  theme(plot.title = element_text(size = 25, face = "bold")) + 
  theme(legend.title = element_text(size = 30), legend.text = element_text(size = 30)) 
p + ggsave("./figures/2.6.1/barplot_mets_cells.pdf", width = 12, height = 10)
p + ggsave("./figures/2.6.1/barplot_mets_cells.svg", width = 12, height = 10)

########## Subset samples with at least 20 metastatic cells ##########
Idents(Injected) <- "orig.ident"
metastasis <- subset(Injected, idents = c("4M1","M5","6M1"))

########## Add close and far identity depending on the presence of metastatic cells #####
###subset metastatic cells and generate vector with sphere-BCs of spheres containing metastatic cells 
Idents(metastasis) <- "annotation"
sub_mets <- subset(metastasis, idents = c("Metastasis"))
mets_spheres <- as.data.frame(table(sub_mets$sphere))$Var1
all_spheres <- as.data.frame(table(metastasis$sphere))$Var1
no_mets_spheres <- setdiff(all_spheres,mets_spheres)

###add meta data with the identity of close and far #####
metastasis$Mets_distance <- NA
metastasis@meta.data <- metastasis@meta.data %>%
  mutate(Mets_distance = case_when(
    sphere %in% mets_spheres ~ "close",
    sphere %in% no_mets_spheres ~ "far",
    TRUE ~ NA_character_))

###save object with new meta-data column 
saveRDS(metastasis,"./data_files_generated/LiverMerged_afterBC_anno_BS_5cells_zC_lobules_mets_distance.Rda")

########## Make split plot of cells within close and far area ##########
###colours 
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

p <- DimPlot(metastasis, reduction = "umap", label = TRUE, pt.size = 0.5, group.by = "annotation",split.by = "Mets_distance",
             cols = c("#EF975B","#F4741E","#56595B","#91C6C4",
                      "#315B5A","#EDB2D4",
                      "#8E5229","#F46042","#E20FE8",
                      "#C491C6","#F4E740", "#F4CB1C", "#19E80F",
                      "#566B44","#F90606","#B5EDB2","#AEAEAF",
                      "#95A0C6","#2C91E5",  "#7AEDD9","#C1A08A", "#2323E5","#9A9AE5")) + 
  theme(title = element_text(size = 25))+ theme(axis.text = element_text(size = 30)) +
  ggtitle("Annotation")
p + ggsave("./figures/2.6.1/Split_plot_mets_noMets.pdf", width = 15, height = 10)
p + ggsave("./figures/2.6.1/Split_plot_mets_noMets.svg", width = 15, height = 10)










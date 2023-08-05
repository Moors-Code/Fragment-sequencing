########## Part 2.6.1: Metastatic distance classification of liver samples ##########
#This part groups fragments into proximal and distal distances to metastatic sites
#proximal = fragments with metastatic cells, distal = fragments without metastatic cells 

########## Prepare environment ##########
###Setting the working directory 
setwd("/mnt/khandler/R_projects/Fragment-sequencing/Fragment-seq_analysis/")

###Load packages and functions 
source("./functions_and_packages/1.Packages.R")

###load R object
liverSpS5C <- readRDS("./data_files_generated/LiverMerged_afterBC_anno_BS_5cells_zC_lobules.Rda")

#Injected samples  
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

########## Add proximal and distal identity depending on the presence of metastatic cells ##########
###subset metastatic cells and generate vector with fragment-BCs of fragments containing metastatic cells 
Idents(metastasis) <- "annotation"
sub_mets <- subset(metastasis, idents = c("Metastasis"))
mets_fragments <- as.data.frame(table(sub_mets$fragment))$Var1
all_fragments <- as.data.frame(table(metastasis$fragment))$Var1
no_mets_fragments <- setdiff(all_fragments,mets_fragments)

###add meta data with the identity of proximal and distal 
metastasis$Mets_distance <- NA
metastasis@meta.data <- metastasis@meta.data %>%
  mutate(Mets_distance = case_when(
    fragment %in% mets_fragments ~ "proximal",
    fragment %in% no_mets_fragments ~ "distal",
    TRUE ~ NA_character_))

###save object with new meta-data column 
saveRDS(metastasis,"./data_files_generated/LiverMerged_afterBC_anno_BS_5cells_zC_lobules_mets_distance.Rda")

########## Make split plot of cells within distal and proximal area ##########
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



########## Metastatic cells per fragment ##########
###subset proximal area 
Idents(metastasis) <- "Mets_distance"
mets_fragments <- subset(metastasis,idents = "proximal")

###subset metastasis cells from proximal area
Idents(mets_fragments) <- "annotation"
mets_fragments_mets_cells <- subset(mets_fragments, idents = "Metastasis")

###plot counts in boxplot in decreasing order
mets_fragments_mets_cells_df <- as.data.frame(table(mets_fragments_mets_cells$fragment))
p <- ggplot(mets_fragments_mets_cells_df, aes( y=Freq, x=reorder(Var1,-Freq), Var1)) + theme_classic() +
  geom_bar( stat="identity", fill = "#F90606") + 
  xlab("fragment") + ylab("Cell count") +ggtitle("Amount of metastastatic cells per fragment") + 
  theme(axis.text = element_text(size = 10)) +  theme(axis.text.x = element_text(angle = 90)) +
  theme(axis.title= element_text(size = 25)) + 
  theme(plot.title = element_text(size = 25, face = "bold")) + 
  theme(legend.title = element_text(size = 30), legend.text = element_text(size = 30)) 
p + ggsave("./figures/2.6.1/barplot_mets_cells_per_fragment.pdf", width = 12, height = 10)
p + ggsave("./figures/2.6.1/barplot_mets_cells_per_fragment.svg", width = 12, height = 10)

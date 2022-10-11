########## Part 2.6.2: Cell type abundance between proximal and distal distances from metastatic sites##########
#This part compares cell type abundance of proximal and distal distances from metastatic sites 

########## Prepare environment ##########
###Setting the working directory 
setwd("/mnt/khandler/R_projects/Sphere-sequencing/Sphere-seq_analysis/")

###Load packages and functions 
source("./functions_and_packages/1.Packages.R")
source("./functions_and_packages/6.Functions_Cell_type_abundance.R")

###load R object
metastasis <- readRDS("./data_files_generated/LiverMerged_afterBC_anno_BS_5cells_zC_lobules_mets_distance.Rda")

########## Differential abundance analysis of all cell types ##########
###run DA analysis 
DA_analysis_cell_type_abundance_Mets_distance(metastasis,metastasis$annotation.broad,"./figures/2.6.2/","all_cells")

###plot results of DA analysis in vulcano plot 
top <- read.csv("./figures/2.6.2/Mets_distance_all_cells_cell_type_abundance_Sphere_seq.csv")

top <- top %>% 
  mutate(
    Expression = case_when(logFC >= 0.3 & FDR <= 0.05 ~ "High in proximal",
                           logFC <= -0.3 & FDR <= 0.05 ~ "High in distal",
                           TRUE ~ "Non sig.")
  )

p <- ggplot(top, aes(x=logFC, y=-log10(FDR))) +
  geom_point(aes(color = Expression),size=5) +
  geom_text(data=top[top$FDR<1 & abs(top$logFC) > 0,], aes(label=Gene),size=8) +
  xlab("logFC") + 
  ylab("-log10(FDR)") + ggtitle("Cell type prop - Sphere-seq (FDR ≤ 0.05, logFC >0.3") + 
  scale_color_manual(values = c("firebrick3","dodgerblue3", "gray50"),guide = "none") + theme_classic() + 
  theme(axis.title= element_text(size = 25)) + theme(axis.text = element_text(size = 30))  + 
  theme(plot.title = element_text(size = 25, face = "bold"))  + 
  guides(colour = guide_legend(override.aes = list(size=1.5))) + 
  theme(legend.title = element_text(size = 15), legend.text = element_text(size = 15))  + 
  geom_vline(xintercept = c(-0.3,0.3), linetype = "dashed", color = "gray50") 
p + ggsave("./figures/2.6.2/Cell_type_abundance_mets_distance_allCells_sphereSeq.pdf",width = 12, height = 10)
p + ggsave("./figures/2.6.2/Cell_type_abundance_mets_distance_allCells_sphereSeq.svg",width = 12, height = 10)

########## Differential abundance analysis of monocytes subtypes ##########
###subset monocytes
Idents(metastasis) <- "annotation.broad"
mono <- subset(metastasis, idents = "Monocytes")

###run DA analysis 
DA_analysis_cell_type_abundance_Mets_distance(mono,mono$annotation,"./figures/2.6.2/","mono_subtypes")

###read result for p-values in boxplot 
top <- read.csv("./figures/2.6.2/Mets_distance_mono_subtypes_cell_type_abundance_Sphere_seq.csv")

########## Plot proportions in boxplots proximal vs. distal of cell types ##########
###create table of cell type proportions per sample 
Idents(metastasis) <- "orig.ident"
S4 <- subset(metastasis, idents = "4M1")
S6 <- subset(metastasis, idents = "M5")
S7 <- subset(metastasis, idents = "6M1")

create_table_cell_type_prop(S4,"Mets_distance", "annotation.broad","./figures/2.6.2/","S4")
create_table_cell_type_prop(S6,"Mets_distance", "annotation.broad","./figures/2.6.2/","S6")
create_table_cell_type_prop(S7,"Mets_distance", "annotation.broad","./figures/2.6.2/","S7")

#read csv files of cell type proportions per sample and generate one data frame with all 
df_S4 <- read.csv("./figures/2.6.2/S4_proportions_Mets_distance_annotation.broad.csv")
df_S6 <- read.csv("./figures/2.6.2/S6_proportions_Mets_distance_annotation.broad.csv")
df_S7 <- read.csv("./figures/2.6.2/S7_proportions_Mets_distance_annotation.broad.csv")

df_S4$sample <- "S4"
df_S6$sample <- "S6"
df_S7$sample <- "S7"

df_all <- rbind.fill(df_S4, df_S6)
df_all <- rbind.fill(df_all, df_S7)

###make boxplot with appropriate color
#Kupffer cells
boxplot_cell_prop(df_all, "Kupffer","#E20FE8","./figures/2.6.2/",".pdf")
boxplot_cell_prop(df_all, "Kupffer","#E20FE8","./figures/2.6.2/",".svg")

#LECs
boxplot_cell_prop(df_all, "LECs","#F4E740","./figures/2.6.2/",".pdf")
boxplot_cell_prop(df_all, "LECs","#F4E740","./figures/2.6.2/",".svg")

#Monocytes
boxplot_cell_prop(df_all, "Monocytes","#19E80F","./figures/2.6.2/",".pdf")
boxplot_cell_prop(df_all, "Monocytes","#19E80F","./figures/2.6.2/",".svg")

#Metastatic cells
boxplot_cell_prop(df_all, "Metastasis","#F90606","./figures/2.6.2/",".pdf")
boxplot_cell_prop(df_all, "Metastasis","#F90606","./figures/2.6.2/",".svg")

########## Plot proportions in boxplots proximal vs. distal of monocytes subtypes ##########
###create table of cell type proportions per sample 
Idents(mono) <- "orig.ident"
S4 <- subset(mono, idents = "4M1")
S6 <- subset(mono, idents = "M5")
S7 <- subset(mono, idents = "6M1")

create_table_cell_type_prop(S4, "Mets_distance","annotation","./figures/2.6.2/","S4_mono")
create_table_cell_type_prop(S6, "Mets_distance","annotation","./figures/2.6.2/","S6_mono")
create_table_cell_type_prop(S7, "Mets_distance","annotation","./figures/2.6.2/","S7_mono")

#read csv files of cell type proportions per sample and generate one data frame with all 
df_S4 <- read.csv("./figures/2.6.2/S4_mono_proportions_Mets_distance_annotation.csv")
df_S6 <- read.csv("./figures/2.6.2/S6_mono_proportions_Mets_distance_annotation.csv")
df_S7 <- read.csv("./figures/2.6.2/S7_mono_proportions_Mets_distance_annotation.csv")

df_S4$sample <- "S4"
df_S6$sample <- "S6"
df_S7$sample <- "S7"

df_all <- rbind.fill(df_S4, df_S6)
df_all <- rbind.fill(df_all, df_S7)

###make boxplot with appropriate color
#C1q+ monocytes
boxplot_cell_prop(df_all, "Mac_C1q","#19E80F","./figures/2.6.2/",".pdf")
boxplot_cell_prop(df_all, "Mac_C1q","#19E80F","./figures/2.6.2/",".svg")

#Ly6c+ monocytes
boxplot_cell_prop(df_all, "Mac_Ly6c","#566B44","./figures/2.6.2/",".pdf")
boxplot_cell_prop(df_all, "Mac_Ly6c","#566B44","./figures/2.6.2/",".svg")

#Patrolling monocytes
boxplot_cell_prop(df_all, "Mono_patrolling","#B5EDB2","./figures/2.6.2/",".pdf")
boxplot_cell_prop(df_all, "Mono_patrolling","#B5EDB2","./figures/2.6.2/",".svg")




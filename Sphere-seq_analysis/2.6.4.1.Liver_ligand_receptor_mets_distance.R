########## Part 2.6.4.1: Ligand-receptor interaction analysis between proximal and distal areas from metastatic sites ##########
#This part first generates input files for CellPhoneDB for proximal and distal areas
#Then significant interactions are compared between areas to find specific ones 

########## Prepare environment ##########
###Setting the working directory 
setwd("/mnt/khandler/R_projects/Sphere-sequencing/Sphere-seq_analysis/")

###Load packages and functions 
source("./functions_and_packages/1.Packages.R")
source("./functions_and_packages/7.Functions_ligand_receptor_analysis.R")

###load R object
metastasis <- readRDS("./data_files_generated/LiverMerged_afterBC_anno_BS_5cells_zC_lobules_mets_distance.Rda")

########## CellPhoneDB input file generation per Mets distance and per mouse ##########
#files are then used in code 2.6.4.2.CellPhoneDB_mets_distance.sh
#we need CellPhoneDB analysis of proximal and distal over all samples and of proximal/distal from each individual mouse for statistical analysis 

#load human and mouse ensemble symbols
human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/") 
mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")

#subset samples and Mets distances 
Idents(metastasis) <- "Mets_distance"
proximal <- subset(metastasis, idents = "proximal")
distal <- subset(metastasis, idents = "distal")

Idents(proximal) <- "orig.ident"
proximal_S4 <- subset(proximal, idents = "4M1")
proximal_S6 <- subset(proximal, idents = "M5")
proximal_S7 <- subset(proximal, idents = "6M1")

Idents(distal) <- "orig.ident"
distal_S4 <- subset(distal, idents = "4M1")
distal_S6 <- subset(distal, idents = "M5")
distal_S7 <- subset(distal, idents = "6M1")

#input file generation 
Input_files_CellPhoneDB_generation(proximal,'annotation.broad',"proximal","./CellPhoneDB/Mets_distance/") 
Input_files_CellPhoneDB_generation(proximal_S4,'annotation.broad',"proximal_S4","./CellPhoneDB/Mets_distance/") 
Input_files_CellPhoneDB_generation(proximal_S6,'annotation.broad',"proximal_S6","./CellPhoneDB/Mets_distance/") 
Input_files_CellPhoneDB_generation(proximal_S7,'annotation.broad',"proximal_S7","./CellPhoneDB/Mets_distance/") 


Input_files_CellPhoneDB_generation(distal,'annotation.broad',"distal","./CellPhoneDB/Mets_distance/") 
Input_files_CellPhoneDB_generation(distal_S4,'annotation.broad',"distal_S4","./CellPhoneDB/Mets_distance/") 
Input_files_CellPhoneDB_generation(distal_S6,'annotation.broad',"distal_S6","./CellPhoneDB/Mets_distance/") 
Input_files_CellPhoneDB_generation(distal_S7,'annotation.broad',"distal_S7","./CellPhoneDB/Mets_distance/") 

########## Run CellPhoneDB_mets_distance.sh ##########
#/mnt/khandler/R_projects/Sphere-sequencing/Sphere-seq_analysis/2.6.4.2.CellPHoneDB_mets_distance.sh

########## Comparison of L-R interactions between proximal and distal areas of Monocytes and T cells ##########
###Matrix generation of interacting cells of interest 
Score_matrix_generation_interacting_pair_oi_per_mouse(
  "./CellPhoneDB/Mets_distance/",
  c("proximal_S4","proximal_S6","proximal_S7"),"Monocytes.T","./CellPhoneDB/Mets_distance/","proximal") 
Score_matrix_generation_interacting_pair_oi_per_mouse(
  "./CellPhoneDB/Mets_distance/",
  c("distal_S4","distal_S6","distal_S7"),"Monocytes.T","./CellPhoneDB/Mets_distance/","distal") 

###DE analysis between proximal and distal of interacting cells of interest 
DE_CellPhoneDB("./CellPhoneDB/Mets_distance/proximal_per_mouse_Monocytes.T.csv","./CellPhoneDB/Mets_distance/distal_per_mouse_Monocytes.T.csv",
               c("proximal_S4","proximal_S6","proximal_S7"),
               c("far_S4","far_S6","far_S7"),
               "proximal","distal","./CellPhoneDB/Mets_distance/","Monocytes.T","mets_distance")

########## Plotting of DE analysis in Barplot with proximal on the left side and distal on the right side ##########
#include Pvalue from DE analysis
top <- read.csv("./CellPhoneDB/Mets_distance/Monocytes.T_mets_distance_comp_LM.csv")
top <- top[, colnames(top) %in% c("X","PValue")]
colnames(top) <- c("interacting_pair","PValue")

###compare sig scores of proximal and distal 
file1 <- "./CellPhoneDB/Mets_distance/proximal/significant_means.txt"
file2 <- "./CellPhoneDB/Mets_distance/distal/significant_means.txt"

L_R_CellPhoneDB_comp_2samples("Monocytes.T",file1,file2,"./CellPhoneDB/Mets_distance/", "sig_means_mets_distance")

df <- read.csv("./CellPhoneDB/Mets_distance/interactions_comp_two_cond_Monocytes.T_sig_means_mets_distance.csv")

df$X <- NULL
colnames(df) <- c("interacting_pair","proximal","distal")

df_proximal <- df[,c("interacting_pair","proximal")]
df_proximal$Mets_distance <- "proximal"
colnames(df_proximal) <- c("interacting_pair","interaction_score","condition")
#make minus to put on one axis, has to be changed to + scores later in Illustrator 
df_proximal$interaction_score <- -(df_proximal$interaction_score)
df_proximal2 <- merge(df_proximal,top)

df_distal <- df[,c("interacting_pair","distal")]
df_distal$Mets_distance <- "far"
colnames(df_distal) <- c("interacting_pair","interaction_score","condition")
df_distal2 <- merge(df_distal,top)

df1 <- rbind(df_proximal2,df_distal2)

##make barplot 
#plot until p val 0.5
df2 <- df1[df1$PValue <=0.5,]

p <- ggplot(df2, aes(x=interaction_score, y=reorder(interacting_pair,-PValue),fill = PValue)) + theme_classic() +
  geom_bar(stat = "identity",position = "identity") + 
  scale_fill_gradientn(colours = c("darkblue","grey","darkred"),limits = c(0,0.5),breaks = c(0,0.05,0.1,0.2,0.3,0.5)) + 
  theme(axis.text = element_text(size = 15))  + theme(axis.text.x = element_text(angle = 90)) +
  theme(axis.title= element_text(size = 25)) + 
  ggtitle("L-R interactions - CellPhoneDB (Monocytes|T)") + xlab("Significant interaction score") + 
  ylab("L-R interaction")  +
  theme(plot.title = element_text(size = 15, face = "bold")) + theme(legend.text = element_text(size = 22),
                                                                     legend.title= element_text(size = 22)) + 
  guides(fill=guide_legend(title="P-value of enrichment")) 
p + ggsave("./figures/2.6.4/Interaction_Mono_T_Mets_distance_0.5.pdf", width = 12, height = 10)
p + ggsave("./figures/2.6.4/Interaction_Mono_T_Mets_distance_0.5.svg", width = 12, height = 10)




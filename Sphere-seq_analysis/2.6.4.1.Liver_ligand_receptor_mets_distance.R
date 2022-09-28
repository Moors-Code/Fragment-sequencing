########## Part 2.6.4.1: Ligand-receptor interaction analysis between close and far areas from metastatic sites ##########
#This part first generates input files for CellPhoneDB for close and far areas
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
#files are then used in code 2.6.4.2.CellPhoneDB_veins.sh
#we need CellPhoneDB analysis of close and far over all samples and of close/far from each individual mouse for statistical analysis 

#load human and mouse ensemble symbols
human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/") 
mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")

#subset samples and Mets distances 
Idents(metastasis) <- "Mets_distance"
close <- subset(metastasis, idents = "close")
far <- subset(metastasis, idents = "far")

Idents(close) <- "orig.ident"
close_S4 <- subset(cv, idents = "4M1")
close_S6 <- subset(cv, idents = "M5")
close_S7 <- subset(cv, idents = "6M1")

Idents(far) <- "orig.ident"
far_S4 <- subset(far, idents = "4M1")
far_S6 <- subset(far, idents = "M5")
far_S7 <- subset(far, idents = "6M1")

#input file generation 
Input_files_CellPhoneDB_generation(close,'annotation.broad',"close","./CellPhoneDB/Mets_distance/") 
Input_files_CellPhoneDB_generation(close_S4,'annotation.broad',"close_S4","./CellPhoneDB/Mets_distance/") 
Input_files_CellPhoneDB_generation(close_S6,'annotation.broad',"close_S6","./CellPhoneDB/Mets_distance/") 
Input_files_CellPhoneDB_generation(close_S7,'annotation.broad',"close_S7","./CellPhoneDB/Mets_distance/") 


Input_files_CellPhoneDB_generation(far,'annotation.broad',"far","./CellPhoneDB/Mets_distance/") 
Input_files_CellPhoneDB_generation(far_S4,'annotation.broad',"far_S4","./CellPhoneDB/Mets_distance/") 
Input_files_CellPhoneDB_generation(far_S6,'annotation.broad',"far_S6","./CellPhoneDB/Mets_distance/") 
Input_files_CellPhoneDB_generation(far_S7,'annotation.broad',"far_S7","./CellPhoneDB/Mets_distance/") 


########## Comparison of L-R interactions between close and far areas of Monocytes and T cells ##########
###Matrix generation of interacting cells of interest 
Score_matrix_generation_interacting_pair_oi_per_mouse(
  "./CellPhoneDB/Mets_distance/",
  c("close_S4","close_S6","close_S7"),"Monocytes.T","./CellPhoneDB/Mets_distance/","close") 
Score_matrix_generation_interacting_pair_oi_per_mouse(
  "./CellPhoneDB/Mets_distance/",
  c("far_S4","far_S6","far_S7"),"Monocytes.T","./CellPhoneDB/Mets_distance/","far") 

###DE analysis between close and far of interacting cells of interest 
DE_CellPhoneDB("./CellPhoneDB/Mets_distance/close_per_mouse_Monocytes.T.csv","./CellPhoneDB/Mets_distance/far_per_mouse_Monocytes.T.csv",
               c("close_S4","close_S6","close_S7"),
               c("far_S4","far_S6","far_S7"),
               "close","far","./CellPhoneDB/Mets_distance/","Monocytes.T","mets_distance")

########## Plotting of DE analysis in Barplot with close on the left side and far on the right side ##########
#include Pvalue from DE analysis
top <- read.csv("./CellPhoneDB/Mets_distance/Monocytes.T_mets_distance_comp_LM.csv")
top <- top[, colnames(top) %in% c("X","PValue")]
colnames(top) <- c("interacting_pair","PValue")

###compare sig scores of close and far 
file1 <- "./CellPhoneDB/Mets_distance/close/significant_means.txt"
file2 <- "./CellPhoneDB/Mets_distance/far/significant_means.txt"

L_R_CellPhoneDB_comp_2samples("Monocytes.T",file1,file2,"./CellPhoneDB/Mets_distance/", "sig_means_mets_distance")

df <- read.csv("./CellPhoneDB/Mets_distance/interactions_comp_two_cond_Monocytes.T_sig_means_mets_distance.csv")

df$X <- NULL
colnames(df) <- c("interacting_pair","close","far")

df_close <- df[,c("interacting_pair","close")]
df_close$Mets_distance <- "close"
colnames(df_close) <- c("interacting_pair","interaction_score","condition")
#make minus to put on one axis, has to be changed to + scores later in Illustrator 
df_close$interaction_score <- -(df_close$interaction_score)
df_close2 <- merge(df_close,top)

df_far <- df[,c("interacting_pair","far")]
df_far$Mets_distance <- "far"
colnames(df_far) <- c("interacting_pair","interaction_score","condition")
df_far2 <- merge(df_far,top)

df1 <- rbind(df_close2,df_far2)

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




########## Part 2.5.3.1: Ligand-receptor interaction analysis between central and portal veins of injected samples ##########
#This part first generates input files for CellPhoneDB for CV and PV vein areas 
#Then significant interactions are compared between vein areas to find specific ones 

########## Prepare environment ##########
###Setting the working directory 
setwd("/mnt/khandler/R_projects/Sphere-sequencing/Sphere-seq_analysis/")

###Load packages and functions 
source("./functions_and_packages/1.Packages.R")
source("./functions_and_packages/7.Functions_ligand_receptor_analysis.R")

###load R object
liverSpS5C <- readRDS("./data_files_generated/LiverMerged_afterBC_anno_BS_5cells_zC_lobules.Rda")

#Injected samples  
Idents(liverSpS5C) <- "orig.ident"
Injected <- subset(liverSpS5C, idents = c("M1","M3","M4","4M1","4M2","M5","6M1","6M2","6M3"))

########## CellPhoneDB input file generation per vein and per mouse ##########
#files are then used in code 2.5.3.1.CellPhoneDB_veins.sh
#we need CellPhoneDB analysis of CV and PV over all samples and of CV/PV from each individual mouse for statistical analysis 

#load human and mouse ensemble symbols
human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/") 
mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")

#subset samples and veins 
Idents(Injected) <- "vein"
cv <- subset(Injected, idents = "CV")
pv <- subset(Injected, idents = "PV")

Idents(cv) <- "orig.ident"
cv_S1 <- subset(cv, idents = "M1")
cv_S2 <- subset(cv, idents = "M3")
cv_S3 <- subset(cv, idents = "M4")
cv_S4 <- subset(cv, idents = "4M1")
cv_S5 <- subset(cv, idents = "4M2")
cv_S6 <- subset(cv, idents = "M5")
cv_S7 <- subset(cv, idents = "6M1")
cv_S8 <- subset(cv, idents = "6M2")
cv_S9 <- subset(cv, idents = "6M3")

Idents(pv) <- "orig.ident"
pv_S1 <- subset(pv, idents = "M1")
pv_S2 <- subset(pv, idents = "M3")
pv_S3 <- subset(pv, idents = "M4")
pv_S4 <- subset(pv, idents = "4M1")
pv_S5 <- subset(pv, idents = "4M2")
pv_S6 <- subset(pv, idents = "M5")
pv_S7 <- subset(pv, idents = "6M1")
pv_S8 <- subset(pv, idents = "6M2")
pv_S9 <- subset(pv, idents = "6M3")

#input file generation 
Input_files_CellPhoneDB_generation(cv,'annotation.broad',"cv","./CellPhoneDB/Vein/") 
Input_files_CellPhoneDB_generation(cv_S1,'annotation.broad',"cv_S1","./CellPhoneDB/Vein/") 
Input_files_CellPhoneDB_generation(cv_S2,'annotation.broad',"cv_S2","./CellPhoneDB/Vein/") 
Input_files_CellPhoneDB_generation(cv_S3,'annotation.broad',"cv_S3","./CellPhoneDB/Vein/") 
Input_files_CellPhoneDB_generation(cv_S4,'annotation.broad',"cv_S4","./CellPhoneDB/Vein/") 
Input_files_CellPhoneDB_generation(cv_S5,'annotation.broad',"cv_S5","./CellPhoneDB/Vein/") 
Input_files_CellPhoneDB_generation(cv_S6,'annotation.broad',"cv_S6","./CellPhoneDB/Vein/") 
Input_files_CellPhoneDB_generation(cv_S7,'annotation.broad',"cv_S7","./CellPhoneDB/Vein/") 
Input_files_CellPhoneDB_generation(cv_S8,'annotation.broad',"cv_S8","./CellPhoneDB/Vein/") 
Input_files_CellPhoneDB_generation(cv_S9,'annotation.broad',"cv_S9","./CellPhoneDB/Vein/") 

Input_files_CellPhoneDB_generation(pv,'annotation.broad',"pv","./CellPhoneDB/Vein/") 
Input_files_CellPhoneDB_generation(pv_S1,'annotation.broad',"pv_S1","./CellPhoneDB/Vein/") 
Input_files_CellPhoneDB_generation(pv_S2,'annotation.broad',"pv_S2","./CellPhoneDB/Vein/") 
Input_files_CellPhoneDB_generation(pv_S3,'annotation.broad',"pv_S3","./CellPhoneDB/Vein/") 
Input_files_CellPhoneDB_generation(pv_S4,'annotation.broad',"pv_S4","./CellPhoneDB/Vein/") 
Input_files_CellPhoneDB_generation(pv_S5,'annotation.broad',"pv_S5","./CellPhoneDB/Vein/") 
Input_files_CellPhoneDB_generation(pv_S6,'annotation.broad',"pv_S6","./CellPhoneDB/Vein/") 
Input_files_CellPhoneDB_generation(pv_S7,'annotation.broad',"pv_S7","./CellPhoneDB/Vein/") 
Input_files_CellPhoneDB_generation(pv_S8,'annotation.broad',"pv_S8","./CellPhoneDB/Vein/") 
Input_files_CellPhoneDB_generation(pv_S9,'annotation.broad',"pv_S9","./CellPhoneDB/Vein/") 

########## Run CellPhoneDB_veins.sh ##########
#/mnt/khandler/R_projects/Sphere-sequencing/Sphere-seq_analysis/2.5.3.2.CellPHoneDB_veins.sh

########## Comparison of L-R interactions between CV and PV of Kupffer cells and T cells ##########
###Matrix generation of interacting cells of interest 
Score_matrix_generation_interacting_pair_oi_per_mouse(
  "./CellPhoneDB/Vein/",
  c("cv_S1","cv_S2","cv_S3","cv_S4","cv_S5","cv_S6","cv_S7","cv_S8","cv_S9"),"Kupffer.T","./CellPhoneDB/Vein/","CV") 
Score_matrix_generation_interacting_pair_oi_per_mouse(
  "./CellPhoneDB/Vein/",
  c("pv_S1","pv_S2","pv_S3","pv_S4","pv_S5","pv_S6","pv_S7","pv_S8","pv_S9"),"Kupffer.T","./CellPhoneDB/Vein/","PV") 

###DE analysis between CV and PV of interacting cells of interest 
DE_CellPhoneDB("./CellPhoneDB/Vein/CV_per_mouse_Kupffer.T.csv","./CellPhoneDB/Vein/PV_per_mouse_Kupffer.T.csv",
               c("cv_S1","cv_S2","cv_S3","cv_S4","cv_S5","cv_S6","cv_S7","cv_S8","cv_S9"),
               c("pv_S1","pv_S2","pv_S3","pv_S4","pv_S5","pv_S6","pv_S7","pv_S8","pv_S9"),
               "CV","PV","./CellPhoneDB/Vein/","Kupffer.T","vein")

########## Plotting of DE analysis in Barplot with CV on the left side and PV on the right side ##########
#include Pvalue FDR from DE analysis
top <- read.csv("./CellPhoneDB/Vein/Kupffer.T_vein_comp_LM.csv")
top <- top[, colnames(top) %in% c("X","FDR")]
colnames(top) <- c("interacting_pair","FDR")

###compare sig scores of CV and PV 
file1 <- "./CellPhoneDB/Vein/cv/significant_means.txt"
file2 <- "./CellPhoneDB/Vein/pv/significant_means.txt"

L_R_CellPhoneDB_comp_2samples("Kupffer.T",file1,file2,"./CellPhoneDB/Vein/", "sig_means_vein")

df <- read.csv("./CellPhoneDB/Vein/interactions_comp_two_cond_Kupffer.T_sig_means_vein.csv")

df$X <- NULL
colnames(df) <- c("interacting_pair","CV","PV")

df_cv <- df[,c("interacting_pair","CV")]
df_cv$vein <- "CV"
colnames(df_cv) <- c("interacting_pair","interaction_score","condition")
#make minus to put on one axis, has to be changed to + scores later in Illustrator 
df_cv$interaction_score <- -(df_cv$interaction_score)
df_cv2 <- merge(df_cv,top)

df_pv <- df[,c("interacting_pair","PV")]
df_pv$vein <- "PV"
colnames(df_pv) <- c("interacting_pair","interaction_score","condition")
df_pv2 <- merge(df_pv,top)

df1 <- rbind(df_cv2,df_pv2)

##make barplot 
#plot until p val 0.7
df2 <- df1[df1$FDR <=0.7,]

p <- ggplot(df2, aes(x=interaction_score, y=reorder(interacting_pair,-FDR),fill = FDR)) + theme_classic() +
  geom_bar(stat = "identity",position = "identity") + 
  scale_fill_gradientn(colours = c("darkblue","grey","darkred"),limits = c(0,1),breaks = c(0,0.1,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)) + 
  theme(axis.text = element_text(size = 15))  + theme(axis.text.x = element_text(angle = 90)) +
  theme(axis.title= element_text(size = 25)) + 
  ggtitle("L-R interactions - CellPhoneDB (Kupffer|T)") + xlab("Significant interaction score") + 
  ylab("L-R interaction")  +
  theme(plot.title = element_text(size = 15, face = "bold")) + theme(legend.text = element_text(size = 22),
                                                                     legend.title= element_text(size = 22)) + 
  guides(fill=guide_legend(title="P-value of enrichment")) 
p + ggsave("./figures/2.5.3/Interaction_KC_T_vein_0.5.pdf", width = 12, height = 4)
p + ggsave("./figures/2.5.3/Interaction_KC_T_vein_0.5.svg", width = 12, height = 4)








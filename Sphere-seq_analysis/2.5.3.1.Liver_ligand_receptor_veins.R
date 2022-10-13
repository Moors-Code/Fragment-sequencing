########## Part 2.5.3.1: Ligand-receptor interaction analysis between central and portal veins of injected samples ##########
#This part first generates input files for CellPhoneDB for CV and PV vein areas 
#Significant interactions are compared between vein areas to find specific ones 

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

#subset veins 
Idents(Injected) <- "vein"
cv <- subset(Injected, idents = "CV")
pv <- subset(Injected, idents = "PV")

#input file generation 
Input_files_CellPhoneDB_generation(cv,'annotation.broad',"cv","./CellPhoneDB/Vein/") 
Input_files_CellPhoneDB_generation(pv,'annotation.broad',"pv","./CellPhoneDB/Vein/") 

########## Run CellPhoneDB_veins.sh ##########
#/mnt/khandler/R_projects/Sphere-sequencing/Sphere-seq_analysis/2.5.3.2.CellPHoneDB_veins.sh
file1_mean <- "./CellPhoneDB/Vein/cv/significant_means.txt"
file1_pval <- "./CellPhoneDB/Vein/cv/pvalues.txt"

file2_mean <- "./CellPhoneDB/Vein/pv/significant_means.txt"
file2_pval <- "./CellPhoneDB/Vein/pv/pvalues.txt"

###Kupffer T cell interactions CV compared to PV 
#combine Kupffer and T cell interactions from CV and PV analysis 
L_R_CellPhoneDB_comp_2samples("Kupffer.T",file1_mean,file1_pval,file2_mean,file2_pval,"./figures/2.5.3/","vein")

##Plotting for L-R pairs 
df <- read.csv("./figures/2.5.3/interactions_comp_two_cond_Kupffer.T_vein.csv")

#calculate difference between distal and proximal to order interacting pairs with that 
#make values positive 
df$diff <- df$sample1_mean - df$sample2_mean
df$diff <- abs(df$diff)

##split the plot so you can put cv and pv scores in one column for interaction_score and one for pvals
#df <- df[df$sample1_mean == 0 | df$sample2_mean == 0,]

df_cv <- df[,c("interacting_pair", "sample1_mean","sample1_pval","diff")]
df_cv$vein <- "CV"
colnames(df_cv) <- c("interacting_pair","interaction_score","pvalue","diff", "condition")
#make minus to put on one axis, has to be changed to + scores later in Illustrator 
df_cv$interaction_score <- -(df_cv$interaction_score)

df_pv <- df[,c("interacting_pair", "sample2_mean","sample2_pval","diff")]
df_pv$vein <- "PV"
colnames(df_pv) <- c("interacting_pair","interaction_score","pvalue", "diff", "condition")

df1 <- rbind(df_cv,df_pv)

df1 <- df1[df1$pval < 0.05,]

#make barplot 
p <- ggplot(df1, aes(x=interaction_score, y=reorder(interacting_pair,+interaction_score),fill = pvalue)) + theme_classic() +
  geom_bar(stat = "identity",position = "identity") + 
  scale_fill_gradientn(colours = colorRampPalette(rev(brewer.pal(n = 15, name ="RdBu")))(256)) + 
  theme(axis.text = element_text(size = 15))  + theme(axis.text.x = element_text(angle = 90)) +
  theme(axis.title= element_text(size = 25)) + 
  ggtitle("L-R interactions - CellPhoneDB (Kupffer|T)") + xlab("Interaction score") + 
  ylab("L-R interaction")  +
  theme(plot.title = element_text(size = 15, face = "bold")) + theme(legend.text = element_text(size = 22),
                                                                     legend.title= element_text(size = 22)) + 
  guides(fill=guide_legend(title="P-value of enrichment")) 
p + ggsave("./figures/2.5.3/Interaction_KC_T_vein.pdf", width = 12, height = 10)
p + ggsave("./figures/2.5.3/Interaction_KC_T_vein.svg", width = 12, height = 10)


###Kupffer B cell interactions CV compared to PV 
#combine Kupffer and B cell interactions from CV and PV analysis 
L_R_CellPhoneDB_comp_2samples("Kupffer.B",file1_mean,file1_pval,file2_mean,file2_pval,"./figures/2.5.3/","vein")

##Plotting for L-R pairs 
df <- read.csv("./figures/2.5.3/interactions_comp_two_cond_Kupffer.B_vein.csv")

#calculate difference between distal and proximal to order interacting pairs with that 
#make values positive 
df$diff <- df$sample1_mean - df$sample2_mean
df$diff <- abs(df$diff)

##split the plot so you can put cv and pv scores in one column for interaction_score and one for pvals
#df <- df[df$sample1_mean == 0 | df$sample2_mean == 0,]

df_cv <- df[,c("interacting_pair", "sample1_mean","sample1_pval","diff")]
df_cv$vein <- "CV"
colnames(df_cv) <- c("interacting_pair","interaction_score","pvalue","diff", "condition")
#make minus to put on one axis, has to be changed to + scores later in Illustrator 
df_cv$interaction_score <- -(df_cv$interaction_score)

df_pv <- df[,c("interacting_pair", "sample2_mean","sample2_pval","diff")]
df_pv$vein <- "PV"
colnames(df_pv) <- c("interacting_pair","interaction_score","pvalue", "diff", "condition")

df1 <- rbind(df_cv,df_pv)

df1 <- df1[df1$pval < 0.05,]

#make barplot 
p <- ggplot(df1, aes(x=interaction_score, y=reorder(interacting_pair,+interaction_score))) + theme_classic() +
  geom_bar(stat = "identity",position = "identity", fill = "#1B3360") + 
  scale_fill_gradientn(colours = colorRampPalette(rev(brewer.pal(n = 15, name ="RdBu")))(256)) + 
  theme(axis.text = element_text(size = 15))  + theme(axis.text.x = element_text(angle = 90)) +
  theme(axis.title= element_text(size = 25)) + 
  ggtitle("L-R interactions - CellPhoneDB (Kupffer|B)") + xlab("Interaction score") + 
  ylab("L-R interaction")  +
  theme(plot.title = element_text(size = 15, face = "bold")) + theme(legend.text = element_text(size = 22),
                                                                     legend.title= element_text(size = 22)) + 
  guides(fill=guide_legend(title="P-value of enrichment")) 
p + ggsave("./figures/2.5.3/Interaction_KC_B_vein.pdf", width = 12, height = 10)
p + ggsave("./figures/2.5.3/Interaction_KC_B_vein.svg", width = 12, height = 10)






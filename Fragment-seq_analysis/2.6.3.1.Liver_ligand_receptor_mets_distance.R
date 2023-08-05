########## Part 2.6.3.1: Ligand-receptor interaction analysis between proximal and distal areas from metastatic sites ##########
#This part first generates input files for CellPhoneDB for proximal and distal areas
#Then significant interactions are compared between areas to find specific ones 

########## Prepare environment ##########
###Setting the working directory 
setwd("/mnt/khandler/R_projects/Fragment-sequencing/Fragment-seq_analysis/")

###Load packages and functions 
source("./functions_and_packages/1.Packages.R")
source("./functions_and_packages/7.Functions_ligand_receptor_analysis.R")

###load R object
metastasis <- readRDS("./data_files_generated/LiverMerged_afterBC_anno_BS_5cells_zC_lobules_mets_distance.Rda")

########## CellPhoneDB input file generation per Mets distance ##########
#files are then used in code 2.6.4.2.CellPhoneDB_mets_distance.sh

#load human and mouse ensemble symbols
human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/") 
mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")

#subset Mets distances 
Idents(metastasis) <- "Mets_distance"
proximal <- subset(metastasis, idents = "proximal")
distal <- subset(metastasis, idents = "distal")

#input file generation 
Input_files_CellPhoneDB_generation(proximal,'annotation.broad',"proximal","./CellPhoneDB/Mets_distance/") 
Input_files_CellPhoneDB_generation(distal,'annotation.broad',"distal","./CellPhoneDB/Mets_distance/") 

########## Run CellPhoneDB_mets_distance.sh ##########
#/mnt/khandler/R_projects/Fragment-sequencing/Fragment-seq_analysis/2.6.3.2.CellPHoneDB_mets_distance.sh

########## Comparison of L-R interactions between proximal and distal areas ##########

###read in ouput files
file1_mean <- "./CellPhoneDB/Mets_distance/distal/significant_means.txt"
file1_pval <- "./CellPhoneDB/Mets_distance/distal/pvalues.txt"

file2_mean <- "./CellPhoneDB/Mets_distance/proximal/significant_means.txt"
file2_pval <- "./CellPhoneDB/Mets_distance/proximal/pvalues.txt"

###Monocytes T cell interactions distal compared to proximal 
#combine Monocytes and T cell interactions from distal and proximal analysis 
L_R_CellPhoneDB_comp_2samples("Monocytes.T",file1_mean,file1_pval,file2_mean,file2_pval,"./figures/2.6.3/","Mets_distance")

##Plotting 
df <- read.csv("./figures/2.6.3/interactions_comp_two_cond_Monocytes.T_Mets_distance.csv")

#calculate difference between distal and proximal to order interacting pairs in decrasing order
#make values positive 
df$diff <- df$sample1_mean - df$sample2_mean
df$diff <- abs(df$diff)

##split the plot so you can put distal and proximal scores in one column for interaction_score and one for pvals
df_dist <- df[,c("interacting_pair", "sample1_mean","sample1_pval","diff")]
df_dist$Mets_distance <- "distal"
colnames(df_dist) <- c("interacting_pair","interaction_score","pvalue","diff", "condition")
#make minus to put on one axis, has to be changed to + scores later in Illustrator 
df_dist$interaction_score <- -(df_dist$interaction_score)

df_prox <- df[,c("interacting_pair", "sample2_mean","sample2_pval","diff")]
df_prox$Mets_distance <- "proximal"
colnames(df_prox) <- c("interacting_pair","interaction_score","pvalue","diff", "condition")

df1 <- rbind(df_dist,df_prox)

df1 <- df1[df1$pval < 0.05,]


#make barplot 
p <- ggplot(df1, aes(x=interaction_score, y=reorder(interacting_pair,+diff),fill = pvalue)) + theme_classic() +
  geom_bar(stat = "identity",position = "identity") + 
  scale_fill_gradientn(colours = colorRampPalette(rev(brewer.pal(n = 15, name ="RdBu")))(256)) + 
  theme(axis.text = element_text(size = 15))  + theme(axis.text.x = element_text(angle = 90)) +
  theme(axis.title= element_text(size = 25)) + 
  ggtitle("L-R interactions - CellPhoneDB (Mono|T)") + xlab("Interaction score") + 
  ylab("L-R interaction")  +
  theme(plot.title = element_text(size = 15, face = "bold")) + theme(legend.text = element_text(size = 22),
                                                                     legend.title= element_text(size = 22)) + 
  guides(fill=guide_legend(title="P-value of enrichment")) 
p + ggsave("./figures/2.6.3/Interaction_Mono_T_vein_0.5.pdf", width = 12, height = 10)
p + ggsave("./figures/2.6.3/Interaction_Mono_T_vein_0.5.svg", width = 12, height = 10)





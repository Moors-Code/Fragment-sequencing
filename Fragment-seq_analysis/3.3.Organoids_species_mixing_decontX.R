########## Part 3.3: Organoids species mixing decontX ##########
#This part applies dexontX to remove cell free RNA because when plotting mouse and human UMI counts per cell there was a lot of cells showing mouse and human UMI counts 
#these we need to remove because they would negatively impact the mouse and human UMI counts per fragment 

########## Prepare environment ##########
###Setting the working directory 
setwd("/mnt/khandler/R_projects/Fragment-sequencing/Fragment-seq_analysis/")

###Load packages and functions 
source("./functions_and_packages/1.Packages.R")

###Load data 
org_mix_ensemble <- readRDS("./data_files_generated/SpS_Organoid_mixing_ensembleIDs_BS_mito0.3.Rda")

########## DecontX ##########
#dexontX needs a single cell Experiment object 
sce <- as.SingleCellExperiment(org_mix_ensemble)
sce.delta <- decontX(sce)

#check decontamination 
plotDecontXContamination(sce.delta)

#decontaminated count matrix 
head(decontXcounts(sce.delta))

########## plot mouse and human UMI counts per cell after decontamination ##########
org_mix_ensemble_AE <- decontXcounts(sce.delta)
org_mix_ensemble_AE_df <- as.data.frame(org_mix_ensemble_AE)

#split data frame in human and mouse 
human.features <- grep(pattern = "^GRCh38-",x = rownames(x = org_mix_ensemble_AE_df), value = TRUE)
org_mix_AE_df_human <- org_mix_ensemble_AE_df[which(row.names(org_mix_ensemble_AE_df) %in% human.features),]

mouse.features <- grep(pattern = "^mm10-",x = rownames(x = org_mix_ensemble_AE_df), value = TRUE)
org_mix_AE_df_mouse <- org_mix_ensemble_AE_df[which(row.names(org_mix_ensemble_AE_df) %in% mouse.features),]

#calculate sum for each dataframe 
#first transpose and then calculate rowSums 
#human
org_mix_AE_df_human_t <- t(org_mix_AE_df_human)
org_mix_AE_df_human_t_df <- as.data.frame(org_mix_AE_df_human_t)
org_mix_AE_df_human_t_df["Total"] <- rowSums(org_mix_AE_df_human_t_df)
org_mix_AE_df_human_t_df_total <- org_mix_AE_df_human_t_df[,"Total"]
org_mix_AE_df_human_t_df_total_df <- as.data.frame(org_mix_AE_df_human_t_df_total)
rownames(org_mix_AE_df_human_t_df_total_df) <- rownames(org_mix_AE_df_human_t_df)

#mouse
org_mix_AE_df_mouse_t <- t(org_mix_AE_df_mouse)
org_mix_AE_df_mouse_t_df <- as.data.frame(org_mix_AE_df_mouse_t)
org_mix_AE_df_mouse_t_df["Total"] <- rowSums(org_mix_AE_df_mouse_t_df)
org_mix_AE_df_mouse_t_df_total <- org_mix_AE_df_mouse_t_df[,"Total"]
org_mix_AE_df_mouse_t_df_total_df <- as.data.frame(org_mix_AE_df_mouse_t_df_total)
rownames(org_mix_AE_df_mouse_t_df_total_df) <- rownames(org_mix_AE_df_mouse_t_df)

#merge both dataframes by rownames
merged_df <- cbind(org_mix_AE_df_human_t_df_total_df,org_mix_AE_df_mouse_t_df_total_df)
colnames(merged_df) <- c("human","mouse")

#scatterplot  --> no contaminated cells left 
p <- ggplot(merged_df, aes(x=human, y=mouse)) + theme_classic() +geom_point() +
  ggtitle("Mouse and human UMI counts per cell after DecontX")+ xlab("Human transcripts") + ylab("Mouse transcripts") +
  theme(title = element_text(size = 25))+
  theme(axis.title= element_text(size = 25)) +  theme(axis.text = element_text(size = 30)) 
p + ggsave("./figures/3.3/Orgnoid_Mixing_species_singleCells_afterDecontX.pdf",width = 12, height = 10)
p + ggsave("./figures/3.3/Orgnoid_Mixing_species_singleCells_afterDecontX.svg",width = 12, height = 10)

########## save R object ##########
#convert to Seurat object and save 
sce.delta.seurat <- convertSCEToSeurat(sce.delta, copyDecontX = TRUE,copyColData = TRUE)
saveRDS(sce.delta.seurat,"./data_files_generated/SpS_Organoid_mixing_ensembleIDs_BS_mito0.3_decontX.Rda")


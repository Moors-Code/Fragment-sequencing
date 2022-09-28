########## Part 3.2: Organoids species mixing quality control ##########
#This part applies QC measures and normalization to CRC Organoid mixing species experiment 

########## Prepare environment ##########
###Setting the working directory 
setwd("/mnt/khandler/R_projects/Sphere-sequencing/Sphere-seq_analysis/")

###Load packages and functions 
source("./functions_and_packages/1.Packages.R")

###Load data 
org_mix <- readRDS("./data_files_generated/SpS_Organoid_mixing_BS.Rda")
org_mix_ensemble <- readRDS("./data_files_generated/SpS_Organoid_mixing_ensembleIDs_BS.Rda")

########## remove Doublets and Negatives ##########
spheres <- as.data.frame(table(org_mix$sphere))$Var1
spheres <- as.character(spheres)
spheres <- spheres[!spheres %in% c("Negative","Doublet")]
Idents(org_mix) <- "sphere"
org_mix <- subset(org_mix, idents = spheres)

spheres <- as.data.frame(table(org_mix_ensemble$sphere))$Var1
spheres <- as.character(spheres)
spheres <- spheres[!spheres %in% c("Negative","Doublet")]
Idents(org_mix_ensemble) <- "sphere"
org_mix_ensemble <- subset(org_mix_ensemble, idents = spheres)

########## apply mitochondrial QC cutoff to remove low quality cells ##########
mito.features <- grep(pattern = "^MT-|^mt-", x = rownames(x = org_mix), value = TRUE)
percent.mito <- Matrix::colSums(x = GetAssayData(object = org_mix, slot = 'counts')[mito.features, ]) / Matrix::colSums(x = GetAssayData(object = org_mix, slot = 'counts'))
org_mix <- AddMetaData(object = org_mix, metadata = percent.mito, col.name = "percent.mito")
Idents(org_mix) <- "orig.ident"
VlnPlot(org_mix, feature = "percent.mito")
org_mix <- subset(org_mix, percent.mito < 0.3 )

########## check mouse and human UMI counts per cell to assess general quality ##########
#only use cells that passed mitochondrial cutoff 
#extract cell IDs that passed mito cutoff from org_mix_ensemble object to then plot human and mouse UMI counts 
org_mix_ensemble$cell_bc <- rownames(org_mix_ensemble@meta.data)
good_bc <- rownames(org_mix@meta.data)
Idents(org_mix_ensemble) <- "cell_bc"
org_mix_ensemble <- subset(org_mix_ensemble, idents = good_bc)

org_mix_ensemble_AE <-org_mix_ensemble@assays$RNA@counts
org_mix_ensemble_AE_df <- as.data.frame(org_mix_ensemble_AE)

#split data frame in human and mouse 
human.features <- grep(pattern = "^GRCh38-",x = rownames(x = org_mix_ensemble_AE_df), value = TRUE)
org_mix_ensemble_AE_df_human <- org_mix_ensemble_AE_df[which(row.names(org_mix_ensemble_AE_df) %in% human.features),]

mouse.features <- grep(pattern = "^mm10-",x = rownames(x = org_mix_ensemble_AE_df), value = TRUE)
org_mix_ensemble_AE_df_mouse <- org_mix_ensemble_AE_df[which(row.names(org_mix_ensemble_AE_df) %in% mouse.features),]

#calculate sum for each dataframe 
#first transpose and then calculate rowSums 
#human
org_mix_ensemble_AE_df_human_t <- t(org_mix_ensemble_AE_df_human)
org_mix_ensemble_AE_df_human_t_df <- as.data.frame(org_mix_ensemble_AE_df_human_t)
org_mix_ensemble_AE_df_human_t_df["Total"] <- rowSums(org_mix_ensemble_AE_df_human_t_df)
org_mix_ensemble_AE_df_human_t_df_total <- org_mix_ensemble_AE_df_human_t_df[,"Total"]
org_mix_ensemble_AE_df_human_t_df_total_df <- as.data.frame(org_mix_ensemble_AE_df_human_t_df_total)
rownames(org_mix_ensemble_AE_df_human_t_df_total_df) <- rownames(org_mix_ensemble_AE_df_human_t_df)

#mouse
org_mix_ensemble_AE_df_mouse_t <- t(org_mix_ensemble_AE_df_mouse)
org_mix_ensemble_AE_df_mouse_t_df <- as.data.frame(org_mix_ensemble_AE_df_mouse_t)
org_mix_ensemble_AE_df_mouse_t_df["Total"] <- rowSums(org_mix_ensemble_AE_df_mouse_t_df)
org_mix_ensemble_AE_df_mouse_t_df_total <- org_mix_ensemble_AE_df_mouse_t_df[,"Total"]
org_mix_ensemble_AE_df_mouse_t_df_total_df <- as.data.frame(org_mix_ensemble_AE_df_mouse_t_df_total)
rownames(org_mix_ensemble_AE_df_mouse_t_df_total_df) <- rownames(org_mix_ensemble_AE_df_mouse_t_df)

#merge both dataframes by rownames
merged_df <- cbind(org_mix_ensemble_AE_df_human_t_df_total_df,org_mix_ensemble_AE_df_mouse_t_df_total_df)
colnames(merged_df) <- c("human","mouse")

#scatterplot --> there is still quite some cell showing human and mouse UMI counts, therefore we apply dexontX to remove cell free RNA in 3.3.
p <- ggplot(merged_df, aes(x=human, y=mouse)) + theme_classic() + geom_point() +
  ggtitle("Mouse and human UMI counts per cell")+ xlab("Human transcripts") + ylab("Mouse transcripts") +
  theme(title = element_text(size = 25))+
  theme(axis.title= element_text(size = 25)) +  theme(axis.text = element_text(size = 30)) 
p + ggsave("./figures/3.2/Orgnoid_Mixing_species_singleCells_preDecontX.pdf",width = 12, height = 10)
p + ggsave("./figures/3.2/Orgnoid_Mixing_species_singleCells_preDecontX.svg",width = 12, height = 10)

########## save R objects ##########
saveRDS(org_mix, file = "./data_files_generated/SpS_Organoid_mixing_BS_mito0.3.Rda")
saveRDS(org_mix_ensemble, file = "./data_files_generated/SpS_Organoid_mixing_ensembleIDs_BS_mito0.3.Rda")




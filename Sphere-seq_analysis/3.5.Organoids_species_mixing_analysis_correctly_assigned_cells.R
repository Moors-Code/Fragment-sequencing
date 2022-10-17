########## Part 3.5: Organoids species mixing - analysis correctly assigned cells ##########
#This part investigates the fraction of correctly assigned cells by matching the cell type annotation with well information of species 

########## Prepare environment ##########
###Setting the working directory 
setwd("/mnt/khandler/R_projects/Sphere-sequencing/Sphere-seq_analysis/")

###Load packages and functions 
source("./functions_and_packages/1.Packages.R")

###Load data 
org_mix_ensemble <- readRDS("./data_files_generated/SpS_Organoid_mixing_ensembleIDs_BS_mito0.3_decontX_annotated.Rda")

########## calculate correctly assigned cells ##########
###annotated human and mouse wells based on sorting and based on prior knowledge that 19,36 and 106 contain 100% mouse cells 
#even though they were GFP negative (sorting error)
human_wells <- paste("Bar",c(1:18,20:35,37:48,97:105,107:144,194:240), sep = "")
mouse_wells <- paste("Bar",c(49:96,145:193,241:288,19,36,106), sep = "")

#add meta data column with well information 
org_mix_ensemble$species_well <- NA
org_mix_ensemble@meta.data <- org_mix_ensemble@meta.data %>%
  mutate(species_well = case_when(
    sphere %in% human_wells ~ "human_well",
    sphere %in% mouse_wells ~ "mouse_well",
    TRUE ~ NA_character_))

#subset by well and check annotation identity
Idents(org_mix_ensemble) <- "species_well"
org_mix_ensemble_human <- subset(org_mix_ensemble,idents = "human_well")
org_mix_ensemble_mouse <- subset(org_mix_ensemble,idents = "mouse_well")

#extract wrongly assigned cells from human and mouse objects
wrong_mouse_in_human <- row.names(org_mix_ensemble_human@meta.data)[which(org_mix_ensemble_human$annotation %in% "mouse")]
wrong_human_in_mouse <- row.names(org_mix_ensemble_mouse@meta.data)[which(org_mix_ensemble_mouse$annotation %in% "human")]
wrong_total <- c(wrong_mouse_in_human,wrong_human_in_mouse)

print(length(org_mix_ensemble@active.ident)) #8272
print(length(wrong_total)) #410
#5% wrong 

########## Plot results in scatterplot ##########
#x and y are human and mouse UMI counts per sphere 
#in color plot the percentage of wrongly assigned cells 

###extract barcodes of human cells within mouse wells = wrong
Idents(org_mix_ensemble) <- "species_well"
mouse_wells <- subset(org_mix_ensemble, idents = "mouse_well")
Idents(mouse_wells) <- "annotation"
mouse_wells_human_cells <- subset(mouse_wells, idents = "human")

#get total number of spheres of wells with wrong cells 
wrong_bc_mouse_wells <- as.data.frame(table(mouse_wells_human_cells$sphere))$Var1
Idents(org_mix_ensemble) <- "sphere"
wrong_spheres_all_cell_number <- subset(org_mix_ensemble, idents = wrong_bc_mouse_wells)
table(wrong_spheres_all_cell_number$sphere)

total_number_wrong_bc_mouse_wells <- as.data.frame(table(wrong_spheres_all_cell_number$sphere))
colnames(total_number_wrong_bc_mouse_wells) <- c("sphere","total")
wrong_number_wrong_bc_mouse_wells <- as.data.frame(table(mouse_wells_human_cells$sphere))
colnames(wrong_number_wrong_bc_mouse_wells) <- c("sphere","wrong")

wrong_mouse_wells_df <- cbind(total_number_wrong_bc_mouse_wells,wrong_number_wrong_bc_mouse_wells)
wrong_mouse_wells_df$percent.wrong <- (wrong_mouse_wells_df$wrong * 100 )/wrong_mouse_wells_df$total

###extract barcodes of mouse cells within human wells and calculate percentage = wrong
Idents(org_mix_ensemble) <- "species_well"
human_wells <- subset(org_mix_ensemble, idents = "human_well")
Idents(human_wells) <- "annotation"
human_wells_mouse_cells <- subset(human_wells, idents = "mouse")

#get total number of spheres of wells with wrong cells 
wrong_bc_human_wells <- as.data.frame(table(human_wells_mouse_cells$sphere))$Var1
Idents(org_mix_ensemble) <- "sphere"
wrong_spheres_all_cell_number <- subset(org_mix_ensemble, idents = wrong_bc_human_wells)
table(wrong_spheres_all_cell_number$sphere)

total_number_wrong_bc_human_wells <- as.data.frame(table(wrong_spheres_all_cell_number$sphere))
colnames(total_number_wrong_bc_human_wells) <- c("sphere","total")
wrong_number_wrong_bc_human_wells <- as.data.frame(table(human_wells_mouse_cells$sphere))
colnames(wrong_number_wrong_bc_human_wells) <- c("sphere","wrong")

wrong_human_wells_df <- cbind(total_number_wrong_bc_human_wells,wrong_number_wrong_bc_human_wells)
wrong_human_wells_df$percent.wrong <- (wrong_human_wells_df$wrong * 100 )/wrong_human_wells_df$total

###extract the rest of the sphere-bcs that are not within the wrong mouse or human wells, these are spheres with 0% wrong cells
all_bc <- as.data.frame(table(org_mix_ensemble$sphere))$Var1
rest_sphere_bc <- setdiff(all_bc,wrong_human_wells_df$sphere )
rest_sphere_bc <- setdiff(rest_sphere_bc,wrong_mouse_wells_df$sphere )

#add 0 to the rest barcodes, these are spheres with 0% wrong cells 
values <- rep(0,length(rest_sphere_bc))
rest_bc_df <- data.frame(rest_sphere_bc,values)
colnames(rest_bc_df) <- c("sphere","percent.wrong")

###make dataframe with all sphere barcodes and percentage of wrong cells per sphere
wrong_mouse_wells_df <- wrong_mouse_wells_df[,c(1,5)]
wrong_human_wells_df <- wrong_human_wells_df[,c(1,5)]

#combine
all_spheres_df <- rbind(rest_bc_df,wrong_mouse_wells_df,wrong_human_wells_df)

###make a data frame with counts 
org_mix_ensemble_AE <- AverageExpression(org_mix_ensemble,assays = "RNA", slot = "count", group.by = "sphere")
org_mix_ensemble_AE_df <- as.data.frame(org_mix_ensemble_AE)

human.features <- grep(pattern = "^GRCh38-",x = rownames(x = org_mix_ensemble_AE_df), value = TRUE)
org_mix_ensemble_AE_df_human <- org_mix_ensemble_AE_df[which(row.names(org_mix_ensemble_AE_df) %in% human.features),]

mouse.features <- grep(pattern = "^mm10-",x = rownames(x = org_mix_ensemble_AE_df), value = TRUE)
org_mix_ensemble_AE_df_mouse <- org_mix_ensemble_AE_df[which(row.names(org_mix_ensemble_AE_df) %in% mouse.features),]

###calculate sum for each dataframe 
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

###merge both dataframes by rownames
merged_df <- cbind(org_mix_ensemble_AE_df_human_t_df_total_df,org_mix_ensemble_AE_df_mouse_t_df_total_df)
colnames(merged_df) <- c("human","mouse")

merged_df$sphere <- str_replace(rownames(merged_df),pattern = "RNA.",replacement = "")

#match merged_df with all_spheres_df to get percent of wrong cells per sphere 
merged_df_percent <- merge(merged_df, all_spheres_df)

###scatterplot 
p <- ggplot(merged_df_percent, aes(x=human, y=mouse)) + theme_classic() +
  geom_point(size = 5,aes(colour = percent.wrong))  + 
  scale_color_gradient(low = "blue", high = "yellow") +
  ggtitle("Mouse and human UMI counts per sphere after DeconX")+ 
  theme(title = element_text(size = 25))  + 
  theme(axis.title= element_text(size = 25)) +  theme(axis.text = element_text(size = 30)) 
p + ggsave("./figures/3.5/Orgnoid_Mixing_species_per_sphere.pdf",width = 12, height = 10)
p + ggsave("./figures/3.5/Orgnoid_Mixing_species_per_sphere.svg",width = 12, height = 10)


########## Barplot of fraction of cells that are wrongly and correctly assigned ##########
all_cells_amount <- length(org_mix_ensemble@active.ident)
wrong <- length(wrong_total)
correct <- all_cells_amount - wrong
condition <- c("correct","wrong")
cell_number <- c(correct,wrong)

df <- data.frame(condition,cell_number)

#barplot 
p <- ggplot(df, aes( y=cell_number, x = condition)) + 
  geom_bar( stat="identity") + 
  xlab("Condition") + ylab("Cell number") +ggtitle("Wrongly and correctly assigned cells") + 
  theme(axis.text = element_text(size = 30)) + 
  theme(axis.title= element_text(size = 25)) + 
  theme(plot.title = element_text(size = 25, face = "bold")) + 
  theme(legend.title = element_text(size = 30), legend.text = element_text(size = 30)) 
p + ggsave("./figures/3.5/barplot_wrong_and_correct_cells.pdf", width = 8, height = 8)
p + ggsave("./figures/3.5/barplot_wrong_and_correct_cells.svg", width = 8, height = 8)



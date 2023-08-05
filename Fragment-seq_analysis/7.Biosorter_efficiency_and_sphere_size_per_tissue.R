########## Part 7: Biosorter efficiency and sphere size per tissue ##########
#This part produces plots to show sorting efficiency after plate sorting, imaging and counting of wells with one sphere, 0 and multiples (2) 
#In this part also a plot is generated that compares sphere size distributions between differen tissues 

########## Prepare environment ##########
###Setting the working directory 
setwd("/mnt/khandler/R_projects/Sphere-sequencing/Sphere-seq_analysis/")

###Load packages and functions 
source("./functions_and_packages/1.Packages.R")

########## counts of wells with 0 spheres and wells with 1 sphere ##########
#plate 1
wells_0_spheres <- 7
wells_1_spheres <- 89
well_2_spheres <- 0
counts <- c(wells_0_spheres, wells_1_spheres, well_2_spheres)
df1 <- as.data.frame(counts)
rownames(df1) <- c("0_spheres","1_spheres","2_spheres")
df1$plate <- "plate1"
df1$condition <- c("0_wells","1_wells","2_wells")
#plate 2
wells_0_spheres <- 10
wells_1_spheres <- 85
well_2_spheres <- 1
counts <- c(wells_0_spheres, wells_1_spheres, well_2_spheres)
df2 <- as.data.frame(counts)
rownames(df2) <- c("0_spheres","1_spheres","2_spheres")
df2$plate <- "plate2"
df2$condition <- c("0_wells","1_wells","2_wells")
#plate 3
wells_0_spheres <- 5
wells_1_spheres <- 91
well_2_spheres <- 0
counts <- c(wells_0_spheres, wells_1_spheres, well_2_spheres)
df3 <- as.data.frame(counts)
rownames(df3) <- c("0_spheres","1_spheres","2_spheres")
df3$plate <- "plate3"
df3$condition <- c("0_wells","1_wells","2_wells")
#plate 4
wells_0_spheres <- 4
wells_1_spheres <- 92
well_2_spheres <- 0
counts <- c(wells_0_spheres, wells_1_spheres, well_2_spheres)
df4 <- as.data.frame(counts)
rownames(df4) <- c("0_spheres","1_spheres","2_spheres")
df4$plate <- "plate4"
df4$condition <- c("0_wells","1_wells","2_wells")
#plate 5
wells_0_spheres <- 5
wells_1_spheres <- 91
well_2_spheres <- 0
counts <- c(wells_0_spheres, wells_1_spheres, well_2_spheres)
df5 <- as.data.frame(counts)
rownames(df5) <- c("0_spheres","1_spheres","2_spheres")
df5$plate <- "plate5"
df5$condition <- c("0_wells","1_wells","2_wells")
#plate 6
wells_0_spheres <- 11
wells_1_spheres <- 85
well_2_spheres <- 0
counts <- c(wells_0_spheres, wells_1_spheres, well_2_spheres)
df6 <- as.data.frame(counts)
rownames(df6) <- c("0_spheres","1_spheres","2_spheres")
df6$plate <- "plate6"
df6$condition <- c("0_wells","1_wells","2_wells")
#plate 7
wells_0_spheres <- 6
wells_1_spheres <- 90
well_2_spheres <- 0
counts <- c(wells_0_spheres, wells_1_spheres, well_2_spheres)
df7 <- as.data.frame(counts)
rownames(df7) <- c("0_spheres","1_spheres","2_spheres")
df7$plate <- "plate7"
df7$condition <- c("0_wells","1_wells","2_wells")
#plate 8
wells_0_spheres <- 7
wells_1_spheres <- 88
well_2_spheres <- 1
counts <- c(wells_0_spheres, wells_1_spheres, well_2_spheres)
df8 <- as.data.frame(counts)
rownames(df8) <- c("0_spheres","1_spheres","2_spheres")
df8$plate <- "plate8"
df8$condition <- c("0_wells","1_wells","2_wells")
#plate 9
wells_0_spheres <- 10
wells_1_spheres <- 86
well_2_spheres <- 0
counts <- c(wells_0_spheres, wells_1_spheres, well_2_spheres)
df9 <- as.data.frame(counts)
rownames(df9) <- c("0_spheres","1_spheres","2_spheres")
df9$plate <- "plate9"
df9$condition <- c("0_wells","1_wells","2_wells")

df <- rbind(df1, df2)
df <- rbind(df, df3)
df <- rbind(df, df4)
df <- rbind(df, df5)
df <- rbind(df, df6)
df <- rbind(df, df7)
df <- rbind(df, df8)
df <- rbind(df, df9)

p <- ggplot(df, aes(x = condition, y = counts, fill = plate, colour = plate))  +theme_classic() +
  geom_bar(stat = "identity", position = "dodge") + scale_fill_brewer(palette="Set1") +
  theme(axis.text = element_text(size = 30))   + 
  ggtitle("Sorting efficiency per plate") + xlab("Sphere count per well") + 
  ylab("Well Counts") + theme(axis.title= element_text(size = 25)) + 
  theme(plot.title = element_text(size = 25, face = "bold"))  + 
  theme(legend.title = element_text(size = 30), legend.text = element_text(size = 30)) +
  theme(legend.title = element_text(size = 30), legend.text = element_text(size = 30))
p + ggsave("./figures/7/Sort_efficiency_per_plate.pdf",width = 12, height = 10)  
p + ggsave("./figures/7/Sort_efficiency_per_plate.svg",width = 12, height = 10)  


########## plot sphere sizes from different biological systems in boxplot ##########
########## Linear model of standard sized beads ##########
#TOF from standard sized beads (60µm, 125µm and 175µm) was acquired and a linear model was generated 

###load acquired data of standard sized beads 
standard60 <- read.table("/mnt/khandler/R_projects/Sphere-sequencing/Sphere-seq_analysis/data/20200824_60beads_for_standard_curve.txt", sep = "\t", 
                         header = TRUE, blank.lines.skip = FALSE, stringsAsFactors = FALSE)
#remove header 
standard60 <- standard60[1:(which(standard60$Id == "")[1] - 1),]
standard125 <- read.table("/mnt/khandler/R_projects/Sphere-sequencing/Sphere-seq_analysis/data/20200824_125beads_for_standard_curve.txt", sep = "\t", 
                          header = TRUE, blank.lines.skip = FALSE, stringsAsFactors = FALSE)
standard125 <- standard125[1:(which(standard125$Id == "")[1] - 1),]
standard175 <- read.table("/mnt/khandler/R_projects/Sphere-sequencing/Sphere-seq_analysis/data/20200824_175beads_for_standard_curve.txt", sep = "\t", 
                          header = TRUE, blank.lines.skip = FALSE, stringsAsFactors = FALSE)
standard175 <- standard175[1:(which(standard175$Id == "")[1] - 1),]

#concatenate individual acquisitions
bead_df <- rbind(standard60, standard125, standard175)
bead_df$bead_size <- c(rep(60, nrow(standard60)), rep(125, nrow(standard125)), rep(175, nrow(standard175)))

###fit a linear model with TOF 
lm.bead_tof <- lm(bead_size ~ TOF, data = bead_df)

##### Spleen samples 
### Read data from biosorter outputs 
## Sample 1 
#read in acquired data 
S1_Plate1 <- read.table("/mnt/khandler/RStudio_Data/Biosorter/Spleen_BS/plate1_25nM.txt", sep = "\t", 
                        header = TRUE, blank.lines.skip = FALSE, stringsAsFactors = FALSE)
S1_Plate1 <- S1_Plate1[1:(which(S1_Plate1$Id == "")[1] - 1),]
S1_Plate2 <- read.table("/mnt/khandler/RStudio_Data/Biosorter/Spleen_BS/plate1_25nM2.txt", sep = "\t", 
                        header = TRUE, blank.lines.skip = FALSE, stringsAsFactors = FALSE)
S1_Plate2 <- S1_Plate2[1:(which(S1_Plate2$Id == "")[1] - 1),]
S1_Plate3 <- read.table("/mnt/khandler/RStudio_Data/Biosorter/Spleen_BS/plate1_25nM3.txt", sep = "\t", 
                        header = TRUE, blank.lines.skip = FALSE, stringsAsFactors = FALSE)
S1_Plate3 <- S1_Plate3[1:(which(S1_Plate3$Id == "")[1] - 1),]
S1_Plate4 <- read.table("/mnt/khandler/RStudio_Data/Biosorter/Spleen_BS/plate1_25nM4.txt", sep = "\t", 
                        header = TRUE, blank.lines.skip = FALSE, stringsAsFactors = FALSE)
S1_Plate4 <- S1_Plate4[1:(which(S1_Plate4$Id == "")[1] - 1),]

#add plate ID 
S1_Plate1$plate_id <- paste0("plate_", c(rep(1, nrow(S1_Plate1))))
S1_Plate2$plate_id <- paste0("plate_", c(rep(1, nrow(S1_Plate2))))
S1_Plate3$plate_id <- paste0("plate_", c(rep(1, nrow(S1_Plate3))))
S1_Plate4$plate_id <- paste0("plate_", c(rep(1, nrow(S1_Plate4))))

#merge files
S1 <- full_join(S1_Plate1, S1_Plate2)
S1 <- full_join(S1,S1_Plate3)
S1 <- full_join(S1,S1_Plate4)

###Sample 2 
S2_Plate1 <- read.table("/mnt/khandler/RStudio_Data/Biosorter/Spleen_BS/plate1_50nM.txt", sep = "\t", 
                        header = TRUE, blank.lines.skip = FALSE, stringsAsFactors = FALSE)
S2_Plate1 <- S2_Plate1[1:(which(S2_Plate1$Id == "")[1] - 1),]
S2_Plate2 <- read.table("/mnt/khandler/RStudio_Data/Biosorter/Spleen_BS/plate1_50nM2.txt", sep = "\t", 
                        header = TRUE, blank.lines.skip = FALSE, stringsAsFactors = FALSE)
S2_Plate2 <- S2_Plate2[1:(which(S2_Plate2$Id == "")[1] - 1),]
S2_Plate3 <- read.table("/mnt/khandler/RStudio_Data/Biosorter/Spleen_BS/plate1_50nM3.txt", sep = "\t", 
                        header = TRUE, blank.lines.skip = FALSE, stringsAsFactors = FALSE)
S2_Plate3 <- S2_Plate3[1:(which(S2_Plate3$Id == "")[1] - 1),]

S2_Plate1$plate_id <- paste0("plate_", c(rep(1, nrow(S2_Plate1))))
S2_Plate2$plate_id <- paste0("plate_", c(rep(1, nrow(S2_Plate2))))
S2_Plate3$plate_id <- paste0("plate_", c(rep(1, nrow(S2_Plate3))))

S2 <- full_join(S2_Plate1, S2_Plate2)
S2 <- full_join(S2,S2_Plate3)

### Sphere size prediction and integration into Seurat object 
##map output of Biosorter with sphere (MULTI-seq) ID 
S1_2 <- Biosorter_output_managing(S1,"_1")
S2_2 <- Biosorter_output_managing(S2,"_2")

#merge all samples 
biosorter_all <- full_join(S1_2, S2_2)

##Integrate sphere-size information with Seurat object 
Biosorter_data_seurat_integration(
  "./data_files_generated/Spleen_merged.Rda",
  biosorter_all,
  "./data_files_generated/Spleen_merged_BS.Rda")

##Read in Seurat objects 
spleen <- readRDS("./data_files_generated/Spleen_merged_BS.Rda")
liver <- readRDS("/mnt/khandler/R_projects/Sphere-sequencing/Sphere-seq_analysis/data_files_generated/LiverMerged_afterBC_anno_BS_5cells_zC_lobules.Rda")
organoids <- readRDS("./data_files_generated/SpS_Organoid_mixing_ensembleIDs_BS_mito0.3_decontX_annotated.Rda")

spleen_size <- BS_df_for_boxplot_per_sample(spleen,"sphere","sphere_size","spleen")
liver_size <- BS_df_for_boxplot_per_sample(liver,"sphere","sphere_size","liver")
organoids_size <- BS_df_for_boxplot_per_sample(organoids,"sphere","sphere_size","organoids")

#merge data frames
size_merged_df <- rbind(spleen_size,liver_size)
size_merged_df <- rbind(size_merged_df,organoids_size)

###Plot in boxplot 
p <- ggplot(size_merged_df,aes(x = sample,y = sphere_size, fill = sample)) +theme_classic() +
  geom_boxplot() +
  geom_jitter(position = position_jitter(seed = 1, width =0.4),size = 0.5) + 
  theme(axis.text = element_text(size = 30))  +
  ggtitle("Sphere size per sample and sphere") + xlab("Sample") + 
  ylab("Size (µm)") + theme(axis.title= element_text(size = 25)) + 
  theme(plot.title = element_text(size = 25, face = "bold")) + 
  theme(legend.title = element_text(size = 30), legend.text = element_text(size = 30)) 
p + ggsave("./figures/7/Size_boxplot_perTissue.pdf",width = 15, height = 10)
p + ggsave("./figures/7/Size_boxplot_perTissue.svg",width = 15, height = 10)

########## Part 3.1: Organoids species mixing sphere size and GFP integration ##########
#This part integrates size and GFP signal from sorted spheres with Seurat object of single cells 
#GFP signal gets normalized by sphere-size because the larger a sphere the higher the autofluorescence 

########## Prepare environment ##########
###Setting the working directory 
setwd("/mnt/khandler/R_projects/Sphere-sequencing/Sphere-seq_analysis/")

###Load packages and functions 
source("./functions_and_packages/1.Packages.R")
source("./functions_and_packages/3.Functions_sphere_size_GFP_integration.R")

###Load data 
org_mix <- readRDS("./data_files_generated/SpS_Organoid_mixing.Rda")
org_mix_ensemble <- readRDS("./data_files_generated/SpS_Organoid_mixing_ensembleIDs.Rda")

########## Linear model of standard sized beads ##########
#TOF from standard sized beads (60µm, 125µm and 175µm) was acquired and a linear model was generated 

###load acquired data of standard sized beads 
standard60 <- read.table("./data/20200824_60beads_for_standard_curve.txt", sep = "\t", 
                         header = TRUE, blank.lines.skip = FALSE, stringsAsFactors = FALSE)
#remove header 
standard60 <- standard60[1:(which(standard60$Id == "")[1] - 1),]
standard125 <- read.table("./data/20200824_125beads_for_standard_curve.txt", sep = "\t", 
                          header = TRUE, blank.lines.skip = FALSE, stringsAsFactors = FALSE)
standard125 <- standard125[1:(which(standard125$Id == "")[1] - 1),]
standard175 <- read.table("./data/20200824_175beads_for_standard_curve.txt", sep = "\t", 
                          header = TRUE, blank.lines.skip = FALSE, stringsAsFactors = FALSE)
standard175 <- standard175[1:(which(standard175$Id == "")[1] - 1),]

#concatenate individual acquisitions
bead_df <- rbind(standard60, standard125, standard175)
bead_df$bead_size <- c(rep(60, nrow(standard60)), rep(125, nrow(standard125)), rep(175, nrow(standard175)))

###fit a linear model with TOF 
lm.bead_tof <- lm(bead_size ~ TOF, data = bead_df)

########## Read data from biosorter outputs ##########
Plate1 <- read.table("./data/OrgMix_plate1_h_2.txt", sep = "\t", 
                     header = TRUE, blank.lines.skip = FALSE, stringsAsFactors = FALSE)
Plate1 <- Plate1[1:(which(Plate1$Id == "")[1] - 1),]
Plate2 <- read.table("./data/OrgMix_plate1_h.txt", sep = "\t", 
                     header = TRUE, blank.lines.skip = FALSE, stringsAsFactors = FALSE)
Plate2 <- Plate2[1:(which(Plate2$Id == "")[1] - 1),]
Plate3 <- read.table("./data/OrgMix_plate1_m.txt", sep = "\t", 
                     header = TRUE, blank.lines.skip = FALSE, stringsAsFactors = FALSE)
Plate3 <- Plate3[1:(which(Plate3$Id == "")[1] - 1),]
Plate4 <- read.table("./data/OrgMix_plate2_h.txt", sep = "\t", 
                     header = TRUE, blank.lines.skip = FALSE, stringsAsFactors = FALSE)
Plate4 <- Plate4[1:(which(Plate4$Id == "")[1] - 1),]
Plate5 <- read.table("./data/OrgMix_plate2_m.txt", sep = "\t", 
                     header = TRUE, blank.lines.skip = FALSE, stringsAsFactors = FALSE)
Plate5 <- Plate5[1:(which(Plate5$Id == "")[1] - 1),]
Plate6 <- read.table("./data/OrgMix_plate3_h.txt", sep = "\t", 
                     header = TRUE, blank.lines.skip = FALSE, stringsAsFactors = FALSE)
Plate6 <- Plate6[1:(which(Plate6$Id == "")[1] - 1),]
Plate7 <- read.table("./data/OrgMix_plate3_m.txt", sep = "\t", 
                     header = TRUE, blank.lines.skip = FALSE, stringsAsFactors = FALSE)
Plate7 <- Plate7[1:(which(Plate7$Id == "")[1] - 1),]

#add plate ID 
Plate1$plate_id <- "plate_1"
Plate2$plate_id <- "plate_1"
Plate3$plate_id <- "plate_1"
Plate4$plate_id <- "plate_2"
Plate5$plate_id <- "plate_2"
Plate6$plate_id <- "plate_3"
Plate7$plate_id <- "plate_3"

#merge files 
all_plates <- full_join(Plate1, Plate2)
all_plates <- full_join(all_plates,Plate3)
all_plates <- full_join(all_plates,Plate4)
all_plates <- full_join(all_plates,Plate5)
all_plates <- full_join(all_plates,Plate6)
all_plates <- full_join(all_plates,Plate7)

########## matching of biosorter output with Seurat object ##########
all_plates <- Biosorter_output_managing(all_plates,"")

#Seurat object with gene IDs 
Biosorter_data_seurat_integration(
  "./data_files_generated/SpS_Organoid_mixing.Rda",
  all_plates,
  "./data_files_generated/SpS_Organoid_mixing_BS.Rda")

#Seurat object with ensemble IDs 
Biosorter_data_seurat_integration(
  "./data_files_generated/SpS_Organoid_mixing_ensembleIDs.Rda",
  all_plates,
  "./data_files_generated/SpS_Organoid_mixing_ensembleIDs_BS.Rda")


########## plot GFP signal per sphere and species ID ##########
###load R object 
org_mix <- readRDS("./data_files_generated/SpS_Organoid_mixing_BS.Rda")

###annotated human and mouse wells based on sorting and based on prior knowledge that Bar19,36 and 106 contain 100% mouse cells 
#even though they were GFP negative (sorting error)
human_wells <- paste("Bar",c(1:18,20:35,37:48,97:105,107:144,194:240), sep = "")
mouse_wells <- paste("Bar",c(49:96,145:193,241:288,19,36,106), sep = "")

#add meta data column with well information 
org_mix$species_well <- NA
org_mix@meta.data <- org_mix@meta.data %>%
  mutate(species_well = case_when(
    sphere %in% human_wells ~ "human_well",
    sphere %in% mouse_wells ~ "mouse_well",
    TRUE ~ NA_character_))

###remove doublets and negatives 
spheres <- as.data.frame(table(org_mix$sphere))$Var1
spheres <- as.character(spheres)
spheres <- spheres[!spheres %in% c("Negative",
                                   "Doublet")]
Idents(org_mix) <- "sphere"
org_mix <- subset(org_mix, idents = spheres)

###plot boxplot 
org_mix_df <- as.data.frame(org_mix@meta.data)

org_mix_df2 <- org_mix_df[,c("species_well","GFP.norm","sphere")]
rownames(org_mix_df2) <- NULL
org_mix_df3 <- distinct(org_mix_df2)

p <- ggplot(org_mix_df3,aes(x = species_well,y = GFP.norm, fill = species_well)) +theme_classic() +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(position = position_jitter(seed = 1, width =0.4),size = 1) + 
  theme(axis.text = element_text(size = 30))  + scale_y_continuous(limits=c(0, 2.5)) + 
  ggtitle("Norm. GFP signal per cell from human and mouse wells") + xlab("species well") + 
  ylab("Normalized signal per sphere") + theme(axis.title= element_text(size = 25)) + 
  theme(plot.title = element_text(size = 25, face = "bold"))  + scale_colour_discrete("Species") +
  theme(legend.title = element_text(size = 30), legend.text = element_text(size = 30)) +
  theme(legend.title = element_text(size = 30), legend.text = element_text(size = 30))  +
  ggsignif::geom_signif(comparisons = list(c("human_well", "mouse_well")),
                        textsize=10, y_position = 2,
                        test = "wilcox.test",map_signif_level = c("***"=0.001,"**"=0.01,"*"=0.05))
p + ggsave("./figures/3.1/Green_boxplot_perSphere_human_mouse_wells.pdf",width = 12, height = 10)  
p + ggsave("./figures/3.1/Green_boxplot_perSphere_human_mouse_wells.svg",width = 12, height = 10)  






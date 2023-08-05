########## Part 2.3: Sphere size integration of liver samples ##########
#This part integrates sphere size with Seurat object
#Sphere size is calculated using acquired TOF measurements of sorted spheres, which is fit into a linear model from standard size beads TOF measurements

########## Prepare environment ##########
###Setting the working directory 
setwd("/mnt/khandler/R_projects/Sphere-sequencing/Sphere-seq_analysis/")

###Load packages and functions 
source("./functions_and_packages/1.Packages.R")
source("./functions_and_packages/3.Functions_sphere_size_GFP_integration.R")

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

###plot linear model 
p <- ggplot(bead_df, aes(bead_size, TOF)) + geom_point() +theme_classic() +
  theme(axis.text.y = element_text(size = 30), 
        axis.text.x = element_text(size = 30),
        axis.title.x = element_text(size = 25, margin=margin(5,0,0,0)), 
        axis.title.y = element_text(size = 25, margin=margin(0,5,0,0))) + 
  ggtitle("Linear regression model") + 
  theme(plot.title = element_text(size = 25, face = "bold", margin=margin(0,0,20,0))) + xlab("Bead size (µm)") +
  ylab("TOF") +
  stat_smooth(method = "lm", se=TRUE)  
p + ggsave("./figures/2.3/linear_model.pdf",width = 12, height = 10)
p + ggsave("./figures/2.3/linear_model.svg",width = 12, height = 10)

########## Read data from biosorter outputs #####
### Sample 1 
#read in acquired data 
S1_Plate1 <- read.table("./data/Liver_sample1_plate1_1_small.txt", sep = "\t", 
                             header = TRUE, blank.lines.skip = FALSE, stringsAsFactors = FALSE)
S1_Plate1 <- S1_Plate1[1:(which(S1_Plate1$Id == "")[1] - 1),]
S1_Plate2 <- read.table("./data/Liver_sample1_plate1_2_large.txt", sep = "\t", 
                             header = TRUE, blank.lines.skip = FALSE, stringsAsFactors = FALSE)
S1_Plate2 <- S1_Plate2[1:(which(S1_Plate2$Id == "")[1] - 1),]
S1_Plate3 <- read.table("./data/Liver_sample1_plate2_largeR19.txt", sep = "\t", 
                             header = TRUE, blank.lines.skip = FALSE, stringsAsFactors = FALSE)
S1_Plate3 <- S1_Plate3[1:(which(S1_Plate3$Id == "")[1] - 1),]
S1_Plate4 <- read.table("./data/Liver_sample1_plate3_smallR16.txt", sep = "\t", 
                             header = TRUE, blank.lines.skip = FALSE, stringsAsFactors = FALSE)
S1_Plate4 <- S1_Plate4[1:(which(S1_Plate4$Id == "")[1] - 1),]

#add plate ID 
S1_Plate1$plate_id <- paste0("plate_", c(rep(1, nrow(S1_Plate1))))
S1_Plate2$plate_id <- paste0("plate_", c(rep(1, nrow(S1_Plate2))))
S1_Plate3$plate_id <- paste0("plate_", c(rep(2, nrow(S1_Plate3))))
S1_Plate4$plate_id <- paste0("plate_", c(rep(3, nrow(S1_Plate4))))

#merge files
S1 <- full_join(S1_Plate1, S1_Plate2)
S1 <- full_join(S1,S1_Plate3)
S1 <- full_join(S1,S1_Plate4)

###Sample 2 
S2_Plate1 <- read.table("./data/Liver_sample2_plate1_largeHalf.txt", sep = "\t", 
                             header = TRUE, blank.lines.skip = FALSE, stringsAsFactors = FALSE)
S2_Plate1 <- S2_Plate1[1:(which(S2_Plate1$Id == "")[1] - 1),]
S2_Plate2 <- read.table("./data/Liver_sample2_plate1_largeHalf2.txt", sep = "\t", 
                             header = TRUE, blank.lines.skip = FALSE, stringsAsFactors = FALSE)
S2_Plate2 <- S2_Plate2[1:(which(S2_Plate2$Id == "")[1] - 1),]
S2_Plate3 <- read.table("./data/Liver_sample2_plate1_smallHalfLower.txt", sep = "\t", 
                             header = TRUE, blank.lines.skip = FALSE, stringsAsFactors = FALSE)
S2_Plate3 <- S2_Plate3[1:(which(S2_Plate3$Id == "")[1] - 1),]
S2_Plate4 <- read.table("./data/Liver_sample2_plate2_smallR16.txt", sep = "\t", 
                             header = TRUE, blank.lines.skip = FALSE, stringsAsFactors = FALSE)
S2_Plate4 <- S2_Plate4[1:(which(S2_Plate4$Id == "")[1] - 1),]
S2_Plate5 <- read.table("./data/Liver_sample2_plate3large_secondtrial_splilled.txt", sep = "\t", 
                             header = TRUE, blank.lines.skip = FALSE, stringsAsFactors = FALSE)
S2_Plate5 <- S2_Plate5[1:(which(S2_Plate5$Id == "")[1] - 1),]

S2_Plate1$plate_id <- paste0("plate_", c(rep(1, nrow(S2_Plate1))))
S2_Plate2$plate_id <- paste0("plate_", c(rep(1, nrow(S2_Plate2))))
S2_Plate3$plate_id <- paste0("plate_", c(rep(1, nrow(S2_Plate3))))
S2_Plate4$plate_id <- paste0("plate_", c(rep(2, nrow(S2_Plate4))))
S2_Plate5$plate_id <- paste0("plate_", c(rep(3, nrow(S2_Plate5))))

S2 <- full_join(S2_Plate1, S2_Plate2)
S2 <- full_join(S2,S2_Plate3)
#in plate 1 the plate was turned while pipetting the sphere BC, therefore switch Row and columns
#A=H, B=G, C=F, D=E, E=D, F=C, G=B, H=A
#1=12,2=11,3=10,4=9,5=8,6=7,7=6,8=5,9=4,10=3,11=2,12=1
numbers <- c(12:1,1:12,12:1,1:12,12:1,1:12,12:1,1:12)
S2$Column <- numbers
letter <- c(rep("H",12),rep("G",12),rep("F",12),rep("E",12),rep("D",12),rep("C",12),rep("B",12),rep("A",12))
S2$Row <- letter
S2 <- full_join(S2,S2_Plate4)
S2 <- full_join(S2,S2_Plate5)

###Sample 3
S3_Plate1 <- read.table("./data/Liver_sample3_Plate1largeR19_2.txt", sep = "\t", 
                             header = TRUE, blank.lines.skip = FALSE, stringsAsFactors = FALSE)
S3_Plate1 <- S3_Plate1[1:(which(S3_Plate1$Id == "")[1] - 1),]
S3_Plate2 <- read.table("./data/Liver_sample3_Plate1largeR25.txt", sep = "\t", 
                             header = TRUE, blank.lines.skip = FALSE, stringsAsFactors = FALSE)
S3_Plate2 <- S3_Plate2[1:(which(S3_Plate2$Id == "")[1] - 1),]
S3_Plate3 <- read.table("./data/Liver_sample3_Plate1small_R16.txt", sep = "\t", 
                             header = TRUE, blank.lines.skip = FALSE, stringsAsFactors = FALSE)
S3_Plate3 <- S3_Plate3[1:(which(S3_Plate3$Id == "")[1] - 1),]
S3_Plate4 <- read.table("./data/Liver_sample3_Plate2smallR16_2.txt", sep = "\t", 
                             header = TRUE, blank.lines.skip = FALSE, stringsAsFactors = FALSE)
S3_Plate4 <- S3_Plate4[1:(which(S3_Plate4$Id == "")[1] - 1),]
S3_Plate5 <- read.table("./data/Liver_sample3_Plate2smallR24.txt", sep = "\t", 
                             header = TRUE, blank.lines.skip = FALSE, stringsAsFactors = FALSE)
S3_Plate5 <- S3_Plate5[1:(which(S3_Plate5$Id == "")[1] - 1),]
S3_Plate6 <- read.table("./data/Liver_sample3_Plate3largeR25.txt", sep = "\t", 
                             header = TRUE, blank.lines.skip = FALSE, stringsAsFactors = FALSE)
S3_Plate6 <- S3_Plate6[1:(which(S3_Plate6$Id == "")[1] - 1),]
S3_Plate7 <- read.table("./data/Liver_sample3_Plate3smallR24_2.txt", sep = "\t", 
                             header = TRUE, blank.lines.skip = FALSE, stringsAsFactors = FALSE)
S3_Plate7 <- S3_Plate7[1:(which(S3_Plate7$Id == "")[1] - 1),]

S3_Plate1$plate_id <- paste0("plate_", c(rep(1, nrow(S3_Plate1))))
S3_Plate2$plate_id <- paste0("plate_", c(rep(1, nrow(S3_Plate2))))
S3_Plate3$plate_id <- paste0("plate_", c(rep(1, nrow(S3_Plate3))))
S3_Plate4$plate_id <- paste0("plate_", c(rep(2, nrow(S3_Plate4))))
S3_Plate5$plate_id <- paste0("plate_", c(rep(2, nrow(S3_Plate5))))
S3_Plate6$plate_id <- paste0("plate_", c(rep(3, nrow(S3_Plate6))))
S3_Plate7$plate_id <- paste0("plate_", c(rep(3, nrow(S3_Plate7))))

S3 <- full_join(S3_Plate1, S3_Plate2)
S3 <- full_join(S3,S3_Plate3)
S3 <- full_join(S3,S3_Plate4)
S3 <- full_join(S3,S3_Plate5)
S3 <- full_join(S3,S3_Plate6)
S3 <- full_join(S3,S3_Plate7)

###Sample 4
S4_Plate1 <- read.table("./data/Liver_sample4_Plate1R26.txt", sep = "\t", 
                             header = TRUE, blank.lines.skip = FALSE, stringsAsFactors = FALSE)
S4_Plate1 <- S4_Plate1[1:(which(S4_Plate1$Id == "")[1] - 1),]
S4_Plate2 <- read.table("./data/Liver_sample4_Plate2R29_2.txt", sep = "\t", 
                             header = TRUE, blank.lines.skip = FALSE, stringsAsFactors = FALSE)
S4_Plate2 <- S4_Plate2[1:(which(S4_Plate2$Id == "")[1] - 1),]
S4_Plate3 <- read.table("./data/Liver_sample4_Plate2R29.txt", sep = "\t", 
                             header = TRUE, blank.lines.skip = FALSE, stringsAsFactors = FALSE)
S4_Plate3 <- S4_Plate3[1:(which(S4_Plate3$Id == "")[1] - 1),]
S4_Plate4 <- read.table("./data/Liver_sample4_Plate3_2R26.txt", sep = "\t", 
                             header = TRUE, blank.lines.skip = FALSE, stringsAsFactors = FALSE)
S4_Plate4 <- S4_Plate4[1:(which(S4_Plate4$Id == "")[1] - 1),]
S4_Plate5 <- read.table("./data/Liver_sample4_Plate3R29.txt", sep = "\t", 
                             header = TRUE, blank.lines.skip = FALSE, stringsAsFactors = FALSE)
S4_Plate5 <- S4_Plate5[1:(which(S4_Plate5$Id == "")[1] - 1),]

S4_Plate1$plate_id <- paste0("plate_", c(rep(1, nrow(S4_Plate1))))
S4_Plate2$plate_id <- paste0("plate_", c(rep(1, nrow(S4_Plate2))))
S4_Plate3$plate_id <- paste0("plate_", c(rep(1, nrow(S4_Plate3))))
S4_Plate4$plate_id <- paste0("plate_", c(rep(2, nrow(S4_Plate4))))
S4_Plate5$plate_id <- paste0("plate_", c(rep(3, nrow(S4_Plate5))))

S4 <- full_join(S4_Plate1, S4_Plate2)
S4 <- full_join(S4,S4_Plate3)
S4 <- full_join(S4,S4_Plate4)
S4 <- full_join(S4,S4_Plate5)

###Sample 5
S5_Plate1 <- read.table("./data/Liver_sample5_plate1_R16.txt", sep = "\t", 
                             header = TRUE, blank.lines.skip = FALSE, stringsAsFactors = FALSE)
S5_Plate1 <- S5_Plate1[1:(which(S5_Plate1$Id == "")[1] - 1),]
S5_Plate2 <- read.table("./data/Liver_sample5_plate2_2R24.txt", sep = "\t", 
                             header = TRUE, blank.lines.skip = FALSE, stringsAsFactors = FALSE)
S5_Plate2 <- S5_Plate2[1:(which(S5_Plate2$Id == "")[1] - 1),]
S5_Plate3 <- read.table("./data/Liver_sample5_plate2_R16.txt", sep = "\t", 
                             header = TRUE, blank.lines.skip = FALSE, stringsAsFactors = FALSE)
S5_Plate3 <- S5_Plate3[1:(which(S5_Plate3$Id == "")[1] - 1),]
S5_Plate4 <- read.table("./data/Liver_sample5_plate3_2R29.txt", sep = "\t", 
                             header = TRUE, blank.lines.skip = FALSE, stringsAsFactors = FALSE)
S5_Plate4 <- S5_Plate4[1:(which(S5_Plate4$Id == "")[1] - 1),]
S5_Plate5 <- read.table("./data/Liver_sample5_plate3_3R25.txt", sep = "\t", 
                             header = TRUE, blank.lines.skip = FALSE, stringsAsFactors = FALSE)
S5_Plate5 <- S5_Plate5[1:(which(S5_Plate5$Id == "")[1] - 1),]
S5_Plate6 <- read.table("./data/Liver_sample5_plate3_R30.txt", sep = "\t", 
                             header = TRUE, blank.lines.skip = FALSE, stringsAsFactors = FALSE)
S5_Plate6 <- S5_Plate6[1:(which(S5_Plate6$Id == "")[1] - 1),]

S5_Plate1$plate_id <- paste0("plate_", c(rep(1, nrow(S5_Plate1))))
S5_Plate2$plate_id <- paste0("plate_", c(rep(1, nrow(S5_Plate2))))
S5_Plate3$plate_id <- paste0("plate_", c(rep(1, nrow(S5_Plate3))))
S5_Plate4$plate_id <- paste0("plate_", c(rep(2, nrow(S5_Plate4))))
S5_Plate5$plate_id <- paste0("plate_", c(rep(3, nrow(S5_Plate5))))
S5_Plate6$plate_id <- paste0("plate_", c(rep(3, nrow(S5_Plate6))))

S5 <- full_join(S5_Plate1, S5_Plate2)
S5 <- full_join(S5,S5_Plate3)
S5 <- full_join(S5,S5_Plate4)
S5 <- full_join(S5,S5_Plate5)
S5 <- full_join(S5,S5_Plate6)

###Sample 6
S6_Plate1 <- read.table("./data/Liver_sample6_plate1R19.txt", sep = "\t", 
                           header = TRUE, blank.lines.skip = FALSE, stringsAsFactors = FALSE)
S6_Plate1 <- S6_Plate1[1:(which(S6_Plate1$Id == "")[1] - 1),]
S6_Plate2 <- read.table("./data/Liver_sample6_plate2R16.txt", sep = "\t", 
                           header = TRUE, blank.lines.skip = FALSE, stringsAsFactors = FALSE)
S6_Plate2 <- S6_Plate2[1:(which(S6_Plate2$Id == "")[1] - 1),]
S6_Plate3 <- read.table("./data/Liver_sample6_plate2R19.txt", sep = "\t", 
                           header = TRUE, blank.lines.skip = FALSE, stringsAsFactors = FALSE)
S6_Plate3 <- S6_Plate3[1:(which(S6_Plate3$Id == "")[1] - 1),]
S6_Plate4 <- read.table("./data/Liver_sample6_plate2R26.txt", sep = "\t", 
                           header = TRUE, blank.lines.skip = FALSE, stringsAsFactors = FALSE)
S6_Plate4 <- S6_Plate4[1:(which(S6_Plate4$Id == "")[1] - 1),]
S6_Plate5 <- read.table("./data/Liver_sample6_plate2R192.txt", sep = "\t", 
                           header = TRUE, blank.lines.skip = FALSE, stringsAsFactors = FALSE)
S6_Plate5 <- S6_Plate5[1:(which(S6_Plate5$Id == "")[1] - 1),]
S6_Plate6 <- read.table("./data/Liver_sample6_plate3R16.txt", sep = "\t", 
                           header = TRUE, blank.lines.skip = FALSE, stringsAsFactors = FALSE)
S6_Plate6 <- S6_Plate6[1:(which(S6_Plate6$Id == "")[1] - 1),]

S6_Plate1$plate_id <- paste0("plate_", c(rep(1, nrow(S6_Plate1))))
S6_Plate2$plate_id <- paste0("plate_", c(rep(2, nrow(S6_Plate2))))
S6_Plate3$plate_id <- paste0("plate_", c(rep(2, nrow(S6_Plate3))))
S6_Plate4$plate_id <- paste0("plate_", c(rep(3, nrow(S6_Plate4))))
S6_Plate5$plate_id <- paste0("plate_", c(rep(3, nrow(S6_Plate5))))
S6_Plate6$plate_id <- paste0("plate_", c(rep(3, nrow(S6_Plate6))))

S6 <- full_join(S6_Plate1, S6_Plate2)
S6 <- full_join(S6,S6_Plate3)
S6 <- full_join(S6,S6_Plate4)
S6 <- full_join(S6,S6_Plate5)
S6 <- full_join(S6,S6_Plate6)

###Sample 7
S7_Plate1 <- read.table("./data/Liver_sample7_plate1R16.txt", sep = "\t", 
                             header = TRUE, blank.lines.skip = FALSE, stringsAsFactors = FALSE)
S7_Plate1 <- S7_Plate1[1:(which(S7_Plate1$Id == "")[1] - 1),]
S7_Plate2 <- read.table("./data/Liver_sample7_plate1R16_2.txt", sep = "\t", 
                             header = TRUE, blank.lines.skip = FALSE, stringsAsFactors = FALSE)
S7_Plate2 <- S7_Plate2[1:(which(S7_Plate2$Id == "")[1] - 1),]
S7_Plate3 <- read.table("./data/Liver_sample7_plate2R19.txt", sep = "\t", 
                             header = TRUE, blank.lines.skip = FALSE, stringsAsFactors = FALSE)
S7_Plate3 <- S7_Plate3[1:(which(S7_Plate3$Id == "")[1] - 1),]
S7_Plate4 <- read.table("./data/Liver_sample7_plate3R16.txt", sep = "\t", 
                             header = TRUE, blank.lines.skip = FALSE, stringsAsFactors = FALSE)
S7_Plate4 <- S7_Plate4[1:(which(S7_Plate4$Id == "")[1] - 1),]
S7_Plate5 <- read.table("./data/Liver_sample7_plate3R19.txt", sep = "\t", 
                             header = TRUE, blank.lines.skip = FALSE, stringsAsFactors = FALSE)
S7_Plate5 <- S7_Plate5[1:(which(S7_Plate5$Id == "")[1] - 1),]

S7_Plate1$plate_id <- paste0("plate_", c(rep(1, nrow(S7_Plate1))))
S7_Plate2$plate_id <- paste0("plate_", c(rep(1, nrow(S7_Plate2))))
S7_Plate3$plate_id <- paste0("plate_", c(rep(2, nrow(S7_Plate3))))
S7_Plate4$plate_id <- paste0("plate_", c(rep(3, nrow(S7_Plate4))))
S7_Plate5$plate_id <- paste0("plate_", c(rep(3, nrow(S7_Plate5))))

S7 <- full_join(S7_Plate1, S7_Plate2)
S7 <- full_join(S7,S7_Plate3)
S7 <- full_join(S7,S7_Plate4)
S7 <- full_join(S7,S7_Plate5)

###Sample 8
S8_Plate1 <- read.table("./data/Liver_sample8_plate1.txt", sep = "\t", 
                             header = TRUE, blank.lines.skip = FALSE, stringsAsFactors = FALSE)
S8_Plate1 <- S8_Plate1[1:(which(S8_Plate1$Id == "")[1] - 1),]
S8_Plate2 <- read.table("./data/Liver_sample8_plate2.txt", sep = "\t", 
                             header = TRUE, blank.lines.skip = FALSE, stringsAsFactors = FALSE)
S8_Plate2 <- S8_Plate2[1:(which(S8_Plate2$Id == "")[1] - 1),]
S8_Plate3 <- read.table("./data/Liver_sample8_plate3.txt", sep = "\t", 
                             header = TRUE, blank.lines.skip = FALSE, stringsAsFactors = FALSE)
S8_Plate3 <- S8_Plate3[1:(which(S8_Plate3$Id == "")[1] - 1),]

S8_Plate1$plate_id <- paste0("plate_", c(rep(1, nrow(S8_Plate1))))
S8_Plate2$plate_id <- paste0("plate_", c(rep(2, nrow(S8_Plate2))))
S8_Plate3$plate_id <- paste0("plate_", c(rep(3, nrow(S8_Plate3))))

S8 <- full_join(S8_Plate1, S8_Plate2)
S8 <- full_join(S8,S8_Plate3)

###Sample 9
S9_Plate1 <- read.table("./data/Liver_sample9_plate1.txt", sep = "\t", 
                             header = TRUE, blank.lines.skip = FALSE, stringsAsFactors = FALSE)
S9_Plate1 <- S9_Plate1[1:(which(S9_Plate1$Id == "")[1] - 1),]
S9_Plate2 <- read.table("./data/Liver_sample9_plate2.txt", sep = "\t", 
                             header = TRUE, blank.lines.skip = FALSE, stringsAsFactors = FALSE)
S9_Plate2 <- S9_Plate2[1:(which(S9_Plate2$Id == "")[1] - 1),]
S9_Plate3 <- read.table("./data/Liver_sample9_plate3.txt", sep = "\t", 
                             header = TRUE, blank.lines.skip = FALSE, stringsAsFactors = FALSE)
S9_Plate3 <- S9_Plate3[1:(which(S9_Plate3$Id == "")[1] - 1),]

S9_Plate1$plate_id <- paste0("plate_", c(rep(1, nrow(S9_Plate1))))
S9_Plate2$plate_id <- paste0("plate_", c(rep(2, nrow(S9_Plate2))))
S9_Plate3$plate_id <- paste0("plate_", c(rep(3, nrow(S9_Plate3))))

S9 <- full_join(S9_Plate1, S9_Plate2)
S9 <- full_join(S9,S9_Plate3)

###Sample 10
S10_Plate1 <- read.table("./data/Liver_sample10_plate1R16.txt", sep = "\t", 
                        header = TRUE, blank.lines.skip = FALSE, stringsAsFactors = FALSE)
S10_Plate1 <- S10_Plate1[1:(which(S10_Plate1$Id == "")[1] - 1),]
S10_Plate2 <- read.table("./data/Liver_sample10_plate2R16.txt", sep = "\t", 
                        header = TRUE, blank.lines.skip = FALSE, stringsAsFactors = FALSE)
S10_Plate2 <- S10_Plate2[1:(which(S10_Plate2$Id == "")[1] - 1),]
S10_Plate3 <- read.table("./data/Liver_sample10_plate2R19.txt", sep = "\t", 
                        header = TRUE, blank.lines.skip = FALSE, stringsAsFactors = FALSE)
S10_Plate3 <- S10_Plate3[1:(which(S10_Plate3$Id == "")[1] - 1),]
S10_Plate4 <- read.table("./data/Liver_sample10_plate3R19.txt", sep = "\t", 
                        header = TRUE, blank.lines.skip = FALSE, stringsAsFactors = FALSE)
S10_Plate4 <- S10_Plate4[1:(which(S10_Plate4$Id == "")[1] - 1),]

S10_Plate1$plate_id <- paste0("plate_", c(rep(1, nrow(S10_Plate1))))
S10_Plate2$plate_id <- paste0("plate_", c(rep(2, nrow(S10_Plate2))))
S10_Plate3$plate_id <- paste0("plate_", c(rep(2, nrow(S10_Plate3))))
S10_Plate4$plate_id <- paste0("plate_", c(rep(3, nrow(S10_Plate4))))

S10 <- full_join(S10_Plate1, S10_Plate2)
S10 <- full_join(S10,S10_Plate3)
S10 <- full_join(S10,S10_Plate4)


########## Sphere size prediction and integration into Seurat object ##########
###map output of Biosorter with sphere (MULTI-seq) ID 
S1_2 <- Biosorter_output_managing(S1,"_1")
S2_2 <- Biosorter_output_managing(S2,"_2")
S3_2 <- Biosorter_output_managing(S3,"_3")
S4_2 <- Biosorter_output_managing(S4,"_4")
S5_2 <- Biosorter_output_managing(S5,"_5")
S6_2 <- Biosorter_output_managing(S6,"_6")
S7_2 <- Biosorter_output_managing(S7,"_7")
S8_2 <- Biosorter_output_managing(S8,"_8")
S9_2 <- Biosorter_output_managing(S9,"_9")
S10_2 <- Biosorter_output_managing(S10,"_10")

#merge all samples 
biosorter_all <- full_join(S1_2, S2_2)
biosorter_all <- full_join(biosorter_all,S3_2)
biosorter_all <- full_join(biosorter_all,S4_2)
biosorter_all <- full_join(biosorter_all,S5_2)
biosorter_all <- full_join(biosorter_all,S6_2)
biosorter_all <- full_join(biosorter_all,S7_2)
biosorter_all <- full_join(biosorter_all,S8_2)
biosorter_all <- full_join(biosorter_all,S9_2)
biosorter_all <- full_join(biosorter_all,S10_2)

###Integrate sphere-size information with Seurat object 
Biosorter_data_seurat_integration(
  "./data_files_generated/LiverMerged_afterBC_anno.Rda",
  biosorter_all,
  "./data_files_generated/LiverMerged_afterBC_anno_BS.Rda")

###Remove negative and doublets 
liverSpS <- readRDS("./data_files_generated/LiverMerged_afterBC_anno_BS.Rda")
spheres <- as.data.frame(table(liverSpS$sphere))$Var1
spheres <- as.character(spheres)
spheres <- spheres[!spheres %in% c("Negative_1","Negative_2","Negative_3","Negative_4","Negative_5",
                                   "Negative_6","Negative_7","Negative_8","Negative_9","Negative_10",
                                   "Doublet_1","Doublet_2","Doublet_3","Doublet_4","Doublet_5",
                                   "Doublet_6","Doublet_7","Doublet_8","Doublet_9","Doublet_10")]
Idents(liverSpS) <- "sphere"
liverSpS <- subset(liverSpS, idents = spheres)

###save object 
saveRDS(liverSpS,file = "./data_files_generated/LiverMerged_afterBC_anno_BS.Rda")

########## Plot sphere size distribution per sample ##########
###subset each sample 
#M1=S1, M3=S2, M4=S3, 4M1=S4, 4M2=S5, M5=S6, 6M1=S7, 6M2=S8, 6M3=S9, WTCre=S10
Idents(liverSpS) <- "orig.ident"
subS1 <- subset(liverSpS, idents = c("M1"))
subS2 <- subset(liverSpS, idents = c("M3"))
subS3 <- subset(liverSpS, idents = c("M4"))
subS4 <- subset(liverSpS, idents = c("4M1"))
subS5 <- subset(liverSpS, idents = c("4M2"))
subS6 <- subset(liverSpS, idents = c("M5"))
subS7 <- subset(liverSpS, idents = c("6M1"))
subS8 <- subset(liverSpS, idents = c("6M2"))
subS9 <- subset(liverSpS, idents = c("6M3"))
subS10 <- subset(liverSpS, idents = c("WTCre"))

###Make data frames of sphere sizes per sample
dfS1_size <- BS_df_for_boxplot_per_sample(subS1,"sphere","sphere_size","S1")
dfS2_size <- BS_df_for_boxplot_per_sample(subS2,"sphere","sphere_size","S2")
dfS3_size <- BS_df_for_boxplot_per_sample(subS3,"sphere","sphere_size","S3")
dfS4_size <- BS_df_for_boxplot_per_sample(subS4,"sphere","sphere_size","S4")
dfS5_size <- BS_df_for_boxplot_per_sample(subS5,"sphere","sphere_size","S5")
dfS6_size <- BS_df_for_boxplot_per_sample(subS6,"sphere","sphere_size","S6")
dfS7_size <- BS_df_for_boxplot_per_sample(subS7,"sphere","sphere_size","S7")
dfS8_size <- BS_df_for_boxplot_per_sample(subS8,"sphere","sphere_size","S8")
dfS9_size <- BS_df_for_boxplot_per_sample(subS9,"sphere","sphere_size","S9")
dfS10_size <- BS_df_for_boxplot_per_sample(subS10,"sphere","sphere_size","S10")

#merge data frames
size_merged_df <- rbind(dfS1_size,dfS2_size)
size_merged_df <- rbind(size_merged_df,dfS3_size)
size_merged_df <- rbind(size_merged_df,dfS4_size)
size_merged_df <- rbind(size_merged_df,dfS5_size)
size_merged_df <- rbind(size_merged_df,dfS6_size)
size_merged_df <- rbind(size_merged_df,dfS7_size)
size_merged_df <- rbind(size_merged_df,dfS8_size)
size_merged_df <- rbind(size_merged_df,dfS9_size)
size_merged_df <- rbind(size_merged_df,dfS10_size)

###Plot in boxplot 
p <- ggplot(size_merged_df,aes(x = sample,y = sphere_size, fill = sample)) +theme_classic() +
  geom_boxplot() +
  geom_jitter(position = position_jitter(seed = 1, width =0.4),size = 1) + 
  theme(axis.text = element_text(size = 30))  +
  ggtitle("Sphere size per sample and sphere") + xlab("Sample") + 
  ylab("Size (µm)") + theme(axis.title= element_text(size = 25)) + 
  theme(plot.title = element_text(size = 25, face = "bold")) + 
  theme(legend.title = element_text(size = 30), legend.text = element_text(size = 30)) 
p + ggsave("./figures/2.3/Size_boxplot_perSample.pdf",width = 15, height = 10)
p + ggsave("./figures/2.3/Size_boxplot_perSample.svg",width = 15, height = 10)



########## Part 2: Pre-processing and feature area integration ##########
#This part does pre-processing on single cell object after cell segmentation 
#then x and y coordinates of manually drawn feature areas (pv, cv and mets) were integrated with x and y coordinates of Seurat object of segmented single cells 

########## Prepare environment ##########
###Setting the working directory 
setwd("/mnt/khandler/R_projects/Sphere-sequencing/Highly_multiplexed_FISH_analysis/")

###Load packages and functions 
source("./functions_and_packages/1.Packages.R")
source("./functions_and_packages/2.Functions_ImageJ_connection.R")

###Load data 
resolveSCE <- readRDS("./data/SCE.RDS")

########## QC ##########
resolveSCE <- addPerCellQCMetrics(resolveSCE)
###Add which ones have visible metastasis 
resolveSCE$Mets <- grepl("A1|B1",resolveSCE$Slide)

###remove low quality cells 
#remove everything with 0 counts 
keep <- resolveSCE$sum > 0
keep <- keep & resolveSCE$detected > 0 
###remove nuclei that are too large (4MAD from median)
resolveSCE$keep <- keep & !isOutlier(resolveSCE$Area,nmads=4)
resolveSCE <- resolveSCE[,keep]

###add log counts 
counts <- assay(resolveSCE, "counts")
libsizes <- colSums(counts)
size.factors <- libsizes/mean(libsizes)
logcounts(resolveSCE) <- log2(t(t(counts)/size.factors) + 1)

###make Seurat object
resolveS <- as.Seurat(resolveSCE)

########## match X and Y coordinate with manually drawn areas ##########
###combine X and Y coordinate to one column 
resolveS$imageJ_coordinate <- NA
resolveS$imageJ_coordinate <- paste0(round(resolveS$X),",", round(resolveS$Y))

###make separate objects per slide 
Idents(resolveS) <- "Slide"
#Samples with visible metastasis
A1_1 <- subset(resolveS, idents = "Slide1_A1-1")
B1_1 <- subset(resolveS, idents = "Slide1_B1-1")
B1_2 <- subset(resolveS, idents = "Slide1_B1-2")

#Samples with no visible metastasis
A2_1 <- subset(resolveS, idents = "Slide1_A2-1")
A2_2 <- subset(resolveS, idents = "Slide1_A2-2")
B2_1 <- subset(resolveS, idents = "Slide1_B2-1")

###A1_1
#pv surrounding mets 
pv_sm_1 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/A1_1/A1_1_pv_sm_5_1.csv")
pv_sm_2 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/A1_1/A1_1_pv_sm_5_2.csv")
pv_sm_3 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/A1_1/A1_1_pv_sm_3_1.csv")
pv_sm_4 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/A1_1/A1_1_pv_sm_3_2.csv")
pv_sm_5 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/A1_1/A1_1_pv_sm_4_1.csv")
pv_sm_6 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/A1_1/A1_1_pv_sm_4_3_1.csv")

#cv surrounding mets 
cv_sm_1 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/A1_1/A1_1_cv_sm_1_1.csv")
cv_sm_2 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/A1_1/A1_1_cv_sm_1_2.csv")
cv_sm_3 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/A1_1/A1_1_cv_sm_1_3_1.csv")
cv_sm_4 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/A1_1/A1_1_cv_sm_1_3.csv")
cv_sm_5 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/A1_1/A1_1_cv_sm_1_4.csv")
cv_sm_6 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/A1_1/A1_1_cv_sm_1_5.csv")
cv_sm_7 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/A1_1/A1_1_cv_sm_2_1.csv")
cv_sm_8 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/A1_1/A1_1_cv_sm_3_1.csv")
cv_sm_9 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/A1_1/A1_1_cv_sm_4_1.csv")
cv_sm_10 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/A1_1/A1_1_cv_sm_4_2.csv")
cv_sm_11 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/A1_1/A1_1_cv_sm_4_3_1.csv")
cv_sm_12 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/A1_1/A1_1_cv_sm_5_1.csv")
cv_sm_13 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/A1_1/A1_1_cv_sm_5_2.csv")

#metastasis 
m1 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/A1_1/A1_1_mets1.csv")
m2 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/A1_1/A1_1_mets2.csv")
m3 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/A1_1/A1_1_mets3.csv")
m4 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/A1_1/A1_1_mets4.csv")
m5 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/A1_1/A1_1_mets5.csv")

#cv no metastasis 
cv_nm_1 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/A1_1/A1_1_cv_nm_1.csv")
cv_nm_2 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/A1_1/A1_1_cv_nm_2.csv")
cv_nm_3 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/A1_1/A1_1_cv_nm_3.csv")
cv_nm_4 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/A1_1/A1_1_cv_nm_4.csv")
cv_nm_5 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/A1_1/A1_1_cv_nm_5.csv")
cv_nm_6 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/A1_1/A1_1_cv_nm_6.csv")
cv_nm_7 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/A1_1/A1_1_cv_nm_7.csv")
cv_nm_8 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/A1_1/A1_1_cv_nm_8.csv")
cv_nm_9 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/A1_1/A1_1_cv_nm_9.csv")
cv_nm_10 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/A1_1/A1_1_cv_nm_10.csv")
cv_nm_11 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/A1_1/A1_1_cv_nm_11.csv")

#pv no metastasis 
pv_nm_1 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/A1_1/A1_1_pv_nm_1.csv")
pv_nm_2 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/A1_1/A1_1_pv_nm_2.csv")
pv_nm_3 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/A1_1/A1_1_pv_nm_3.csv")
pv_nm_4 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/A1_1/A1_1_pv_nm_4.csv")
pv_nm_5 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/A1_1/A1_1_pv_nm_5.csv")
pv_nm_6 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/A1_1/A1_1_pv_nm_6.csv")
pv_nm_7 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/A1_1/A1_1_pv_nm_7.csv")
pv_nm_8 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/A1_1/A1_1_pv_nm_8.csv")
pv_nm_9 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/A1_1/A1_1_pv_nm_9.csv")
pv_nm_10 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/A1_1/A1_1_pv_nm_10.csv")
pv_nm_11 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/A1_1/A1_1_pv_nm_11.csv")


A1_1_mets_pv_sm <- c(pv_sm_1, pv_sm_2,pv_sm_3,pv_sm_4, pv_sm_5,pv_sm_6)
A1_1_mets_cv_sm <- c(cv_sm_1,cv_sm_2,cv_sm_3,cv_sm_4,cv_sm_5,cv_sm_6,cv_sm_7,cv_sm_8,cv_sm_9,cv_sm_10,
                     cv_sm_11,cv_sm_12,cv_sm_13)
A1_1_mets_pv_nm <- c(pv_nm_1,pv_nm_2,pv_nm_3,pv_nm_4,pv_nm_5,pv_nm_6,pv_nm_7,pv_nm_8,pv_nm_9,pv_nm_10,
                     pv_nm_11)
A1_1_mets_cv_nm <- c(cv_nm_1,cv_nm_2,cv_nm_3,cv_nm_4,cv_nm_5,cv_nm_6,cv_nm_7,cv_nm_8,cv_nm_9,cv_nm_10,
                     cv_nm_11)
A1_1_mets <- c(m1,m2,m3,m4,m5)

#add spatial info 
A1_1$spatial_feature <- NA
A1_1@meta.data <- A1_1@meta.data %>%
  mutate(spatial_feature = case_when(
    imageJ_coordinate %in% A1_1_mets_pv_sm ~ "pv_sm",
    imageJ_coordinate %in% A1_1_mets_cv_sm ~ "cv_sm",
    imageJ_coordinate %in% A1_1_mets_pv_nm ~ "pv_nm",
    imageJ_coordinate %in% A1_1_mets_cv_nm ~ "cv_nm",
    imageJ_coordinate %in% A1_1_mets ~ "mets",
    TRUE ~ NA_character_))

#remove NA 
Idents(A1_1) <- "spatial_feature"
A1_1 <- subset(A1_1, idents = c("cv_nm", "cv_sm" , "mets", "pv_nm" ,"pv_sm"))

#add spatial info id number
A1_1$spatial_feature_number <- NA
A1_1@meta.data <- A1_1@meta.data %>%
  mutate(spatial_feature_number = case_when(
    imageJ_coordinate %in% pv_sm_1 ~ "s1_pv_sm_1",
    imageJ_coordinate %in% pv_sm_2 ~ "s1_pv_sm_2",
    imageJ_coordinate %in% pv_sm_3 ~ "s1_pv_sm_3",
    imageJ_coordinate %in% pv_sm_4 ~ "s1_pv_sm_4",
    imageJ_coordinate %in% pv_sm_5 ~ "s1_pv_sm_5",
    imageJ_coordinate %in% pv_sm_6 ~ "s1_pv_sm_6",
    imageJ_coordinate %in% cv_sm_1 ~ "s1_cv_sm_1",
    imageJ_coordinate %in% cv_sm_2 ~ "s1_cv_sm_2",
    imageJ_coordinate %in% cv_sm_3 ~ "s1_cv_sm_3",
    imageJ_coordinate %in% cv_sm_4 ~ "s1_cv_sm_4",
    imageJ_coordinate %in% cv_sm_5 ~ "s1_cv_sm_5",
    imageJ_coordinate %in% cv_sm_6 ~ "s1_cv_sm_6",
    imageJ_coordinate %in% cv_sm_7 ~ "s1_cv_sm_7",
    imageJ_coordinate %in% cv_sm_8 ~ "s1_cv_sm_8",
    imageJ_coordinate %in% cv_sm_9 ~ "s1_cv_sm_9",
    imageJ_coordinate %in% cv_sm_10 ~ "s1_cv_sm_10",
    imageJ_coordinate %in% cv_sm_11 ~ "s1_cv_sm_11",
    imageJ_coordinate %in% cv_sm_12 ~ "s1_cv_sm_12",
    imageJ_coordinate %in% cv_sm_13 ~ "s1_cv_sm_13",
    imageJ_coordinate %in% pv_nm_1 ~ "s1_pv_nm_1",
    imageJ_coordinate %in% pv_nm_2 ~ "s1_pv_nm_1",
    imageJ_coordinate %in% pv_nm_3 ~ "s1_pv_nm_3",
    imageJ_coordinate %in% pv_nm_4 ~ "s1_pv_nm_4",
    imageJ_coordinate %in% pv_nm_5 ~ "s1_pv_nm_5",
    imageJ_coordinate %in% pv_nm_6 ~ "s1_pv_nm_6",
    imageJ_coordinate %in% pv_nm_7 ~ "s1_pv_nm_7",
    imageJ_coordinate %in% pv_nm_8 ~ "s1_pv_nm_8",
    imageJ_coordinate %in% pv_nm_9 ~ "s1_pv_nm_9",
    imageJ_coordinate %in% pv_nm_10 ~ "s1_pv_nm_10",
    imageJ_coordinate %in% pv_nm_11 ~ "s1_pv_nm_11",
    imageJ_coordinate %in% cv_nm_1 ~ "s1_cv_nm_1",
    imageJ_coordinate %in% cv_nm_2 ~ "s1_cv_nm_2",
    imageJ_coordinate %in% cv_nm_3 ~ "s1_cv_nm_3",
    imageJ_coordinate %in% cv_nm_4 ~ "s1_cv_nm_4",
    imageJ_coordinate %in% cv_nm_5 ~ "s1_cv_nm_5",
    imageJ_coordinate %in% cv_nm_6 ~ "s1_cv_nm_6",
    imageJ_coordinate %in% cv_nm_7 ~ "s1_cv_nm_7",
    imageJ_coordinate %in% cv_nm_8 ~ "s1_cv_nm_8",
    imageJ_coordinate %in% cv_nm_9 ~ "s1_cv_nm_9",
    imageJ_coordinate %in% cv_nm_10 ~ "s1_cv_nm_10",
    imageJ_coordinate %in% cv_nm_10 ~ "s1_cv_nm_11",
    imageJ_coordinate %in% m1 ~ "s1_m1",
    imageJ_coordinate %in% m2 ~ "s1_m2",
    imageJ_coordinate %in% m3 ~ "s1_m3",
    imageJ_coordinate %in% m4 ~ "s1_m4",
    imageJ_coordinate %in% m5 ~ "s1_m5",
    TRUE ~ NA_character_))

##save objects 
saveRDS(A1_1, file = "./data_files_generated/mets_A1_1.rds")

##A2_1
#pv 
pv_1 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/A2_1/A2_1_pv_1.csv")
pv_2 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/A2_1/A2_1_pv_2.csv")
pv_3 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/A2_1/A2_1_pv_3.csv")
pv_4 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/A2_1/A2_1_pv_4.csv")
pv_5 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/A2_1/A2_1_pv_5.csv")
pv_6 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/A2_1/A2_1_pv_6.csv")

#cv 
cv_1 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/A2_1/A2_1_cv_1.csv")
cv_2 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/A2_1/A2_1_cv_2.csv")
cv_3 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/A2_1/A2_1_cv_3.csv")
cv_4 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/A2_1/A2_1_cv_4.csv")
cv_5 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/A2_1/A2_1_cv_5.csv")
cv_6 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/A2_1/A2_1_cv_6.csv")
cv_7 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/A2_1/A2_1_cv_7.csv")


A2_1_noMets_pv <- c(pv_1,pv_2,pv_3,pv_4,pv_5,pv_6)
A2_1_noMets_cv <- c(cv_1,cv_2,cv_3,cv_4,cv_5,cv_6,cv_7)

#add spatial info 
A2_1$spatial_feature <- NA
A2_1@meta.data <- A2_1@meta.data %>%
  mutate(spatial_feature = case_when(
    imageJ_coordinate %in% A2_1_noMets_pv ~ "pv",
    imageJ_coordinate %in% A2_1_noMets_cv ~ "cv",
    TRUE ~ NA_character_))

#remove NA 
Idents(A2_1) <- "spatial_feature"
A2_1 <- subset(A2_1, idents = c("cv","pv"))

#add spatial info id number
A2_1$spatial_feature_number <- NA
A2_1@meta.data <- A2_1@meta.data %>%
  mutate(spatial_feature_number = case_when(
    imageJ_coordinate %in% pv_1 ~ "s2_pv_1",
    imageJ_coordinate %in% pv_2 ~ "s2_pv_2",
    imageJ_coordinate %in% pv_3 ~ "s2_pv_3",
    imageJ_coordinate %in% pv_4 ~ "s2_pv_4",
    imageJ_coordinate %in% pv_5 ~ "s2_pv_5",
    imageJ_coordinate %in% pv_6 ~ "s2_pv_6",
    imageJ_coordinate %in% cv_1 ~ "s2_cv_1",
    imageJ_coordinate %in% cv_2 ~ "s2_cv_2",
    imageJ_coordinate %in% cv_3 ~ "s2_cv_3",
    imageJ_coordinate %in% cv_4 ~ "s2_cv_4",
    imageJ_coordinate %in% cv_5 ~ "s2_cv_5",
    imageJ_coordinate %in% cv_6 ~ "s2_cv_6",
    imageJ_coordinate %in% cv_7 ~ "s2_cv_7",
    TRUE ~ NA_character_))


##save objects 
saveRDS(A2_1, file = "./data_files_generated/noMets_A2_1.rds")

##A2_2
#pv 
pv_1 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/A2_2/A2_2_pv_1.csv")
pv_2 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/A2_2/A2_2_pv_2.csv")
pv_3 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/A2_2/A2_2_pv_3.csv")
pv_4 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/A2_2/A2_2_pv_4.csv")
pv_5 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/A2_2/A2_2_pv_5.csv")

#cv 
cv_1 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/A2_2/A2_2_cv_1.csv")
cv_2 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/A2_2/A2_2_cv_2.csv")
cv_3 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/A2_2/A2_2_cv_3.csv")
cv_4 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/A2_2/A2_2_cv_4.csv")

A2_2_noMets_pv <- c(pv_1,pv_2,pv_3,pv_4,pv_5)
A2_2_noMets_cv <- c(cv_1,cv_2,cv_3,cv_4)

#add spatial info 
A2_2$spatial_feature <- NA
A2_2@meta.data <- A2_2@meta.data %>%
  mutate(spatial_feature = case_when(
    imageJ_coordinate %in% A2_2_noMets_pv ~ "pv",
    imageJ_coordinate %in% A2_2_noMets_cv ~ "cv",
    TRUE ~ NA_character_))

#remove NA 
Idents(A2_2) <- "spatial_feature"
A2_2 <- subset(A2_2, idents = c("cv","pv"))

#add spatial info id number
A2_2$spatial_feature_number <- NA
A2_2@meta.data <- A2_2@meta.data %>%
  mutate(spatial_feature_number = case_when(
    imageJ_coordinate %in% pv_1 ~ "s3_pv_1",
    imageJ_coordinate %in% pv_2 ~ "s3_pv_2",
    imageJ_coordinate %in% pv_3 ~ "s3_pv_3",
    imageJ_coordinate %in% pv_4 ~ "s3_pv_4",
    imageJ_coordinate %in% pv_5 ~ "s3_pv_5",
    imageJ_coordinate %in% cv_1 ~ "s3_cv_1",
    imageJ_coordinate %in% cv_2 ~ "s3_cv_2",
    imageJ_coordinate %in% cv_3 ~ "s3_cv_3",
    imageJ_coordinate %in% cv_4 ~ "s3_cv_4",
    TRUE ~ NA_character_))

##save objects 
saveRDS(A2_2, file = "./data_files_generated/noMets_A2_2.rds")

##B1_1
#pv surrounding mets 
pv_sm_1 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/B1_1/B1_1_pv_sm_1_1.csv")
pv_sm_2 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/B1_1/B1_1_pv_sm_1_2.csv")
pv_sm_3 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/B1_1/B1_1_pv_sm_3_1.csv")

#cv surrounding mets 
cv_sm_1 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/B1_1/B1_1_cv_sm_1_1.csv")
cv_sm_2 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/B1_1/B1_1_cv_sm_1_2.csv")
cv_sm_3 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/B1_1/B1_1_cv_sm_1_3.csv")
cv_sm_4 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/B1_1/B1_1_cv_sm_1_4.csv")
cv_sm_5 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/B1_1/B1_1_cv_sm_2_1.csv")
cv_sm_6 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/B1_1/B1_1_cv_sm_2_3_1.csv")
cv_sm_7 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/B1_1/B1_1_cv_sm_3_1.csv")

#metastasis
m1 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/B1_1/B1_1_mets1.csv")
m2 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/B1_1/B1_1_mets2.csv")
m3 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/B1_1/B1_1_mets3.csv")

#cv no metastasis 
cv_nm_1 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/B1_1/B1_1_cv_nm_1.csv")
cv_nm_2 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/B1_1/B1_1_cv_nm_2.csv")
cv_nm_3 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/B1_1/B1_1_cv_nm_3.csv")
cv_nm_4 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/B1_1/B1_1_cv_nm_4.csv")
cv_nm_5 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/B1_1/B1_1_cv_nm_5.csv")
cv_nm_6 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/B1_1/B1_1_cv_nm_6.csv")

#pv no metastasis 
pv_nm_1 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/B1_1/B1_1_pv_nm_1.csv")
pv_nm_2 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/B1_1/B1_1_pv_nm_2.csv")
pv_nm_3 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/B1_1/B1_1_pv_nm_3.csv")
pv_nm_4 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/B1_1/B1_1_pv_nm_4.csv")
pv_nm_5 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/B1_1/B1_1_pv_nm_5.csv")
pv_nm_6 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/B1_1/B1_1_pv_nm_6.csv")

B1_1_mets_pv_sm <- c(pv_sm_1, pv_sm_2,pv_sm_3)
B1_1_mets_cv_sm <- c(cv_sm_1,cv_sm_2,cv_sm_3,cv_sm_4,cv_sm_5,cv_sm_6,cv_sm_7)
B1_1_mets_pv_nm <- c(pv_nm_1,pv_nm_2,pv_nm_3,pv_nm_4,pv_nm_5,pv_nm_6)
B1_1_mets_cv_nm <- c(cv_nm_1,cv_nm_2,cv_nm_3,cv_nm_4,cv_nm_5,cv_nm_6)
B1_1_mets <- c(m1,m2,m3)

#add spatial info 
B1_1$spatial_feature <- NA
B1_1@meta.data <- B1_1@meta.data %>%
  mutate(spatial_feature = case_when(
    imageJ_coordinate %in% B1_1_mets_pv_sm ~ "pv_sm",
    imageJ_coordinate %in% B1_1_mets_cv_sm ~ "cv_sm",
    imageJ_coordinate %in% B1_1_mets_pv_nm ~ "pv_nm",
    imageJ_coordinate %in% B1_1_mets_cv_nm ~ "cv_nm",
    imageJ_coordinate %in% B1_1_mets ~ "mets",
    TRUE ~ NA_character_))

#remove NA 
Idents(B1_1) <- "spatial_feature"
B1_1 <- subset(B1_1, idents = c("cv_nm", "cv_sm" , "mets", "pv_nm" ,"pv_sm"))

#add spatial info id number
B1_1$spatial_feature_number <- NA
B1_1@meta.data <- B1_1@meta.data %>%
  mutate(spatial_feature_number = case_when(
    imageJ_coordinate %in% pv_sm_1 ~ "s4_pv_sm_1",
    imageJ_coordinate %in% pv_sm_2 ~ "s4_pv_sm_2",
    imageJ_coordinate %in% pv_sm_3 ~ "s4_pv_sm_3",
    imageJ_coordinate %in% cv_sm_1 ~ "s4_cv_sm_1",
    imageJ_coordinate %in% cv_sm_2 ~ "s4_cv_sm_2",
    imageJ_coordinate %in% cv_sm_3 ~ "s4_cv_sm_3",
    imageJ_coordinate %in% cv_sm_4 ~ "s4_cv_sm_4",
    imageJ_coordinate %in% cv_sm_5 ~ "s4_cv_sm_5",
    imageJ_coordinate %in% cv_sm_6 ~ "s4_cv_sm_6",
    imageJ_coordinate %in% cv_sm_7 ~ "s4_cv_sm_7",
    imageJ_coordinate %in% pv_nm_1 ~ "s4_pv_nm_1",
    imageJ_coordinate %in% pv_nm_2 ~ "s4_pv_nm_1",
    imageJ_coordinate %in% pv_nm_3 ~ "s4_pv_nm_3",
    imageJ_coordinate %in% pv_nm_4 ~ "s4_pv_nm_4",
    imageJ_coordinate %in% pv_nm_5 ~ "s4_pv_nm_5",
    imageJ_coordinate %in% pv_nm_6 ~ "s4_pv_nm_6",
    imageJ_coordinate %in% cv_nm_1 ~ "s4_cv_nm_1",
    imageJ_coordinate %in% cv_nm_2 ~ "s4_cv_nm_2",
    imageJ_coordinate %in% cv_nm_3 ~ "s4_cv_nm_3",
    imageJ_coordinate %in% cv_nm_4 ~ "s4_cv_nm_4",
    imageJ_coordinate %in% cv_nm_5 ~ "s4_cv_nm_5",
    imageJ_coordinate %in% cv_nm_6 ~ "s4_cv_nm_6",
    imageJ_coordinate %in% m1 ~ "s4_m1",
    imageJ_coordinate %in% m2 ~ "s4_m2",
    imageJ_coordinate %in% m3 ~ "s4_m3",
    TRUE ~ NA_character_))

##save objects 
saveRDS(B1_1, file = "./data_files_generated/mets_B1_1.rds")

##B1_2
#pv surrounding mets 
pv_sm_1 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/B1_2/B1_2_pv_sm_2_1.csv")
pv_sm_2 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/B1_2/B1_2_pv_sm_2_2.csv")
pv_sm_3 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/B1_2/B1_2_pv_sm_2_3_1.csv")

#cv surrounding mets 
cv_sm_1 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/B1_2/B1_2_cv_sm_1_1.csv")
cv_sm_2 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/B1_2/B1_2_cv_sm_1_2.csv")
cv_sm_3 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/B1_2/B1_2_cv_sm_3_1.csv")
cv_sm_4 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/B1_2/B1_2_cv_sm_3_2.csv")
cv_sm_5 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/B1_2/B1_2_cv_sm_3_3.csv")


#metastasis
m1 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/B1_2/B1_2_mets1.csv")
m2 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/B1_2/B1_2_mets2.csv")
m3 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/B1_2/B1_2_mets3.csv")

#cv no metastasis 
cv_nm_1 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/B1_2/B1_2_cv_nm_1.csv")
cv_nm_2 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/B1_2/B1_2_cv_nm_2.csv")

#pv no metastasis 
pv_nm_1 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/B1_2/B1_2_pv_nm_1.csv")

B1_2_mets_pv_sm <- c(pv_sm_1, pv_sm_2,pv_sm_3)
B1_2_mets_cv_sm <- c(cv_sm_1,cv_sm_2,cv_sm_3,cv_sm_4,cv_sm_5)
B1_2_mets_pv_nm <- c(pv_nm_1)
B1_2_mets_cv_nm <- c(cv_nm_1,cv_nm_2)
B1_2_mets <- c(m1,m2,m3)

#add spatial info 
B1_2$spatial_feature <- NA
B1_2@meta.data <- B1_2@meta.data %>%
  mutate(spatial_feature = case_when(
    imageJ_coordinate %in% B1_2_mets_pv_sm ~ "pv_sm",
    imageJ_coordinate %in% B1_2_mets_cv_sm ~ "cv_sm",
    imageJ_coordinate %in% B1_2_mets_pv_nm ~ "pv_nm",
    imageJ_coordinate %in% B1_2_mets_cv_nm ~ "cv_nm",
    imageJ_coordinate %in% B1_2_mets ~ "mets",
    TRUE ~ NA_character_))

#remove NA 
Idents(B1_2) <- "spatial_feature"
B1_2 <- subset(B1_2, idents = c("cv_nm", "cv_sm" , "mets", "pv_nm" ,"pv_sm"))

#add spatial info id number
B1_2$spatial_feature_number <- NA
B1_2@meta.data <- B1_2@meta.data %>%
  mutate(spatial_feature_number = case_when(
    imageJ_coordinate %in% pv_sm_1 ~ "s5_pv_sm_1",
    imageJ_coordinate %in% pv_sm_2 ~ "s5_pv_sm_2",
    imageJ_coordinate %in% pv_sm_3 ~ "s5_pv_sm_3",
    imageJ_coordinate %in% cv_sm_1 ~ "s5_cv_sm_1",
    imageJ_coordinate %in% cv_sm_2 ~ "s5_cv_sm_2",
    imageJ_coordinate %in% cv_sm_3 ~ "s5_cv_sm_3",
    imageJ_coordinate %in% cv_sm_4 ~ "s5_cv_sm_4",
    imageJ_coordinate %in% cv_sm_5 ~ "s5_cv_sm_5",
    imageJ_coordinate %in% pv_nm_1 ~ "s5_pv_nm_1",
    imageJ_coordinate %in% cv_nm_1 ~ "s5_cv_nm_1",
    imageJ_coordinate %in% cv_nm_2 ~ "s5_cv_nm_2",
    imageJ_coordinate %in% m1 ~ "s5_m1",
    imageJ_coordinate %in% m2 ~ "s5_m2",
    imageJ_coordinate %in% m3 ~ "s5_m3",
    TRUE ~ NA_character_))


##save objects 
saveRDS(B1_2, file = "./data_files_generated/mets_B1_2.rds")

##B2_1
#pv 
pv_1 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/B2_1/B2_1_pv_1.csv")
pv_2 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/B2_1/B2_1_pv_2.csv")
pv_3 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/B2_1/B2_1_pv_3.csv")
pv_4 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/B2_1/B2_1_pv_4.csv")
pv_5 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/B2_1/B2_1_pv_5.csv")
pv_6 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/B2_1/B2_1_pv_6.csv")
pv_7 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/B2_1/B2_1_pv_7.csv")
pv_8 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/B2_1/B2_1_pv_8.csv")
pv_9 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/B2_1/B2_1_pv_9.csv")
pv_10 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/B2_1/B2_1_pv_10.csv")
pv_11 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/B2_1/B2_1_pv_11.csv")
pv_12 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/B2_1/B2_1_pv_12.csv")
pv_13 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/B2_1/B2_1_pv_13.csv")
pv_14 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/B2_1/B2_1_pv_14.csv")
pv_15 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/B2_1/B2_1_pv_15.csv")
pv_16 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/B2_1/B2_1_pv_16.csv")
pv_17 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/B2_1/B2_1_pv_17.csv")
pv_18 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/B2_1/B2_1_pv_18.csv")
pv_19 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/B2_1/B2_1_pv_19.csv")
pv_20 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/B2_1/B2_1_pv_20.csv")
pv_21 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/B2_1/B2_1_pv_21.csv")
pv_22 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/B2_1/B2_1_pv_22.csv")
pv_23 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/B2_1/B2_1_pv_23.csv")
pv_24 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/B2_1/B2_1_pv_24.csv")
pv_25 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/B2_1/B2_1_pv_25.csv")
pv_26 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/B2_1/B2_1_pv_26.csv")
pv_27 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/B2_1/B2_1_pv_27.csv")
pv_28 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/B2_1/B2_1_pv_28.csv")

#cv 
cv_1 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/B2_1/B2_1_cv_1.csv")
cv_2 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/B2_1/B2_1_cv_2.csv")
cv_3 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/B2_1/B2_1_cv_3.csv")
cv_4 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/B2_1/B2_1_cv_4.csv")
cv_5 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/B2_1/B2_1_cv_5.csv")
cv_6 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/B2_1/B2_1_cv_6.csv")
cv_7 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/B2_1/B2_1_cv_7.csv")
cv_8 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/B2_1/B2_1_cv_8.csv")
cv_9 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/B2_1/B2_1_cv_9.csv")
cv_10 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/B2_1/B2_1_cv_10.csv")
cv_11 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/B2_1/B2_1_cv_11.csv")
cv_12 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/B2_1/B2_1_cv_12.csv")
cv_13 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/B2_1/B2_1_cv_13.csv")
cv_14 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/B2_1/B2_1_cv_14.csv")
cv_15 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/B2_1/B2_1_cv_15.csv")
cv_16 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/B2_1/B2_1_cv_16.csv")
cv_17 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/B2_1/B2_1_cv_17.csv")
cv_18 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/B2_1/B2_1_cv_18.csv")
cv_19 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/B2_1/B2_1_cv_19.csv")
cv_20 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/B2_1/B2_1_cv_20.csv")
cv_21 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/B2_1/B2_1_cv_21.csv")
cv_22 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/B2_1/B2_1_cv_22.csv")
cv_23 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/B2_1/B2_1_cv_23.csv")
cv_24 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/B2_1/B2_1_cv_24.csv")
cv_25 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/B2_1/B2_1_cv_25.csv")
cv_26 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/B2_1/B2_1_cv_26.csv")
cv_27 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/B2_1/B2_1_cv_27.csv")
cv_28 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/B2_1/B2_1_cv_28.csv")
cv_29 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/B2_1/B2_1_cv_29.csv")
cv_30 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/B2_1/B2_1_cv_30.csv")
cv_31 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/B2_1/B2_1_cv_31.csv")
cv_32 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/B2_1/B2_1_cv_32.csv")
cv_33 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/B2_1/B2_1_cv_33.csv")
cv_34 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/B2_1/B2_1_cv_34.csv")
cv_35 <- extracting_x_y_info_Resolve("./data/spatial_features_coordinates/B2_1/B2_1_cv_35.csv")

B2_1_noMets_pv <- c(pv_1,pv_2,pv_3,pv_4,pv_5,pv_6,pv_7,pv_8,pv_9,pv_10,
                    pv_11,pv_12,pv_13,pv_14,pv_15,pv_16,pv_17,pv_18,pv_19,pv_20,
                    pv_21,pv_22,pv_23,pv_24,pv_25,pv_26,pv_27,pv_28)
B2_1_noMets_cv <- c(cv_1,cv_2,cv_3,cv_4,cv_5,cv_6,cv_7,cv_8,cv_9,cv_10,
                    cv_11,cv_12,cv_13,cv_14,cv_15,cv_16,cv_17,cv_18,cv_19,cv_20,
                    cv_21,cv_22,cv_23,cv_24,cv_25,cv_26,cv_27,cv_28,cv_29,cv_30,
                    cv_31,cv_32,cv_33,cv_34,cv_35)

#add spatial info 
B2_1$spatial_feature <- NA
B2_1@meta.data <- B2_1@meta.data %>%
  mutate(spatial_feature = case_when(
    imageJ_coordinate %in% B2_1_noMets_pv ~ "pv",
    imageJ_coordinate %in% B2_1_noMets_cv ~ "cv",
    TRUE ~ NA_character_))

#remove NA 
Idents(B2_1) <- "spatial_feature"
B2_1 <- subset(B2_1, idents = c("cv","pv"))

#add spatial info id number
B2_1$spatial_feature_number <- NA
B2_1@meta.data <- B2_1@meta.data %>%
  mutate(spatial_feature_number = case_when(
    imageJ_coordinate %in% pv_1 ~ "s6_pv_1",
    imageJ_coordinate %in% pv_2 ~ "s6_pv_2",
    imageJ_coordinate %in% pv_3 ~ "s6_pv_3",
    imageJ_coordinate %in% pv_4 ~ "s6_pv_4",
    imageJ_coordinate %in% pv_5 ~ "s6_pv_5",
    imageJ_coordinate %in% pv_6 ~ "s6_pv_6",
    imageJ_coordinate %in% pv_7 ~ "s6_pv_7",
    imageJ_coordinate %in% pv_8 ~ "s6_pv_8",
    imageJ_coordinate %in% pv_9 ~ "s6_pv_9",
    imageJ_coordinate %in% pv_10 ~ "s6_pv_10",
    imageJ_coordinate %in% pv_11 ~ "s6_pv_11",
    imageJ_coordinate %in% pv_11 ~ "s6_pv_12",
    imageJ_coordinate %in% pv_13 ~ "s6_pv_13",
    imageJ_coordinate %in% pv_14 ~ "s6_pv_14",
    imageJ_coordinate %in% pv_15 ~ "s6_pv_15",
    imageJ_coordinate %in% pv_16 ~ "s6_pv_16",
    imageJ_coordinate %in% pv_17 ~ "s6_pv_17",
    imageJ_coordinate %in% pv_18 ~ "s6_pv_18",
    imageJ_coordinate %in% pv_19 ~ "s6_pv_19",
    imageJ_coordinate %in% pv_20 ~ "s6_pv_20",
    imageJ_coordinate %in% pv_21 ~ "s6_pv_21",
    imageJ_coordinate %in% pv_22 ~ "s6_pv_22",
    imageJ_coordinate %in% pv_23 ~ "s6_pv_23",
    imageJ_coordinate %in% pv_24 ~ "s6_pv_24",
    imageJ_coordinate %in% pv_25 ~ "s6_pv_25",
    imageJ_coordinate %in% pv_26 ~ "s6_pv_26",
    imageJ_coordinate %in% pv_27 ~ "s6_pv_27",
    imageJ_coordinate %in% pv_28 ~ "s6_pv_28",
    imageJ_coordinate %in% cv_1 ~ "s6_cv_1",
    imageJ_coordinate %in% cv_2 ~ "s6_cv_2",
    imageJ_coordinate %in% cv_3 ~ "s6_cv_3",
    imageJ_coordinate %in% cv_4 ~ "s6_cv_4",
    imageJ_coordinate %in% cv_5 ~ "s6_cv_5",
    imageJ_coordinate %in% cv_6 ~ "s6_cv_6",
    imageJ_coordinate %in% cv_7 ~ "s6_cv_7",
    imageJ_coordinate %in% cv_8 ~ "s6_cv_8",
    imageJ_coordinate %in% cv_9 ~ "s6_cv_9",
    imageJ_coordinate %in% cv_10 ~ "s6_cv_10",
    imageJ_coordinate %in% cv_11 ~ "s6_cv_11",
    imageJ_coordinate %in% cv_12 ~ "s6_cv_12",
    imageJ_coordinate %in% cv_13 ~ "s6_cv_13",
    imageJ_coordinate %in% cv_14 ~ "s6_cv_14",
    imageJ_coordinate %in% cv_15 ~ "s6_cv_15",
    imageJ_coordinate %in% cv_16 ~ "s6_cv_16",
    imageJ_coordinate %in% cv_17 ~ "s6_cv_17",
    imageJ_coordinate %in% cv_18 ~ "s6_cv_18",
    imageJ_coordinate %in% cv_19 ~ "s6_cv_19",
    imageJ_coordinate %in% cv_20 ~ "s6_cv_20",
    imageJ_coordinate %in% cv_21 ~ "s6_cv_21",
    imageJ_coordinate %in% cv_22 ~ "s6_cv_22",
    imageJ_coordinate %in% cv_23 ~ "s6_cv_23",
    imageJ_coordinate %in% cv_24 ~ "s6_cv_24",
    imageJ_coordinate %in% cv_25 ~ "s6_cv_25",
    imageJ_coordinate %in% cv_26 ~ "s6_cv_26",
    imageJ_coordinate %in% cv_27 ~ "s6_cv_27",
    imageJ_coordinate %in% cv_28 ~ "s6_cv_28",
    imageJ_coordinate %in% cv_29 ~ "s6_cv_29",
    imageJ_coordinate %in% cv_30 ~ "s6_cv_30",
    imageJ_coordinate %in% cv_31 ~ "s6_cv_31",
    imageJ_coordinate %in% cv_32 ~ "s6_cv_32",
    imageJ_coordinate %in% cv_33 ~ "s6_cv_33",
    imageJ_coordinate %in% cv_34 ~ "s6_cv_34",
    imageJ_coordinate %in% cv_35 ~ "s6_cv_35",
    TRUE ~ NA_character_))

##save objects 
saveRDS(B2_1, file = "./data_files_generated/noMets_B2_1.rds")

########## Combine all slides with added spatial area information ##########
###Read in objects 
B1_1 <- readRDS(file = "./data_files_generated/mets_B1_1.rds")
B1_2 <- readRDS(file = "./data_files_generated/mets_B1_2.rds")
A1_1 <- readRDS(file = "./data_files_generated/mets_A1_1.rds")

A2_1 <- readRDS(file = "./data_files_generated/noMets_A2_1.rds")
A2_2 <- readRDS(file = "./data_files_generated/noMets_A2_2.rds")
B2_1 <- readRDS(file = "./data_files_generated/noMets_B2_1.rds")

###add prefix to the sphere name 
mets_samples <- merge(B1_2,c(B1_1,A1_1))
noMets_samples <- merge(B2_1,c(A2_1,A2_2))

###save objects 
saveRDS(mets_samples, file = "./data_files_generated/mets_QC_imageJcoordinate.rds")
saveRDS(noMets_samples,file = "./data_files_generated/noMets_QC_imageJcoordinate.rds")






########## Part 4: Visium public data analysis ##########
#This part uses public data from livercellatlas.org (Guilliams et al. 2022, Cell) to test newly found zonation specific genes 
#in wild type mouse liver data and NAFLD mouse liver data

########## Prepare environment ##########
###Setting the working directory 
setwd("/mnt/khandler/R_projects/Sphere-sequencing/Visium_analysis/")

###Load packages and functions 
source("./functions_and_packages/1.Packages.R")

###load R object 
WT <- Read10X("./data/Visium_public/rawData_mouseStStVisium/countTable_mouseStStVisium",gene.column = 1)
Nafld <- Read10X("./data/Visium_public/rawData_mouseNafldVisium/countTable_mouseNafldVisium/",gene.column = 1)

#create Seruat objects
WT <- CreateSeuratObject(WT)
Nafld <- CreateSeuratObject(Nafld)

###load annotation files 
WT_anno <- read.csv("./data/Visium_public/annot_mouseStStVisium.csv")
Nafld_anno <- read.csv("./data/Visium_public/annot_mouseNafldVisium.csv")

########## wild type sample analysis ##########
#match spot with zonationGroup to get annotation in Seurat object 
WT_central_spots <- WT_anno[WT_anno$zonationGroup %in% "Central",]$spot
WT_mid_spots <- WT_anno[WT_anno$zonationGroup %in% "Mid",]$spot
WT_periportal_spots <- WT_anno[WT_anno$zonationGroup %in% "Periportal",]$spot
WT_portal_spots <- WT_anno[WT_anno$zonationGroup %in% "Portal",]$spot

###add meta data with the identity of the vein area 
WT$zonationGroup <- NA
WT@meta.data <- WT@meta.data %>%
  mutate(zonationGroup = case_when(
    rownames(WT@meta.data) %in% WT_central_spots ~ "Central",
    rownames(WT@meta.data) %in% WT_mid_spots ~ "Mid",
    rownames(WT@meta.data) %in% WT_periportal_spots ~ "Periportal",
    rownames(WT@meta.data) %in% WT_portal_spots ~ "Portal",
    TRUE ~ NA_character_))

#remove NA 
Idents(WT) <- "zonationGroup"
WT <- subset(WT, idents = c("Central","Mid","Periportal","Portal"))

###analyse zonated genes found in sphere-seq 
Idents(WT) <- "zonationGroup"
WT@active.ident <- factor(x = WT@active.ident, 
                          levels = c("Central","Mid","Periportal","Portal"))
p <- DotPlot(WT, features = c("Plpp1","Galnt15", 
                              "Vcam1","Itgb1","Ccl3","Ccr5"), dot.scale = 10) +
  theme(legend.title = element_text(size = 22), legend.text = element_text(size = 22)) + 
  theme(title = element_text(size = 25))+ theme(axis.text = element_text(size = 30)) +
  ggtitle("Public WT")  +  theme(axis.text.x = element_text(angle = 90)) 
p + ggsave("./figures/4/WT_public2.pdf",width = 12, height = 10)
p + ggsave("./figures/4/WT_public2.svg",width = 12, height = 10)

########## NAFLD sample analysis ##########
#match spot with zonationGroup to get annotation in Seurat object 
Nafld_central_spots <- Nafld_anno[Nafld_anno$zonationGroup %in% "Central",]$spot
Nafld_mid_spots <- Nafld_anno[Nafld_anno$zonationGroup %in% "Mid",]$spot
Nafld_periportal_spots <- Nafld_anno[Nafld_anno$zonationGroup %in% "Periportal",]$spot
Nafld_portal_spots <- Nafld_anno[Nafld_anno$zonationGroup %in% "Portal",]$spot

###add meta data with the identity of the vein area 
Nafld$zonationGroup <- NA
Nafld@meta.data <- Nafld@meta.data %>%
  mutate(zonationGroup = case_when(
    rownames(Nafld@meta.data) %in% Nafld_central_spots ~ "Central",
    rownames(Nafld@meta.data) %in% Nafld_mid_spots ~ "Mid",
    rownames(Nafld@meta.data) %in% Nafld_periportal_spots ~ "Periportal",
    rownames(Nafld@meta.data) %in% Nafld_portal_spots ~ "Portal",
    TRUE ~ NA_character_))

#remove NA 
Idents(Nafld) <- "zonationGroup"
Nafld <- subset(Nafld, idents = c("Central","Mid","Periportal","Portal"))

###analyse zonated genes found in sphere-seq 
p <- DotPlot(Nafld, features = c("Plpp1","Galnt15", 
                              "Vcam1","Itgb1","Ccl3","Ccr5"), dot.scale = 10) +
  theme(legend.title = element_text(size = 22), legend.text = element_text(size = 22)) + 
  theme(title = element_text(size = 25))+ theme(axis.text = element_text(size = 30)) +
  ggtitle("Public Nafld")  +  theme(axis.text.x = element_text(angle = 90)) 
p + ggsave("./figures/4/Nafld_public2.pdf",width = 12, height = 10)
p + ggsave("./figures/4/Nafld_public2.svg",width = 12, height = 10)

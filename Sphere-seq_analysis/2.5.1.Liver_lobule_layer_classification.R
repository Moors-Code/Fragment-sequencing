########## Part 2.5.1: Liver lobule layer classification ##########
#This part allocates spheres into lobule layers between central and portal vein depending on landmark gene expression in LECs 

########## Prepare environment ##########
###Setting the working directory 
setwd("/mnt/khandler/R_projects/Sphere-sequencing/Sphere-seq_analysis/")

###Load packages and functions 
source("./functions_and_packages/1.Packages.R")
source("./functions_and_packages/5.Functions_lobule_layers.R")

###load R object
liverSpS5C <- readRDS("./data_files_generated/LiverMerged_afterBC_anno_BS_5cells.Rda")

########## Add zonation coordinate (ZC), remove spheres without LECs and group into 8 lobule areas ##########
#lobule areas: L1-L3, L4, L5, L6, L7, L8-L10, not many spheres for most central and most portal therefore group into L1-L3 and L8-L10
add_zonation_coefficient(
  "./data_files_generated/LiverMerged_afterBC_anno_BS_5cells.Rda",
  "./data_files_generated/LiverMerged_afterBC_anno_BS_5cells_zC.Rda")
remove_spheres_noLECs(
  "./data_files_generated/LiverMerged_afterBC_anno_BS_5cells_zC.Rda",
  "./data_files_generated/LiverMerged_afterBC_anno_BS_5cells_zC.Rda")
add_lobule_zones(
  "./data_files_generated/LiverMerged_afterBC_anno_BS_5cells_zC.Rda",
  "./data_files_generated/LiverMerged_afterBC_anno_BS_5cells_zC_lobules.Rda")


########## Check ZC algorithm with marker genes in hepatocytes ##########
liverSpS5C <- readRDS("./data_files_generated/LiverMerged_afterBC_anno_BS_5cells_zC_lobules.Rda")

###subset hepatocytes 
Idents(liverSpS5C) <- "annotation"
hep <- subset(liverSpS5C, idents = "Hepatocytes")

###plot CV (Cyp2e1,"Cyp1a2) and PV (Cyp2f2,Alb) landmark genes in Dotplot 
Idents(hep) <- "vein"
p <- DotPlot(hep, features = c("Cyp2e1","Cyp1a2","Cyp2f2","Alb"), dot.scale = 50) + 
  theme(legend.title = element_text(size = 22), legend.text = element_text(size = 22)) + 
  theme(title = element_text(size = 25))+ theme(axis.text = element_text(size = 30)) +
  ggtitle("Landmark genes in hepatocytes of zonated spheres")  +  theme(axis.text.x = element_text(angle = 90)) 
p + ggsave("./figures/2.5.1/Hep_LM_dotplot_vein.pdf",width = 12, height = 10)
p + ggsave("./figures/2.5.1/Hep_LM_dotplot_vein.svg",width = 12, height = 10)



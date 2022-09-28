########## Part 1.2: Pre-processing - WTA integration with MULTI-seq barcodes ##########
#This part integrates Cell IDs allocated to MULTI-seq barcodes from 1.1. to Cell IDs of Seurat object of whole transcriptome (WTA) analysis

########## Prepare environment ##########
###Setting the working directory 
setwd("/mnt/khandler/R_projects/Sphere-sequencing/Sphere-seq_analysis/")

###Load packages and functions 
source("./functions_and_packages/1.Packages.R")
source("./functions_and_packages/2.Functions_preprocessing.R")

###load human and mouse ensemble ID/Gene ID dataframes 
##mouse gene/ensemble IDs 
ensembl<-useEnsembl(biomart="ensembl")
list<-listDatasets(ensembl)
mart <- useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl",version = 95)
attributes<-listAttributes(mart)
gene_ids <- getBM(attributes = c("ensembl_gene_id_version","external_gene_name"), mart = mart)
write.csv(gene_ids, "./data/gene_ids_ensemble_ids_mouse.csv")

##human gene/ensemble IDs 
mart1 <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl",version = 95)
attributes1<-listAttributes(mart1)
gene_ids1 <- getBM(attributes = c("ensembl_gene_id_version","external_gene_name"), mart = mart1)
write.csv(gene_ids1, "./data/gene_ids_ensemble_ids_human.csv")

gene_ids_human <- read.csv( "./data/gene_ids_ensemble_ids_human.csv")
gene_ids_mouse <- read.csv( "./data/gene_ids_ensemble_ids_mouse.csv")

########## Liver samples ##########
###Sample 1
MULTIseq_WTA_integration_zUMIoutputBD_mouse(
  "./data_files_generated/classified_table_Liver_sample1.Rda",
  "./data/MetsM1.dgecounts.rds",
  "M1",
  "./data_files_generated/SpS_LiverMets_sample1.Rda")

###Sample 2
MULTIseq_WTA_integration_zUMIoutputBD_mouse(
  "./data_files_generated/classified_table_Liver_sample2.Rda",
  "./data/MetsM3.dgecounts.rds",
  "M3",
  "./data_files_generated/SpS_LiverMets_sample2.Rda")

###Sample 3
MULTIseq_WTA_integration_zUMIoutputBD_mouse(
  "./data_files_generated/classified_table_Liver_sample3.Rda",
  "./data/MetsM4.dgecounts.rds",
  "M4",
  "./data_files_generated/SpS_LiverMets_sample3.Rda")

###Sample 4
MULTIseq_WTA_integration_zUMIoutputBD_mouse(
  "./data_files_generated/classified_table_Liver_sample4.Rda",
  "./data/Mets4M1.dgecounts.rds",
  "4M1",
  "./data_files_generated/SpS_LiverMets_sample4.Rda")

###Sample 5
MULTIseq_WTA_integration_zUMIoutputBD_mouse(
  "./data_files_generated/classified_table_Liver_sample5.Rda",
  "./data/Mets4M2.dgecounts.rds",
  "4M2",
  "./data_files_generated/SpS_LiverMets_sample5.Rda")

###Sample 6
MULTIseq_WTA_integration_zUMIoutputBD_mouse(
  "./data_files_generated/classified_table_Liver_sample6.Rda",
  "./data/MetsM5.dgecounts.rds",
  "M5",
  "./data_files_generated/SpS_LiverMets_sample6.Rda")

###Sample 7
MULTIseq_WTA_integration_zUMIoutputBD_mouse(
  "./data_files_generated/classified_table_Liver_sample7.Rda",
  "./data/Mets6M1.dgecounts.rds",
  "6M1",
  "./data_files_generated/SpS_LiverMets_sample7.Rda")

###Sample 8
MULTIseq_WTA_integration_zUMIoutputBD_mouse(
  "./data_files_generated/classified_table_Liver_sample8.Rda",
  "./data/Mets6M2.dgecounts.rds",
  "6M2",
  "./data_files_generated/SpS_LiverMets_sample8.Rda")

###Sample 9
MULTIseq_WTA_integration_zUMIoutputBD_mouse(
  "./data_files_generated/classified_table_Liver_sample9.Rda",
  "./data/Mets6M3.dgecounts.rds",
  "6M3",
  "./data_files_generated/SpS_LiverMets_sample9.Rda")

###Sample 10
MULTIseq_WTA_integration_zUMIoutputBD_mouse(
  "./data_files_generated/classified_table_Liver_sample10.Rda",
  "./data/WT.dgecounts.rds",
  "WTCre",
  "./data_files_generated/SpS_LiverWT_sample10.Rda")

########## CRC Organoid mixing species sample ##########
#Object with Gene IDs 
MULTIseq_WTA_integration_zUMIoutputBD_human_mouse(
  "./data_files_generated/classified_table_Organoid_mixing.Rda",
  "./data/OrgMixing.dgecounts.rds",
  "OrgMix",
  "./data_files_generated/SpS_Organoid_mixing.Rda")

#Object with ensemble IDs (this is needed to differentiate human and mouse UMI counts based on ensemble IDs in Part 3) 
MULTIseq_WTA_integration_zUMIoutputBD_ensembleIDs(
  "./data_files_generated/classified_table_Organoid_mixing.Rda",
  "./data/OrgMixing.dgecounts.rds",
  "OrgMix",
  "./data_files_generated/SpS_Organoid_mixing_ensembleIDs.Rda"
)

########## Spleen samples #####
###Sample 1
spleen1 <- Read10X(data.dir = "./data/Spleen_sample1_filtered_featrue_bc_matrix")
spleen1<- CreateSeuratObject(counts = spleen1, project = "S1", min.cells = 3, min.features = 200)

MULTIseq_WTA_integration_10X(
  "./data_files_generated/classified_table_Spleen_sample1.Rda",
  spleenS1,
  "./data_files_generated/SpS_Spleen_Sample1.Rda"
)

###Sample 2
spleen2 <- Read10X(data.dir = "./data/Spleen_sample2_filtered_feature_bc_matrix")
spleen2<- CreateSeuratObject(counts = spleen2, project = "S2", min.cells = 3, min.features = 200)

MULTIseq_WTA_integration_10X(
  "./data_files_generated/classified_table_Spleen_sample2.Rda",
  spleenS2,
  "./data_files_generated/SpS_Spleen_Sample2.Rda"
)

######### Crohn's biopsy samples ##########
###P387
MULTIseq_WTA_integration_zUMIoutputBD_human(
  "./data_files_generated/classified_table_Crohn_P387.Rda",
  "./data/CrohnP387.dgecounts.rds",
  "P387",
  "./data_files_generated/SpS_Crohn_P387.Rda")

###P393
MULTIseq_WTA_integration_zUMIoutputBD_human(
  "./data_files_generated/classified_table_Crohn_P393.Rda",
  "./data/CrohnP393.dgecounts.rds",
  "P393",
  "./data_files_generated/SpS_Crohn_P393.Rda")

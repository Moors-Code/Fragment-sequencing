###Function to generate MULTI-seq barcode UMI count matrix of BD Rhapsody data 
MULTIseq_UMI_count_matrix_generation_BD <- function (
  MULTIseq_BC_reference_directory,
  zUMIs_output_kept_barcodes_directory,
  MULTIseq_R1_fastq_file_directory,
  MULTIseq_R2_fastq_file_directory
) {
  #load barcode reference file 
  MULTI_seq_barcodes_and_indices <- readxl::read_excel(MULTIseq_BC_reference_directory)
  #extract only barcode sequences 
  bar.ref <- MULTI_seq_barcodes_and_indices$`Barcode Sequence`
  
  #load text file of kept barcodes from zUMI outputs 
  kept_barcodes <- read.csv(zUMIs_output_kept_barcodes_directory)
  
  #only extract column with barcode sequences 
  cell.id.vec <- as.character(kept_barcodes$XC)
  
  #add constant regions to the barcode of zUMI output to match MULTI-seq FASTQ files, 
  #in the end it should look like that BC1_CR1_BC2_CR2_BC3
  constant_region1 <- "ACTGGCCTGCGA"
  constant_region2 <- "GGTAGCGGTGACA"
  zUMI_output_BC_list <- list(cell.id.vec)
  final_BC_list <- c()
  for(BC in zUMI_output_BC_list) {
    BC1 = str_sub(BC, 1,9)
    BC2 = str_sub(BC, 10,18)
    BC3 = str_sub(BC, 19,27)
    final_BC = paste(BC1, constant_region1, BC2, constant_region2, BC3, sep = "")
    final_BC_list <- c(final_BC_list, final_BC) 
  }
  cell.id.vec <- as.character(final_BC_list)
  
  #extract required information from FASTQ files of MULTI-seq library 
  #Cell ID position 1:52, UMI position 53:60, MULTI-seq ID position 1:8 (tag)
  readTable <- MULTIseq.preProcess(R1 = MULTIseq_R1_fastq_file_directory, 
                                   R2 = MULTIseq_R2_fastq_file_directory, 
                                   cellIDs = cell.id.vec, cell = c(1,52), umi = c(53,60), tag = c(1,8))

  #align cell ID information with reads and MULTI-seq reference to generate MULTI-seq barcode UMI matrix per cell 
  bar.table.start <- MULTIseq.align(readTable, cell.id.vec, bar.ref)
  return(bar.table.start)
}

###Function to generate MULTI-seq barcode UMI count matrix of Chromium 10X data 
MULTIseq_bar_table_generation10X <- function (
  MULTIseq_BC_reference_directory,
  feature_bc_matrix_10X_directory,
  MULTIseq_R1_fastq_file_directory,
  MULTIseq_R2_fastq_file_directory
) {
  #load barcode reference file 
  MULTI_seq_barcodes_and_indices <- readxl::read_excel(MULTIseq_BC_reference_directory)
  #extract barcode sequences 
  bar.ref <- MULTI_seq_barcodes_and_indices$`Barcode Sequence`
  #load 10X filtered feature bc matrix 
  data <- Read10X(data.dir = feature_bc_matrix_10X_directory)
  #extract cell IDs 
  cell.id.vec <- data@Dimnames[[2]]
  cell.id.vec <- str_sub(cell.id.vec, 1, str_length(cell.id.vec)-2)
  #extract required information from FASTQ files from MULTI-seq library 
  #Cell ID position 1:16, UMI position 17:26, MULTI-seq ID position 1:8 (tag)
  readTable <- MULTIseq.preProcess(R1 = MULTIseq_R1_fastq_file_directory, 
                                   R2 =  MULTIseq_R2_fastq_file_directory, 
                                   cellIDs = cell.id.vec, cell = c(1,16), umi=c(17,26), tag = c(1,8))
  #align cell ID information with reads and MULTI-seq reference to generate MULTI-seq barcode UMI matrix per cell 
  bar.table.start <- MULTIseq.align(readTable, cell.id.vec, bar.ref)
  return(bar.table.start)
}

###Function to change ensemble IDs of zUMI outputs to gene name IDs - Mouse 
Annotation_mouse <- function(zUMI_output) {
  ##load ensemble gene ID data from mmusculus 
  #ensembl<-useEnsembl(biomart="ensembl")
  #list<-listDatasets(ensembl)
  #mart <- useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl",version = 95)
  #attributes<-listAttributes(mart)
  #add GFP and mCherry, it was added to the STAR index with the ID "CellTagUTR" and "sLPmCherry"
  #rename to "GFP" and "mCherry" 
  gfp_id <- data.frame(ensembl_gene_id_version = "CellTagUTR",external_gene_name= "GFP")
  sLP_mCherry_id <- data.frame(ensembl_gene_id_version = "sLPmCherry",external_gene_name= "mCherry")
  #gene_ids <- getBM(attributes = c("ensembl_gene_id_version","external_gene_name"), mart = mart)
  gene_ids <- gene_ids_mouse 
  gene_ids$X <- NULL
  #add GFP and mCherry to Gene ID dataframe 
  all_ids <- rbind(gene_ids,gfp_id)
  all_ids2 <- rbind(all_ids, sLP_mCherry_id)
  #use all exon UMI counts from zUMI output 
  dataframe<-as.data.frame(as.matrix(zUMI_output$umicount$exon$all))%>%
    as.matrix(.)
  #add column name "ensemble_gene_id_version" to the column of Ensemble IDs in zUMI output count matrix 
  dataframe<-mutate(as.data.frame(dataframe),ensembl_gene_id_version=rownames(dataframe))
  #compare ensemble IDs of zUMI count matrix with gene ID data frame to match ensemble IDs and to rename them to external gene names 
  join<-dataframe%>%
    left_join(dplyr::select(all_ids2,1:2))
  #remove duplicates 
  join<-join[!duplicated(join$external_gene_name),]
  join[is.na(join)]<-0 #make all empty value to zero
  #add external gene names as rownames (instead of ensemble IDs before)
  rownames(join)<-join$external_gene_name
  #remove unnessesary columns 
  join<-dplyr::select(join,-ensembl_gene_id_version,-external_gene_name)
  return(join)
}

###Function to change ensemble IDs of zUMI outputs to gene name IDs - Human 
Annotation_human <- function(zUMI_output) {
  #ensembl<-useEnsembl(biomart="ensembl")
  #list<-listDatasets(ensembl)
  #mart <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl",version = 95)
  #attributes<-listAttributes(mart)
  #gene_ids <- getBM(attributes = c("ensembl_gene_id_version","external_gene_name"), mart = mart)
  gene_ids <- gene_ids_human 
  gene_ids$X <- NULL
  #and also remove the "." prefix of the ENSMUS numbers 
  gene_ids2 <- str_replace(gene_ids$ensembl_gene_id_version,pattern = ".[0-9]+$",replacement = "")
  gene_ids$ensembl_gene_id_version <- gene_ids2
  dataframe<-as.data.frame(as.matrix(zUMI_output$umicount$exon$all))%>%
    as.matrix(.)
  dataframe<-mutate(as.data.frame(dataframe),ensembl_gene_id_version=rownames(dataframe))
  join<-dataframe%>%
    full_join(dplyr::select(gene_ids,1:2))
  length(unique(join$external_gene_name))
  join<-join[!duplicated(join$external_gene_name),]
  join[is.na(join)]<-0 
  rownames(join)<-join$external_gene_name
  join<-dplyr::select(join,-ensembl_gene_id_version,-external_gene_name)
  return(join)
}

###Function to change ensemble IDs of zUMI outputs to gene name IDs - Human/Mouse mixed  
Annotation_human_mouse <- function(zUMI_output) {
  #ensembl<-useEnsembl(biomart="ensembl")
  #list<-listDatasets(ensembl)
  #human
  #mart1 <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl",version = 95)
  #attributes1<-listAttributes(mart1)
  gene_ids1 <- gene_ids_human 
  gene_ids1$X <- NULL
  #gene_ids1 <- getBM(attributes = c("ensembl_gene_id_version","external_gene_name"), mart = mart1)
  gene_ids3 <- str_replace(gene_ids1$ensembl_gene_id_version,pattern = ".[0-9]+$",replacement = "")
  gene_ids1$ensembl_gene_id_version <- gene_ids3
  #mouse
  #mart2 <- useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl",version = 95)
  #attributes2<-listAttributes(mart2)
  #gene_ids2 <- getBM(attributes = c("ensembl_gene_id_version","external_gene_name"), mart = mart2)
  gene_ids2 <- gene_ids_mouse 
  gene_ids2$X <- NULL
  gene_ids4 <- str_replace(gene_ids2$ensembl_gene_id_version,pattern = ".[0-9]+$",replacement = "")
  gene_ids2$ensembl_gene_id_version <- gene_ids4
  all_ids <- rbind(gene_ids1,gene_ids2)
  dataframe<-as.data.frame(as.matrix(zUMI_output$umicount$exon$all))%>%
    as.matrix(.)
  dataframe<-mutate(as.data.frame(dataframe),ensembl_gene_id_version=rownames(dataframe))
  genes <- str_replace(dataframe$ensembl_gene_id_version,pattern = "mm10___|GRCh38_",replacement = "")
  dataframe$ensembl_gene_id_version <- genes
  join<-dataframe%>%
    full_join(dplyr::select(all_ids,1:2))
  join<-join[!duplicated(join$external_gene_name),]
  join[is.na(join)]<-0 
  rownames(join)<-join$external_gene_name
  join<-dplyr::select(join,-ensembl_gene_id_version,-external_gene_name)
  return(join)
}

###Function to integrate MULTIseq BC with WTA of zUMI output of BD Rhapsody data - Mouse
MULTIseq_WTA_integration_zUMIoutputBD_mouse <- function(
  classified_table_directory,
  dgecounts_zUMI_output_directory,
  sample_id,
  Seurat_object_save_directory
){
  classified_table <- readRDS(file = classified_table_directory)
  dge <- readRDS(dgecounts_zUMI_output_directory)
  dgeA <- Annotation_mouse(dge)
  dgeS <- CreateSeuratObject(dgeA, sample_id, min.cells = 3, min.features = 200)
  hash_table <- classified_table[,c("Row.names", "Barcode")]
  rownames(hash_table) <- hash_table$Row.names
  
  #remove the constant regions from row names of hash table
  BC_list <- list(rownames(hash_table))
  new_BC_list <- c()
  for(BC in BC_list) {
    BC1 = str_sub(BC, 1,9)
    BC2 = str_sub(BC, 22,30)
    BC3 = str_sub(BC, 44,52)
    final_BC = paste(BC1, BC2, BC3, sep = "")
    new_BC_list <- c(new_BC_list, final_BC) 
  }

  rownames(hash_table) <- new_BC_list
  hash_table$Row.names <- new_BC_list
  #extract cell IDs that are in classified table but not in zUMI output 
  cells_in_hash_table <- setdiff(hash_table$Row.names,rownames(dgeS@meta.data))
  extra.cells <- data.frame(Row.names = cells_in_hash_table)
  #remove extra cells so you have the same cell IDs 
  hash_table <- hash_table[! hash_table$Row.names %in% extra.cells$Row.names,]
  #extract common barcodes from both 
  common_barcodes <- intersect(rownames(dgeS@meta.data), hash_table$Row.names)
  #match cell IDs from classified table and zUMI output and add sphere BC (MUTLI-seq BC)
  dgeS@meta.data <- cbind(dgeS@meta.data, sphere = hash_table[rownames(dgeS@meta.data), "Barcode"] )
  saveRDS(dgeS,file = Seurat_object_save_directory)
}

###Function to integrate MULTIseq BC with WTA of zUMI output of BD Rhapsody data - Human 
MULTIseq_WTA_integration_zUMIoutputBD_human <- function(
  classified_table_directory,
  dgecounts_zUMI_output_directory,
  sample_id,
  Seurat_object_save_directory
){
  classified_table <- readRDS(file = classified_table_directory)
  dge <- readRDS(dgecounts_zUMI_output_directory)
  dgeA <- Annotation_human(dge)
  dgeS <- CreateSeuratObject(dgeA, sample_id, min.cells = 3, min.features = 200)
  hash_table <- classified_table[,c("Row.names", "Barcode")]
  rownames(hash_table) <- hash_table$Row.names
  
  #remove the constant regions from row names of hash table
  BC_list <- list(rownames(hash_table))
  new_BC_list <- c()
  for(BC in BC_list) {
    BC1 = str_sub(BC, 1,9)
    BC2 = str_sub(BC, 22,30)
    BC3 = str_sub(BC, 44,52)
    final_BC = paste(BC1, BC2, BC3, sep = "")
    new_BC_list <- c(new_BC_list, final_BC) 
  }
  
  rownames(hash_table) <- new_BC_list
  hash_table$Row.names <- new_BC_list
  cells_in_hash_table <- setdiff(hash_table$Row.names,rownames(dgeS@meta.data))
  extra.cells <- data.frame(Row.names = cells_in_hash_table)
  hash_table <- hash_table[! hash_table$Row.names %in% extra.cells$Row.names,]
  common_barcodes <- intersect(rownames(dgeS@meta.data), hash_table$Row.names)
  dgeS@meta.data <- cbind(dgeS@meta.data, sphere = hash_table[rownames(dgeS@meta.data), "Barcode"] )
  saveRDS(dgeS,file = Seurat_object_save_directory)
}

###Function to integrate MULTIseq BC with WTA of zUMI output of BD Rhapsody data - Human/Mouse mixed   
MULTIseq_WTA_integration_zUMIoutputBD_human_mouse <- function(
  classified_table_directory,
  dgecounts_zUMI_output_directory,
  sample_id,
  Seurat_object_save_directory
){
  classified_table <- readRDS(file = classified_table_directory)
  dge <- readRDS(dgecounts_zUMI_output_directory)
  dgeA <- Annotation_human_mouse(dge)
  dgeS <- CreateSeuratObject(dgeA, sample_id, min.cells = 3, min.features = 200)
  hash_table <- classified_table[,c("Row.names", "Barcode")]
  rownames(hash_table) <- hash_table$Row.names
  
  #remove the constant regions from row names of hash table
  BC_list <- list(rownames(hash_table))
  new_BC_list <- c()
  for(BC in BC_list) {
    BC1 = str_sub(BC, 1,9)
    BC2 = str_sub(BC, 22,30)
    BC3 = str_sub(BC, 44,52)
    final_BC = paste(BC1, BC2, BC3, sep = "")
    new_BC_list <- c(new_BC_list, final_BC) 
  }
  
  rownames(hash_table) <- new_BC_list
  hash_table$Row.names <- new_BC_list
  cells_in_hash_table <- setdiff(hash_table$Row.names,rownames(dgeS@meta.data))
  extra.cells <- data.frame(Row.names = cells_in_hash_table)
  hash_table <- hash_table[! hash_table$Row.names %in% extra.cells$Row.names,]
  common_barcodes <- intersect(rownames(dgeS@meta.data), hash_table$Row.names)
  dgeS@meta.data <- cbind(dgeS@meta.data, sphere = hash_table[rownames(dgeS@meta.data), "Barcode"] )
  saveRDS(dgeS,file = Seurat_object_save_directory)
}

###Function to integrate MULTIseq BC with WTA of Cell Ranger output of Chromium 10X 
MULTIseq_WTA_integration_10X <- function(
  classified_table_directory,
  Seurat_object,
  Seurat_object_save_directory
){
  classified_table <- readRDS(file = classified_table_directory)
  hash_table <- classified_table[,c("Row.names", "Barcode")]
  rownames(hash_table) <- hash_table$Row.names
  cells_in_hash_table <- setdiff(rownames(hash_table), str_sub(rownames(Seurat_object@meta.data), 1, str_length(rownames(Seurat_object@meta.data))-2))
  extra.cells <- data.frame(Row.names = cells_in_hash_table)
  hash_table <- hash_table[! hash_table$Row.names %in% extra.cells$Row.names,]
  common_barcodes <- intersect(str_sub(rownames(Seurat_object@meta.data), 1, str_length(rownames(Seurat_object@meta.data))-2), rownames(hash_table))
  Seurat_object@meta.data <- cbind(Seurat_object@meta.data, sphere = hash_table[str_sub(rownames(Seurat_object@meta.data), 1, str_length(rownames(Seurat_object@meta.data))-2), "Barcode"] )
  saveRDS(Seurat_object, file = Seurat_object_save_directory)
}


###Function to integrate MULTIseq BC with WTA of zUMI output of BD Rhapsody data while keeping ensemble IDs 
MULTIseq_WTA_integration_zUMIoutputBD_ensembleIDs <- function(
  classified_table_directory,
  dgecounts_zUMI_output_directory,
  sample_id,
  Seurat_object_save_directory
){
  classified_table <- readRDS(file = classified_table_directory)
  dge <- readRDS(dgecounts_zUMI_output_directory)
  dgeS <- CreateSeuratObject(dge$umicount$exon$all, sample_id, min.cells = 3, min.features = 200)
  hash_table <- classified_table[,c("Row.names", "Barcode")]
  rownames(hash_table) <- hash_table$Row.names
  
  #remove the constant regions from row names of hash table
  BC_list <- list(rownames(hash_table))
  new_BC_list <- c()
  for(BC in BC_list) {
    BC1 = str_sub(BC, 1,9)
    BC2 = str_sub(BC, 22,30)
    BC3 = str_sub(BC, 44,52)
    final_BC = paste(BC1, BC2, BC3, sep = "")
    new_BC_list <- c(new_BC_list, final_BC) 
  }
  
  rownames(hash_table) <- new_BC_list
  hash_table$Row.names <- new_BC_list
  cells_in_hash_table <- setdiff(hash_table$Row.names,rownames(dgeS@meta.data))
  extra.cells <- data.frame(Row.names = cells_in_hash_table)
  hash_table <- hash_table[! hash_table$Row.names %in% extra.cells$Row.names,]
  common_barcodes <- intersect(rownames(dgeS@meta.data), hash_table$Row.names)
  dgeS@meta.data <- cbind(dgeS@meta.data, sphere = hash_table[rownames(dgeS@meta.data), "Barcode"] )
  saveRDS(dgeS,file = Seurat_object_save_directory)
}


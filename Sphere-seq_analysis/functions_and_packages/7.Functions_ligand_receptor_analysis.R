###Function for generation of input files for CellPhoneDB 
#first gene symbols have to be converted from mouse to human because CellPhoneDB database only contains human L-R interactions   
#Conversion of mouse to human genes adapted partially from: https://github.com/CostaLab/CrossTalkeR/blob/master/CellPhoneDB%20Tutorial.md
#and: https://www.cellphonedb.org/faq-and-troubleshooting
#outputs are two text tiles: gene counts and cell annotations (meta)
Input_files_CellPhoneDB_generation <- function(
  seurat_object,
  annotation_column,
  sample_name,
  ouput_file_path
){
  #generating counts file 
  # take raw data and normalize it
  count_raw_meta <- GetAssayData(object = seurat_object, slot = "counts")[,colnames(x = seurat_object)]
  count_norm_meta <- apply(count_raw_meta, 2, function(x) (x/sum(x))*10000)
  genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = rownames(count_norm_meta) , mart = mouse, attributesL = c("hgnc_symbol","hgnc_id",'ensembl_gene_id'), martL = human, uniqueRows=T)
  print(head(genesV2))
  matrixA <- count_norm_meta[match(genesV2$MGI.symbol,rownames(count_norm_meta),nomatch=F),]
  matrixB <- matrixA
  matrixB$gene <- genesV2$Gene.stable.ID
  rownames(matrixA) <- matrixB$gene
  #save count matrix as text file 
  write.table(matrixA, paste0(ouput_file_path,sample_name,"_count.txt"), sep='\t', quote=F, row.names = T)
  # generating meta file based on cell type annotation of Seurat object 
  meta_data_meta <- cbind(rownames(seurat_object@meta.data), seurat_object@meta.data[,annotation_column, drop=F])  
  #save meta file as text file 
  write.table(meta_data_meta, paste0(ouput_file_path,sample_name,"_meta.txt"), sep='\t', quote=F, row.names=F)
}

###Function to generate score matrix of interacting pair of interest per mouse 
Score_matrix_generation_interacting_pair_oi_per_mouse <- function(
  cpdb_outputs_file_path,
  samples,
  interacting_cells_oi_pair,
  ouput_file_path,
  sample_name
){
  #combine significant_means files from these spheres by filling in 0s for no score 
  my_files <- paste0(cpdb_outputs_file_path,
                     samples,"/significant_means.txt")
  
  my_data <- list()
  for (i in seq_along(my_files)) {
    my_data[[i]] <- read.table(file = my_files[i],header = TRUE,sep = "\t")
    my_data[[i]][is.na(my_data[[i]])] <- 0
    my_data[[i]] <- my_data[[i]][,which(names(my_data[[i]]) %in% c("interacting_pair",interacting_cells_oi_pair))]
  }
  
  combined_df <- Reduce(function(x, y) merge(x, y, by=1),my_data)
  
  combined_df$sum <-rowSums(combined_df[sapply(combined_df, is.numeric)], na.rm = TRUE)
  
  #remove rows with 0
  combined_df <- combined_df[combined_df$sum !=0,]
  combined_df$sum <- NULL
  
  #add barcode ID information as colnames
  colnames(combined_df) <- c("interacting_pair",samples)
  
  #sort interactions decreasing by norm_sum 
  write.csv(combined_df, file = paste0(ouput_file_path,sample_name, "_per_mouse_",interacting_cells_oi_pair,".csv"))
}

###Function for DE analysis of L-R interactions scores between two conditions
DE_CellPhoneDB <- function(
  csv_file_path_cond1,
  csv_file_path_cond2,
  samples1,
  samples2,
  condition1,
  condition2,
  ouput_file_path,
  interaction_oi,
  condition
) {
  ###load and prepare data 
  cond1_df <- read.csv(csv_file_path_cond1)
  cond2_df <- read.csv(csv_file_path_cond2)
  
  cond1_df$X <- NULL
  cond2_df$X <- NULL
  
  ###DE analysis 
  samples1 <- as.data.frame(samples1)
  colnames(samples1) <- "sample"
  samples1$condition <- condition1
  samples2 <- as.data.frame(samples2)
  colnames(samples2) <- "sample"
  samples2$condition <- condition2
  
  #combine extra information 
  extra.info <- rbind(samples1,samples2)
  rownames(extra.info) <- extra.info$sample
  extra.info$sample <- NULL
  
  df_cond_both <- merge(cond1_df,cond2_df,by = "interacting_pair",all = TRUE)
  #put 0 instead of NA
  df_cond_both[is.na(df_cond_both)] <- 0
  
  #remove duplicated rows 
  df_cond_both <- df_cond_both[!duplicated(df_cond_both$interacting_pair),]
  
  #remove rows with 0 counts 
  df_cond_both$sum <-rowSums(df_cond_both[sapply(df_cond_both, is.numeric)], na.rm = TRUE)
  df_cond_both <- df_cond_both[df_cond_both$sum !=0,]
  df_cond_both$sum <- NULL
  
  #modify "count" table 
  rownames(df_cond_both) <- df_cond_both$interacting_pair
  df_cond_both$interacting_pair <- NULL
  cluster.counts <- df_cond_both
  cluster.counts <- as(cluster.counts, "matrix")
  
  ##Set up edgeR object
  y.ab <- DGEList(cluster.counts, samples = extra.info)
  
  y.ab$samples$condition <- factor(y.ab$samples$condition)
  
  ##Defining the model matrix
  mdl <- model.matrix(~ factor(condition),y.ab$samples)
  y <- estimateDisp(y.ab, mdl)
  fit <- glmQLFit(y, mdl, robust=TRUE)
  res <- glmQLFTest(fit)
  
  ##save results for plotting 
  top <- topTags(res,n=nrow(y))$table
  top$Interaction <- rownames(top)
  write.csv(top, paste0(ouput_file_path,interaction_oi,"_",condition,"_comp_LM.csv"))
}

###Function to compare two conditions of CellPhoneDB outputs 
L_R_CellPhoneDB_comp_2samples <- function(
  interacting_cells_oi,
  file1,
  file2,
  ouput_file_path,
  sample_name
){
  #read in significant.txt files from two conditions and merge the data frames, if no value = 0 
  df1 <- read.table(file1,header = TRUE, sep = "\t")
  #convert NA to 0 
  df1[is.na(df1)] <- 0
  df2 <- read.table(file2,header = TRUE, sep = "\t")
  df2[is.na(df2)] <- 0
  
  #merge all data frames to one containing interactions from all 
  #substract rows of interest 
  df1 <- df1[,which(names(df1) %in% c("interacting_pair",interacting_cells_oi))]
  df2 <- df2[,which(names(df2) %in% c("interacting_pair",interacting_cells_oi))]
  
  merged <- merge(df1, df2, by="interacting_pair", all=TRUE) 
  names(merged) <- c("interacting_pair","sample1","sample2")
  merged[is.na(merged)] <- 0
  
  merged$sum <-rowSums(merged[sapply(merged, is.numeric)], na.rm = TRUE)
  #remove rows with 0
  merged <-merged[merged$sum !=0,]
  merged$sum <- NULL
  write.csv(merged,paste0(ouput_file_path,"interactions_comp_two_cond_",interacting_cells_oi,"_",sample_name, ".csv"))
}



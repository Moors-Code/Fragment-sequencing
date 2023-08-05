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
  ###generating counts file 
  #take raw data and normalize it
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
  ###generating meta file based on cell type annotation of Seurat object 
  meta_data_meta <- cbind(rownames(seurat_object@meta.data), seurat_object@meta.data[,annotation_column, drop=F])  
  #save meta file as text file 
  write.table(meta_data_meta, paste0(ouput_file_path,sample_name,"_meta.txt"), sep='\t', quote=F, row.names=F)
}

###Function to compare two conditions of CellPhoneDB outputs, makes a combined dataframe with both samples + p-values 
L_R_CellPhoneDB_comp_2samples <- function(
  interacting_cells_oi,
  file1_mean,
  file1_pval,
  file2_mean,
  file2_pval,
  ouput_file_path,
  sample_name
){
  #read in significant.txt and pvalues.txt files from two conditions and merge the data frames, if no value = 0 
  df1_mean <- read.table(file1_mean,header = TRUE, sep = "\t")
  #convert NA to 0 
  df1_mean[is.na(df1_mean)] <- 0
  df2_mean <- read.table(file2_mean,header = TRUE, sep = "\t")
  df2_mean[is.na(df2_mean)] <- 0
  df1_pval <- read.table(file1_pval, header = TRUE, sep = "\t")
  df1_pval[is.na(df1_pval)] <- 0
  df2_pval <- read.table(file2_pval, header = TRUE, sep = "\t")
  df2_pval[is.na(df2_pval)] <- 0
  
  #merge interactions of an interacting cell pair of interest to one data frame 
  #substract rows of interest 
  df1_mean <- df1_mean[,which(names(df1_mean) %in% c("interacting_pair",interacting_cells_oi))]
  df2_mean <- df2_mean[,which(names(df2_mean) %in% c("interacting_pair",interacting_cells_oi))]
  
  df1_pval <- df1_pval[,which(names(df1_pval) %in% c("interacting_pair",interacting_cells_oi))]
  df2_pval <- df2_pval[,which(names(df2_pval) %in% c("interacting_pair",interacting_cells_oi))]
  
  #merge mean files 
  merged <- merge(df1_mean, df2_mean,by="interacting_pair", all=TRUE) 
  names(merged) <- c("interacting_pair","sample1_mean","sample2_mean")
  #remove interactions with 0 
  merged$sum <-rowSums(merged[sapply(merged, is.numeric)], na.rm = TRUE)
  #remove rows with 0
  merged <-merged[merged$sum !=0,]
  merged$sum <- NULL
  
  #merge with pvalue files, remove all=TRUE so only interactions are merged with >0 interaction scores 
  merged <- merge(merged, df1_pval, by="interacting_pair")
  names(merged) <- c("interacting_pair","sample1_mean","sample2_mean","sample1_pval")
  merged <- merge(merged, df2_pval, by="interacting_pair")
  names(merged) <- c("interacting_pair","sample1_mean","sample2_mean","sample1_pval","sample2_pval")
  merged[is.na(merged)] <- 0
  
  #save data frame for plotting 
  write.csv(merged,paste0(ouput_file_path,"interactions_comp_two_cond_",interacting_cells_oi,"_",sample_name, ".csv"))
}

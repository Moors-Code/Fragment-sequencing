###Function to extract info from x y coordinate files from highly multiplexed FISH images 
extracting_x_y_info_Resolve <- function(
  file_path
) {
  csv_file <- read.table(file_path,header = TRUE,sep = "\t")
  csv_file$coordinate <- NA
  csv_file$coordinate <- gsub("(.*),.*","\\1",csv_file$X.Y.Value)      
  output_vector <- csv_file$coordinate
}

###Function to generate text files with cell type annotation per slide for the overlay with the Dapi image in ImageJ 
#it produces text files for each slide or cell type with the cell number ID in the first column and the cell annotation in the second column 
text_file_generation_anno_resolve <- function(
  seurat_object,
  cell_type,
  slide_name,
  sample_name
) {
  Idents(seurat_object) <- "annotation"
  text_file <- subset(seurat_object, idents = cell_type)
  text_file <- as.data.frame(text_file@meta.data)
  text_file <- text_file[,c("Cell","annotation")]
  text_file$Cell <- gsub(slide_name, "Cell", rownames(text_file))
  write.table(text_file, paste0(sample_name, "_", cell_type, ".txt"),row.names = FALSE, quote = F, sep = "\t")
}

text_file_generation_per_slide_resolve <- function(
  seurat_object,
  slide_name,
  slide_name2
) {
  Idents(seurat_object) <- "Slide"
  text_file <- subset(seurat_object, idents = slide_name)
  text_file <- as.data.frame(text_file@meta.data)
  text_file <- text_file[,c(10,24)]
  text_file$Cell <- gsub(slide_name2, "Cell", rownames(text_file))
  write.table(text_file, paste0(slide_name, "_", ".txt"),row.names = FALSE, quote = F, sep = "\t")
}

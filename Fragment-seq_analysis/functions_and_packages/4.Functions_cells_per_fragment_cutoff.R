###Function that applies a custom cutoff to remove spheres with a low cell number 
Fragment_cell_cutoff <- function(
  seurat_object_file_path, 
  cell_cutoff,
  path_to_save_seurat_file
){
  seurat_object <- readRDS(seurat_object_file_path)
  df <- as.data.frame(table(seurat_object$fragment))
  df2 <- df[df$Freq >= cell_cutoff,]
  barcodes <- as.character(df2$Var1)
  Idents(seurat_object) <- "fragment"
  seurat_object2 <- subset(seurat_object, idents = barcodes)
  saveRDS(seurat_object2, file = path_to_save_seurat_file)
}  

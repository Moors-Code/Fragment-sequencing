###Function to extract relevant information from Biosorter output files 
Biosorter_output_managing <- function(
  biosorter_files, 
  seurat_object_prefix
){
  ###predict sphere size based on linear model fit and add to a column in dataframe 
  predicted_size <- predict.lm(lm.bead_tof, biosorter_files)
  #round sizes to integers 
  biosorter_files$sphere_size <- round(predicted_size)
  ###remove columns not needed 
  #Sorted status: 0 = just acquiring, 6= sorted (this we want!), 9= coincidence with previous, 10 = coincidence with following 
  #Clog should be N, remove if not 
  biosorter_files <- biosorter_files[biosorter_files$Sorted.status == "6" & biosorter_files$Clog == "N",]
  biosorter_files <- biosorter_files[,-which(names(biosorter_files) %in% c("Id","Time","Plate","Clog","Sorted.status", "Source.well","In.Regions",
                                                                           "TOF","Extinction","Violet","Red","PH.Extinction",
                                                                           "PW.Extinction","PC.Extinction","PH.Violet",
                                                                           "PW.Violet","PC.Violet","PH.Green","PW.Green",
                                                                           "PC.Green","PH.Red","PW.Red","PC.Red","X"))]
  ###add sphere ID by matching the Row and Column ID with the plate ID of MULTI-seq barcodes 
  sphere1 <- paste0("Bar", c(1:96), seurat_object_prefix)
  sphere2 <- paste0("Bar", c(97:192), seurat_object_prefix)
  sphere3 <- paste0("Bar", c(193:288), seurat_object_prefix)
  biosorter_id <- c(paste0("A ", c(1:12)),paste0("B ", c(1:12)), paste0("C ", c(1:12)),paste0("D ", c(1:12)),
                    paste0("E ", c(1:12)),paste0("F ", c(1:12)),paste0("G ", c(1:12)),paste0("H ", c(1:12)))
  
  plate1 <- data.frame(sphere1,biosorter_id)
  colnames(plate1) <- c("sphere","biosorter_id")
  plate2 <- data.frame(sphere2,biosorter_id)
  colnames(plate2) <- c("sphere","biosorter_id")
  plate3 <- data.frame(sphere3,biosorter_id)
  colnames(plate3) <- c("sphere","biosorter_id")
  
  #combine Row and Column (f.e. A1 in the end) in dataframe to then merge it with the sphere ids 
  biosorter_files$biosorter_id <- paste(biosorter_files$Row, biosorter_files$Column)
  
  #subset per plate 
  biosorter_files_plate1 <- subset(biosorter_files, plate_id == "plate_1")
  biosorter_files_plate2 <- subset(biosorter_files, plate_id == "plate_2")
  biosorter_files_plate3 <- subset(biosorter_files, plate_id == "plate_3")
  
  #add sphere_id to the biosorter data by matching with plate ids 
  biosorter_files_plate1 <- merge(biosorter_files_plate1, plate1)
  biosorter_files_plate2 <- merge(biosorter_files_plate2, plate2)
  biosorter_files_plate3 <- merge(biosorter_files_plate3, plate3)
  
  #merge to one file 
  biosorter_files_all_plates <- full_join(biosorter_files_plate1, biosorter_files_plate2)
  biosorter_files_all_plates <- full_join(biosorter_files_all_plates,biosorter_files_plate3)
  return(biosorter_files_all_plates)
}  

###Function to integrate biosorter sphere-size information into Seurat object 
Biosorter_data_seurat_integration <- function(
  seurat_object_file_path,
  biosorter_files,
  output_Seurat_object_file_path_name
){
  #add an empty column in meta.data for sphere-size 
  seurat_object <- readRDS(file = seurat_object_file_path)
  seurat_object$sphere_size <- NA
  seurat_object$GFP <- NA
  #match with sphere-size and Green (GFP) of biosorter files 
  seurat_object$sphere_size <- biosorter_files$sphere_size[match(seurat_object$sphere,biosorter_files$sphere)]
  seurat_object$GFP <- biosorter_files$Green[match(seurat_object$sphere,biosorter_files$sphere)]
  #normalize GFP with sphere_size 
  seurat_object$GFP <- as.character(seurat_object$GFP)
  seurat_object$GFP <- as.numeric(seurat_object$GFP)
  seurat_object$sphere_size <- as.character(seurat_object$sphere_size)
  seurat_object$sphere_size <- as.numeric(seurat_object$sphere_size)
  seurat_object$GFP.norm <- (seurat_object$GFP / seurat_object$sphere_size)
  #save R object 
  saveRDS(seurat_object, file = output_Seurat_object_file_path_name)
}    

###Function to produce data frames for boxplot plotting of sphere size signal 
BS_df_for_boxplot_per_sample <- function(
  sub_seurat_obj, 
  ident_of_interest1,
  ident_of_interest2,
  sample_name
){
  df <- as.data.frame(table(sub_seurat_obj@meta.data[[ident_of_interest1]],sub_seurat_obj@meta.data[[ident_of_interest2]]))
  df <- filter(df, df$Freq > 1)
  df <- df[,1:2]
  colnames(df) <- c(ident_of_interest1,ident_of_interest2)
  df$sample <- sample_name
  df$sample <- as.factor(df$sample)
  df[,2] <- as.character(df[,2])
  df[,2] <- as.numeric(df[,2])
  return(df)
}


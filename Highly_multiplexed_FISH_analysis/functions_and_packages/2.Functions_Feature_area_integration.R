###Function to extract info from x y coordinate files from highly multiplexed FISH images 
extracting_x_y_info_Resolve <- function(
  file_path
) {
  csv_file <- read.table(file_path,header = TRUE,sep = "\t")
  csv_file$coordinate <- NA
  csv_file$coordinate <- gsub("(.*),.*","\\1",csv_file$X.Y.Value)      
  output_vector <- csv_file$coordinate
}

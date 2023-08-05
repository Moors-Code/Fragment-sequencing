###Function to add zonation coordinate to Seurat object, only take LECs into account
add_zonation_coefficient <- function(
  seurat_object_rds_file_path,
  output_seurat_file_name
){
  ###load Seurat object
  seurat_object <- readRDS(seurat_object_rds_file_path)
  #subset only LECs
  Idents(seurat_object) <- "annotation.broad"
  endo <- subset(seurat_object, idents = c("LECs"))
  ###Generate pseudobulks of spheres with portal and central landmark genes 
  #genes used from Halpern et al 2018, Paired-cell seq paper 
  PV_genes <- c("Sh3rf1","Ppp1r9b","Ctdspl","Serpina3f", "Impdh1","Meis1","Lama4", "Armcx1",
                "Txndc16","Klhl3","Lpcat1","Tnik","Gdf2","Trim30a", "Rnf157", "Lyve1","Pcolce",
                "Osmr","Enpp6","Pla2r1","Pygb", "Pld1","Ntm","Il33","Igfbp3","Hoxb5","Art3","Tgfb1",
                "B4galt4","Klf4","Sod3","Fgfr1", "Ly6a","Hlx","Slc43a2", "Rapgef3","Sept8","Itgb3",
                "Efnb2","Twist1","Sox17","Itga9", "Rasd1", "Slc41a1","Spats2l","Gata2", "Col1a2",
                "Chst2","Nid2","Tmem44","Adam23","Prkcq","Tomm40l", "Dll4", "Il1a","Ntn4", "Msr1")
  
  CV_genes <- c("Wnt9b", "Rspo3", "Cdh13","Thbd","Lmcd1","Ar", "Ier3","Cc2d2a", "Dkk3",
                "Fabp4","Dennd2d","Cln6","Lfng","Cebpd","Tmed8","Rgp1", "Pfkfb3","Wnt2", "Robo1",
                "Fam84b","Rab3b", "Ptgs1","P2ry1","Olfm1", "Amigo2","Trappc2", "Kit",
                "Mpp1","Ap4m1","Nkd1","Sgsh","Zxdb","Mbd5","Atr","Slc48a1", "Gen1","Pigh","Unkl","Itgb3bp",
                "Zmat1","Med9","Anxa3","Kcnb1","Gja4","Eid1","Arl16","Orc2","Hs1bp3","Gorasp1",
                "Skp1a", "Med21", "Ptcd1","Rpa1","Dbp","Tsr2","Diablo")
  
  #Average expression of CV and PV genes in spheres, take raw counts and log normalize them beforehand
  #log normalized counts will be stored in slot "data"
  #raw data in slot "counts"
  endo_norm <- NormalizeData(endo, normalization.method = "LogNormalize",
                             scale.factor = 10000,
                             margin = 1, assay = "RNA")
  endo_av_PV <- AverageExpression(endo_norm, features = PV_genes,assays = "RNA", slot = "data", group.by = "sphere")
  endo_av_CV <- AverageExpression(endo_norm, features = CV_genes,assays = "RNA", slot = "data", group.by = "sphere")
  
  #convert to matrix 
  endo_av_PV_df <- as.data.frame(endo_av_PV)
  endo_av_CV_df <- as.data.frame(endo_av_CV)
  
  ###Calculation of zonation coordinate (ZC)
  #normalize second time by dividing each gene by the maximum level from all spheres, this makes sure that each gene gets the same weight   
  output_endo_av_PV_df <- matrix(0,nrow(endo_av_PV_df), ncol(endo_av_PV_df), dimnames = list(rownames(endo_av_PV_df), colnames(endo_av_PV_df)))  
  for (i in 1:nrow(endo_av_PV_df)){          
    for (j in 1:ncol(endo_av_PV_df)){        
      original.val <- endo_av_PV_df[i,j]         
      max_umi_count <- max(endo_av_PV_df[i,])     
      adapted_value <- original.val/max_umi_count            
      output_endo_av_PV_df[i,j] <- adapted_value                            
    }
  }
  
  output_endo_av_CV_df <- matrix(0,nrow(endo_av_CV_df), ncol(endo_av_CV_df), dimnames = list(rownames(endo_av_CV_df), colnames(endo_av_CV_df)))  
  for (i in 1:nrow(endo_av_CV_df)){          
    for (j in 1:ncol(endo_av_CV_df)){        
      original.val <- endo_av_CV_df[i,j]         
      max_umi_count <- max(endo_av_CV_df[i,])     
      adapted_value <- original.val/max_umi_count            
      output_endo_av_CV_df[i,j] <- adapted_value                            
    }
  }
  
  #calculate colSums to get sum of counts for CV and PV genes from each sphere
  endo_av_PV_df_total <- rbind(output_endo_av_PV_df, Total = colSums(output_endo_av_PV_df))
  endo_av_CV_df_total <- rbind(output_endo_av_CV_df, Total = colSums(output_endo_av_CV_df))
  
  #Calculate 
  Zonation_score_endo_per_sphere <- endo_av_PV_df_total["Total",]/(endo_av_PV_df_total["Total",]+endo_av_CV_df_total["Total",])
  
  #Rescale ZCs to that 0 is most central, 1 is most portal (ZC = (ZC-min/ZC)) / (max(ZC)-min(ZC))
  minimal_zC <- min(Zonation_score_endo_per_sphere)
  maximal_zC <- max(Zonation_score_endo_per_sphere)
  Zonation_score_endo_per_sphere <- sort(Zonation_score_endo_per_sphere)
  Zonation_score_endo_per_sphere.df <- as.data.frame(Zonation_score_endo_per_sphere)
  colnames(Zonation_score_endo_per_sphere.df) <- "zonation_coordinate"
  Zonation_score_endo_per_sphere.df2 <- matrix(0,nrow(Zonation_score_endo_per_sphere.df), ncol(Zonation_score_endo_per_sphere.df), 
                                               dimnames = list(rownames(Zonation_score_endo_per_sphere.df), colnames(Zonation_score_endo_per_sphere.df)))  
  for (i in 1:nrow(Zonation_score_endo_per_sphere.df)){          
    original.val <- Zonation_score_endo_per_sphere.df[i,]         
    adapted_value <- (original.val-minimal_zC)/(maximal_zC-minimal_zC)           
    Zonation_score_endo_per_sphere.df2[i] <- adapted_value                            
  }
  Zonation_score_endo_per_sphere.df2 <- as.data.frame(Zonation_score_endo_per_sphere.df2)
  
  ###Add ZC to Seurat object 
  #remove the 'RNA.' from the Barcode name 
  rownames(Zonation_score_endo_per_sphere.df2) <- gsub(pattern = "RNA.", x = rownames(Zonation_score_endo_per_sphere.df2), replacement = "")
  
  current.cluster.ids <- rownames(Zonation_score_endo_per_sphere.df2)
  new.cluster.ids <- as.numeric(Zonation_score_endo_per_sphere.df2$zonation_coordinate)
  seurat_object$zonation_coordinate <- plyr::mapvalues(x = seurat_object$sphere, from = current.cluster.ids, to = new.cluster.ids)
  #save Seurat object 
  saveRDS(seurat_object,file = output_seurat_file_name)
  return(seurat_object)
}  

###Function to remove spheres without LECs because there is no ZC
remove_spheres_noLECs <- function(
  seurat_object_path,
  seurat_object_new_path
){
  seurat_object <- readRDS(seurat_object_path)
  #check for spheres that don't have zonation coordinate 
  vector_all_bc <- as.character(as.data.frame(table(seurat_object$sphere))$Var1)
  vector_all_zC <- as.character(as.data.frame(table(seurat_object$zonation_coordinate))$Var1)
  vector_no_zC_bc <- grep("^B", vector_all_zC, value=TRUE)
  vector_zC_bc <- setdiff(vector_all_bc,vector_no_zC_bc)
  Idents(seurat_object) <- "sphere"
  sub <- subset(seurat_object, idents =vector_no_zC_bc )
  print(table(sub$annotation))
  seurat_object <- subset(seurat_object, idents = vector_zC_bc)
  saveRDS(seurat_object,seurat_object_new_path)
}

###Function to add lobule layer name to Seurat object 
add_lobule_zones <- function(
  seurat_object_file_path,
  output_seurat_file_name
){
  seurat_object <- readRDS(seurat_object_file_path)
  ###add 10 lobule layers first and then combine L1-L3 and L8-L10
  seurat_object$lobule10 <- NA
  for(i in 1:nrow(seurat_object@meta.data)){
    row <- seurat_object@meta.data[i,]
    if(row$zonation_coordinate <=0.1 & row$zonation_coordinate >=0 ) {
      row$lobule10 <- "L1" 
    }else if (row$zonation_coordinate <=0.2 & row$zonation_coordinate >0.1 ){
      row$lobule10 <- "L2"
    }else if (row$zonation_coordinate <=0.3 & row$zonation_coordinate >0.2 ){
      row$lobule10 <- "L3"
    }else if (row$zonation_coordinate <=0.4 & row$zonation_coordinate >0.3 ){
      row$lobule10 <- "L4"
    }else if (row$zonation_coordinate <=0.5 & row$zonation_coordinate >0.4 ){
      row$lobule10 <- "L5"
    }else if (row$zonation_coordinate <=0.6 & row$zonation_coordinate >0.5 ){
      row$lobule10 <- "L6"
    }else if (row$zonation_coordinate <=0.7 & row$zonation_coordinate >0.6 ){
      row$lobule10 <- "L7"
    }else if (row$zonation_coordinate <=0.8 & row$zonation_coordinate >0.7 ){
      row$lobule10 <- "L8"
    }else if (row$zonation_coordinate <=0.9 & row$zonation_coordinate >0.8 ){
      row$lobule10 <- "L9"
    }else if (row$zonation_coordinate <=1 & row$zonation_coordinate >0.9 ){
      row$lobule10 <- "L10"
    }else{
      row$lobule10 <- NA
    }
    seurat_object@meta.data[i,] <- row
  }
  current.cluster.ids <- c("L1","L2","L3","L4","L5","L6","L7","L8","L9","L10")
  new.cluster.ids <- c("L1-L3","L1-L3","L1-L3","L4","L5","L6","L7","L8-L10","L8-L10","L8-L10")
  seurat_object$lobule_layer <- plyr::mapvalues(x = seurat_object$lobule10, from = current.cluster.ids, to = new.cluster.ids)
  ###Add CV and PV information by combining L1-L5 and L6-L10
  current.cluster.ids <- c("L1","L2","L3","L4","L5","L6","L7","L8","L9","L10")
  new.cluster.ids <- c("CV","CV","CV","CV","CV","PV","PV","PV","PV","PV")
  seurat_object$vein <- plyr::mapvalues(x = seurat_object$lobule10, from = current.cluster.ids, to = new.cluster.ids)
  ###save R object 
  saveRDS(seurat_object,file = output_seurat_file_name)
  return(seurat_object)
}


###Function to analyse differentially expressed genes along the lobule layer axis per cell type of interest 
#adapted from Karsten Bach 
DE_zonated_genes <- function(
  seurat_object,
  ident_column_oi, #annotation or annotation.broad
  ident_oi, #cell type of interest 
  output_file_path,
  sample_name
) {
  ###load and prepare data 
  Idents(seurat_object) <- ident_column_oi
  seurt <- subset(seurat_object, idents = ident_oi)
  
  #convert Seurat object to SCE object with required meta data information 
  m <- GetAssayData(seurt, assay = "RNA", slot = "counts")
  pD <- data.frame("barcode"=colnames(m),
                   "Sphere"=seurt$sphere,
                   "Layer"=seurt$lobule_layer,
                   "Batch"=seurt$orig.ident)
  
  sce <- SingleCellExperiment(assays=list(counts=m),colData=DataFrame(pD))
  
  #Sum across the spheres
  sumd <- aggregateAcrossCells(sce,ids=sce$Sphere)
  #Only consider spheres with at least 5 cells 
  sumd <- sumd[,sumd$ncells >= 5]
  
  ###Set up edgeR object 
  y <- DGEList(counts=counts(sumd),samples=data.frame(colData(sumd)))
  
  #Remove lowly expressed genes
  keep <- filterByExpr(y, group=y$samples$Layer)
  y <- y[keep,]
  
  #Normalization
  y <- calcNormFactors(y)
  
  ###Defining the model matrix
  #First we need to transform the layer column into an ordered factor
  y$samples$Layer <- factor(y$samples$Layer)
  y$samples$Layer <- ordered(y$samples$Layer)
  #We are interested in linear changes that go up/down from L1-L3-L8-L10 
  #account for variability in samples with '~ Batch + ...'
  mdl <- model.matrix(~Batch + Layer,y$samples)
  
  ###Follow the standard edgeR workflow
  y <- estimateDisp(y, mdl)
  fit <- glmQLFit(y, mdl, robust=TRUE)
  
  #Specify liner 
  res <- glmQLFTest(fit, coef="Layer.L")
  
  ###save results for plotting 
  top <- topTags(res,n=nrow(y))$table
  top$Gene <- rownames(top)
  write.csv(top, paste0(output_file_path,sample_name,"_", ident_oi, "_LM_zonated_genes.csv"))
}

###Function to plot genes of interest in LECs across the lobule layer axis in injected samples, you can save as .svg or .pdf 
boxplot_zonation_genes_LECs_injected <- function(
  seurat_object,
  gene_oi,
  color_oi,
  output_file_path,
  sample_name,
  save_as
){
  Idents(seurat_object) <-"annotation.broad"
  cellType <- subset(seurat_object, idents = "LECs")
  #take normalized counts
  cellType_norm <- NormalizeData(cellType, normalization.method = "LogNormalize",
                                 scale.factor = 10000,
                                 margin = 1, assay = "RNA")
  #per zonation coordinate = per sphere 
  AE <- AverageExpression(cellType_norm,assays = "RNA", features = gene_oi,slot = "data", group.by = "zonation_coordinate")
  AE_df <- as.data.frame(AE$RNA)
  AE_df2 <- data.frame(names(AE_df), as.numeric(AE_df[1,]))
  names(AE_df2) <- c("zonation_coordinate",gene_oi)
  AE_df2[is.na(AE_df2)] <- 0
  m <- reshape2::melt(AE_df2, id.vars = "zonation_coordinate", meature.vars = gene_oi)
  #specify the position of different lobule layers based on ZCs and add in an additional column 
  m$lobule_layer <- NA
  m[1:137,]$lobule_layer <- "L1-L3" 
  m[138:351,]$lobule_layer <- "L4" 
  m[352:753,]$lobule_layer <- "L5" 
  m[754:1162,]$lobule_layer <- "L6" 
  m[1163:1358,]$lobule_layer <- "L7" 
  m[1359:1384,]$lobule_layer <- "L8-L10" 
  
  #plot boxplot 
  p <- ggplot(m, aes(x=lobule_layer, y=value, fill=variable)) +  theme_classic() +
    geom_boxplot(fill = color_oi,outlier.shape = NA) + theme(axis.text = element_text(size = 30))  + theme(axis.text.x = element_text(angle = 90)) +
    geom_jitter(color="black", size=1, alpha=0.9) +theme(axis.title= element_text(size = 25)) + 
    ggtitle(paste0(gene_oi," expression in LECs - sphere-seq")) + xlab("Lobule layer (CV-PV)") + 
    ylab("Average expression per sphere") +
    theme(plot.title = element_text(size = 25, face = "bold")) + theme(legend.text = element_text(size = 30),
                                                                       legend.title= element_text(size = 30)) + 
    guides(fill=guide_legend(title="Gene"))  + 
    guides(fill=guide_legend(title="Gene"))
  p + ggsave(paste0(output_file_path,gene_oi,"_" ,sample_name,"_LECs_boxplot",save_as),width = 12, height = 10)
}

###Function to plot genes of interest in Kupffer cells across lobule layer axis in injected samples, you can save as .svg or .pdf 
boxplot_zonation_genes_KCs_injected <- function(
  seurat_object,
  gene_oi,
  color_oi,
  output_file_path,
  sample_name,
  save_as
){
  Idents(seurat_object) <-"annotation.broad"
  cellType <- subset(seurat_object, idents = "Kupffer")
  #take normalized counts
  cellType_norm <- NormalizeData(cellType, normalization.method = "LogNormalize",
                                 scale.factor = 10000,
                                 margin = 1, assay = "RNA")
  #per zonation coordinate = per sphere 
  AE <- AverageExpression(cellType_norm,assays = "RNA", features = gene_oi,slot = "data", group.by = "zonation_coordinate")
  AE_df <- as.data.frame(AE$RNA)
  AE_df2 <- data.frame(names(AE_df), as.numeric(AE_df[1,]))
  names(AE_df2) <- c("zonation_coordinate",gene_oi)
  AE_df2[is.na(AE_df2)] <- 0
  m <- reshape2::melt(AE_df2, id.vars = "zonation_coordinate", meature.vars = gene_oi)
  #specify the position of different lobule layers based on ZCs and add in an additional column 
  m$lobule_layer <- NA
  m[1:115,]$lobule_layer <- "L1-L3" 
  m[116:321,]$lobule_layer <- "L4" 
  m[322:700,]$lobule_layer <- "L5" 
  m[701:1087,]$lobule_layer <- "L6" 
  m[1088:1269,]$lobule_layer <- "L7" 
  m[1270:1295,]$lobule_layer <- "L8-L10" 
  
  #plot boxplot 
  p <- ggplot(m, aes(x=lobule_layer, y=value, fill=variable)) +  theme_classic() +
    geom_boxplot(fill = color_oi,outlier.shape = NA) + theme(axis.text = element_text(size = 30))  + theme(axis.text.x = element_text(angle = 90)) +
    geom_jitter(color="black", size=1, alpha=0.9) +theme(axis.title= element_text(size = 25)) + 
    ggtitle(paste0(gene_oi," expression in KCs - sphere-seq")) + xlab("Lobule layer (CV-PV)") + 
    ylab("Average expression per sphere") +
    theme(plot.title = element_text(size = 25, face = "bold")) + theme(legend.text = element_text(size = 30),
                                                                       legend.title= element_text(size = 30)) + 
    guides(fill=guide_legend(title="Gene"))  + 
    guides(fill=guide_legend(title="Gene"))
  p + ggsave(paste0(output_file_path,gene_oi,"_" ,sample_name,"_KCs_boxplot",save_as),width = 12, height = 10)
}


###Function to plot genes of interest in LECs of WT sample between CV and PV + applying wilcox rank sum test 
#you can save as pdf or svg
boxplot_zonation_genes_LECs_WT <- function(
  seurat_object,
  gene_oi,
  color_oi,
  output_file_path,
  sample_name,
  save_as
){
  Idents(seurat_object) <-"annotation.broad"
  cellType <- subset(seurat_object, idents = "LECs")
  #take normalized counts
  cellType_norm <- NormalizeData(cellType, normalization.method = "LogNormalize",
                                 scale.factor = 10000,
                                 margin = 1, assay = "RNA")
  #per zonation coordinate = per sphere 
  AE <- AverageExpression(cellType_norm,assays = "RNA", features = gene_oi,slot = "data", group.by = "zonation_coordinate")
  AE_df <- as.data.frame(AE$RNA)
  AE_df2 <- data.frame(names(AE_df), as.numeric(AE_df[1,]))
  names(AE_df2) <- c("vein_area",gene_oi)
  AE_df2[is.na(AE_df2)] <- 0
  m <- reshape2::melt(AE_df2, id.vars = "vein_area", meature.vars = gene_oi)
  #specify the position of different vein areas based on ZCs and add in an additional column 
  m$vein <- NA
  m[c(1:76),]$vein <- "CV" 
  m[c(77:159),]$vein <- "PV" 
  
  #plot boxplot with wilcox test
  p <- ggplot(m, aes(x=vein, y=value, fill=variable)) + theme_classic() +
    geom_boxplot(fill = color_oi,outlier.shape = NA) + theme(axis.text = element_text(size = 30))  + theme(axis.text.x = element_text(angle = 90)) +
    geom_jitter(color="black", size=1, alpha=0.9) +
    ggtitle(paste0(gene_oi," expression in LECs - Sphere-seq WT")) + xlab("Lobule layer") + 
    ylab("Average expression per sphere") + theme(axis.title= element_text(size = 25)) +
    theme(plot.title = element_text(size = 25, face = "bold")) + theme(legend.text = element_text(size = 30),
                                                                       legend.title= element_text(size = 30)) + 
    guides(fill=guide_legend(title="Gene")) + ggsignif::geom_signif(comparisons = list(c("CV", "PV")),
                                                                    textsize=7,
                                                                    test = "wilcox.test",map_signif_level = c("***"=0.001,"**"=0.01,"*"=0.05))
  p + ggsave(paste0(output_file_path,gene_oi,"_" ,sample_name,"_LECs_boxplot",save_as),width = 12, height = 10)
}

###Function to plot genes of interest in Kupffer cells of WT sample between CV and PV + applying wilcox rank sum test
#you can save as pdf or svg
boxplot_zonation_genes_KCs_WT <- function(
  seurat_object,
  gene_oi,
  color_oi,
  output_file_path,
  sample_name,
  save_as
){
  Idents(seurat_object) <-"annotation.broad"
  cellType <- subset(seurat_object, idents = "Kupffer")
  #take normalized counts
  cellType_norm <- NormalizeData(cellType, normalization.method = "LogNormalize",
                                 scale.factor = 10000,
                                 margin = 1, assay = "RNA")
  #per zonation coordinate = per sphere 
  AE <- AverageExpression(cellType_norm,assays = "RNA", features = gene_oi,slot = "data", group.by = "zonation_coordinate")
  AE_df <- as.data.frame(AE$RNA)
  AE_df2 <- data.frame(names(AE_df), as.numeric(AE_df[1,]))
  names(AE_df2) <- c("vein_area",gene_oi)
  AE_df2[is.na(AE_df2)] <- 0
  m <- reshape2::melt(AE_df2, id.vars = "vein_area", meature.vars = gene_oi)
  #specify the position of different vein areas based on ZCs and add in an additional column 
  m$vein <- NA
  m[1:67,]$vein <- "CV" 
  m[68:141,]$vein <- "PV" 
  
  #plot boxplot with wilcox test
  p <- ggplot(m, aes(x=vein, y=value, fill=variable)) + theme_classic() +
    geom_boxplot(fill = color_oi,outlier.shape = NA) + theme(axis.text = element_text(size = 30))  + theme(axis.text.x = element_text(angle = 90)) +
    geom_jitter(color="black", size=1, alpha=0.9) +
    ggtitle(paste0(gene_oi," expression in Kupffer - Sphere-seq WT")) + xlab("Lobule layer") + 
    ylab("Average expression per sphere") + theme(axis.title= element_text(size = 25)) +
    theme(plot.title = element_text(size = 25, face = "bold")) + theme(legend.text = element_text(size = 30),
                                                                       legend.title= element_text(size = 30)) + 
    guides(fill=guide_legend(title="Gene")) + ggsignif::geom_signif(comparisons = list(c("CV", "PV")),
                                                                    textsize=7,
                                                                    test = "wilcox.test",map_signif_level = c("***"=0.001,"**"=0.01,"*"=0.05))
  p + ggsave(paste0(output_file_path,gene_oi,"_" ,sample_name,"_KCs_boxplot",save_as),width = 12, height = 10)
}


###Function that used edgeR to test differentially expressed genes between central and portal veins 
#Code adapted from Karsten Bach
DGE_between_veins_hmFISH <- function(
  seurat_object,
  ident_column_oi,
  ident_oi,
  condition,
  output_file_path
) {
  ###load and prepare data 
  Idents(seurat_object) <- ident_column_oi
  seurt <- subset(seurat_object, idents = ident_oi)
  
  #convert Seurat object to SCE object with required meta data information 
  m <- GetAssayData(seurt, assay = "originalexp", slot = "counts")
  pD <- data.frame("barcode"=colnames(m),
                   "spatial_feature_number"=seurt$spatial_feature_number,
                   "Condition"=seurt$vein,
                   "Batch"=seurt$sampleID)
  
  sce <- SingleCellExperiment(assays=list(counts=m),colData=DataFrame(pD))
  
  #Sum across the spatial_feature_number
  sumd <- aggregateAcrossCells(sce,ids=sce$spatial_feature_number)
  #Only consider spheres with at least 5 cells
  sumd <- sumd[,sumd$ncells >= 5]
  
  ###Set up edgeR object
  y <- DGEList(counts=counts(sumd),samples=data.frame(colData(sumd)))
  
  #Remove lowly expressed genes
  keep <- filterByExpr(y, group=y$samples$Condition)
  y <- y[keep,]
  
  #Normalization
  y <- calcNormFactors(y)
  
  ###Defining the model matrix
  y$samples$Condition <- factor(y$samples$Condition)
  #account for variability in samples with '~ Batch + ...'
  mdl <- model.matrix(~Batch + Condition,y$samples)
  
  #Follow the standard edgeR workflow
  y <- estimateDisp(y, mdl)
  fit <- glmQLFit(y, mdl, robust=TRUE)
  res <- glmQLFTest(fit)
  
  ###save results for plotting
  top <- topTags(res,n=nrow(y))$table
  top$Gene <- rownames(top)
  write.csv(top, paste0(output_file_path,condition,"_",ident_oi,"_DGE_hmFISH.csv"))
}

###Function that used edgeR to test differentially expressed genes between proximal and distal areas from metastatic sites  
#Code adapted from Karsten Bach 
DGE_between_Mets_distance_hmFISH <- function(
  seurat_object,
  ident_column_oi,
  ident_oi,
  condition,
  output_file_path
) {
  ###load and prepare data 
  Idents(seurat_object) <- ident_column_oi
  seurt <- subset(seurat_object, idents = ident_oi)
  
  #convert Seurat object to SCE object with required meta data information 
  m <- GetAssayData(seurt, assay = "originalexp", slot = "counts")
  pD <- data.frame("barcode"=colnames(m),
                   "spatial_feature_number"=seurt$spatial_feature_number,
                   "Condition"=seurt$Mets_distance,
                   "Batch"=seurt$sampleID)
  
  sce <- SingleCellExperiment(assays=list(counts=m),colData=DataFrame(pD))
  
  #Sum across the spheres
  sumd <- aggregateAcrossCells(sce,ids=sce$spatial_feature_number)
  #Only consider spheres with at least 5 cells
  sumd <- sumd[,sumd$ncells >= 5]
  
  ###Set up edgeR object
  y <- DGEList(counts=counts(sumd),samples=data.frame(colData(sumd)))
  
  #Remove lowly expressed genes
  keep <- filterByExpr(y, group=y$samples$Condition)
  y <- y[keep,]
  
  #Normalization
  y <- calcNormFactors(y)
  
  ###Defining the model matrix
  y$samples$Condition <- factor(y$samples$Condition)
  #account for variability in samples with '~ Batch + ...'
  mdl <- model.matrix(~Batch + Condition,y$samples)
  
  #Follow the standard edgeR workflow
  y <- estimateDisp(y, mdl)
  fit <- glmQLFit(y, mdl, robust=TRUE)
  res <- glmQLFTest(fit)
  
  ###save results for plotting
  top <- topTags(res,n=nrow(y))$table
  top$Gene <- rownames(top)
  write.csv(top, paste0(output_file_path,condition,"_",ident_oi,"_DGE_hmFISH.csv"))
}


###Function that produces boxplot of LECs 
#with vein on the x axis and on the y axis average gene expression per spatial feature area 
#storage of output as pdf or svg 
boxplot_LECs_per_spatial_area <- function(
  seurat_object, 
  gene_oi,
  color_oi,
  output_file_path,
  save_as
){
  Idents(seurat_object) <-"annotation"
  cellType <- subset(seurat_object, idents = "LECs")
  #take normalized counts
  cellType_norm <- NormalizeData(cellType, normalization.method = "LogNormalize",
                                 scale.factor = 10000,
                                 margin = 1, assay = "originalexp")
  #generate pseudobulks per spatial feature area
  AE <- AverageExpression(cellType_norm,assays = "originalexp", features = gene_oi,slot = "data", group.by = "spatial_feature_number")
  AE_df <- as.data.frame(AE$originalexp)
  AE_df2 <- data.frame(names(AE_df), as.numeric(AE_df[1,]))
  names(AE_df2) <- c("vein_area",gene_oi)
  AE_df2[is.na(AE_df2)] <- 0
  m <- reshape2::melt(AE_df2, id.vars = "vein_area", meature.vars = gene_oi)
  #specify the position of CV and PV and add as condition in an additional column 
  m$vein <- NA
  m[c(1:23,40:46,53:56,62:74,83:89,94:128),]$vein <- "CV" 
  m[c(24:39,47:52,57:61,75:82,90:93,129:155),]$vein <- "PV" 
  
  #plot in boxplot 
  ggplot(m, aes(x=vein, y=value, fill=variable)) + theme_classic() +
    geom_boxplot(fill = color_oi,outlier.shape = NA) + theme(axis.text = element_text(size = 30))  + theme(axis.text.x = element_text(angle = 90)) +
    geom_jitter(color="black", size=1, alpha=0.9) +
    ggtitle(paste0(gene_oi," expression in LECs - Resolve")) + xlab("Lobule layer") + 
    ylab("Average expression per spatial area") + theme(axis.title= element_text(size = 25)) +
    theme(plot.title = element_text(size = 25, face = "bold")) + theme(legend.text = element_text(size = 30),
                                                                       legend.title= element_text(size = 30)) + 
    guides(fill=guide_legend(title="Gene"))  +
    ggsave(paste0(output_file_path,gene_oi,"_LECs_boxplot",save_as), width = 12, height = 10)
}

###Function that produces boxplot of Kupffer cells 
#with vein on the x axis and on the y axis average gene expression per spatial feature area 
#storage of output as pdf or svg 
boxplot_KCs_per_spatial_area <- function(
  seurat_object, 
  gene_oi,
  color_oi,
  output_file_path,
  save_as
){
  Idents(seurat_object) <-"annotation"
  cellType <- subset(seurat_object, idents = "Kupffer")
  #take normalized counts
  cellType_norm <- NormalizeData(cellType, normalization.method = "LogNormalize",
                                 scale.factor = 10000,
                                 margin = 1, assay = "originalexp")
  #generate pseudobulks per spatial feature area
  AE <- AverageExpression(cellType_norm,assays = "originalexp", features = gene_oi,slot = "data", group.by = "spatial_feature_number")
  AE_df <- as.data.frame(AE$originalexp)
  AE_df2 <- data.frame(names(AE_df), as.numeric(AE_df[1,]))
  names(AE_df2) <- c("vein_area",gene_oi)
  AE_df2[is.na(AE_df2)] <- 0
  m <- reshape2::melt(AE_df2, id.vars = "vein_area", meature.vars = gene_oi)
  #specify the position of CV and PV and add as condition in an additional column
  m$vein <- NA
  m[c(1:23,40:46,53:56,62:74,83:89,94:128),]$vein <- "CV" 
  m[c(24:39,47:52,57:61,75:82,90:93,129:155),]$vein <- "PV" 
  
  #plot in boxplot 
  ggplot(m, aes(x=vein, y=value, fill=variable)) + theme_classic() +
    geom_boxplot(fill = color_oi,outlier.shape = NA) + theme(axis.text = element_text(size = 30))  + theme(axis.text.x = element_text(angle = 90)) +
    geom_jitter(color="black", size=1, alpha=0.9) +
    ggtitle(paste0(gene_oi," expression in KCs - Resolve")) + xlab("Lobule layer") + 
    ylab("Average expression per spatial area") + theme(axis.title= element_text(size = 25)) +
    theme(plot.title = element_text(size = 25, face = "bold")) + theme(legend.text = element_text(size = 30),
                                                                       legend.title= element_text(size = 30)) + 
    guides(fill=guide_legend(title="Gene"))  +
    ggsave(paste0(output_file_path,gene_oi,"_KCs_boxplot",save_as), width = 12, height = 10)
}


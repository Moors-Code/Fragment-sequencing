###Function to test cell type abundances between two groups, in this case distal vs. proximal  
#following workflow of: http://bioconductor.org/books/3.15/OSCA.multisample/differential-abundance.html#performing-the-da-analysis 
DA_analysis_cell_type_abundance_Mets_distance_hmFISH <- function(
  seurat_object,
  output_file_path
) {
  ###load and prepare data 
  #add extra meta data column that combines cell type and Mets_distance
  seurat_object$sample_mets_dist <- paste0(seurat_object$sampleID,"_",seurat_object$Mets_distance)
  #convert to SCE object 
  sce <- as.SingleCellExperiment(seurat_object)
  
  cluster.counts <- table(seurat_object$annotation,seurat_object$sample_mets_dist)
  cluster.counts <- as(cluster.counts, "matrix")
  #add extra information 
  extra.info <- colData(sce)[match(colnames(cluster.counts),sce$sample_mets_dist),]
  y.ab <- DGEList(cluster.counts, samples = extra.info)
  
  y.ab$samples$Mets_distance <- factor(y.ab$samples$Mets_distance)
  
  ###Defining the model matrix
  #account for variability in samples with '~ factor(sampleID) + ...'
  mdl <- model.matrix(~factor(sampleID) + factor(Mets_distance),y.ab$samples)
  
  ###Follow the standard edgeR workflow
  y <- estimateDisp(y.ab, mdl)
  fit <- glmQLFit(y, mdl, robust=TRUE)
  res <- glmQLFTest(fit)
  
  ###save results for plotting 
  top <- topTags(res,n=nrow(y))$table
  top$Gene <- rownames(top)
  write.csv(top,paste0(output_file_path, "Cell_type_abundance_Mets_distance_hmFISH.csv"))
}


###Function to get proportions of cell types per spatial feature area  
create_table_cell_type_prop_resolve <- function(
  seurat_object, 
  ident_of_interest1,
  ident_of_interest2,
  output_file_path,
  sample_name
) {
  ann_tab<-  table(seurat_object@meta.data[[ident_of_interest1]], seurat_object@meta.data[[ident_of_interest2]])
  ann_tab <- cbind(ann_tab, Total = rowSums(ann_tab))
  ann_tab <- as.data.frame(ann_tab)
  ann_tab_pct = lapply(ann_tab[,], function(x) {
    x/ann_tab$Total})
  ann_tab_pct <- as.data.frame(ann_tab_pct)
  rownames(ann_tab_pct) <- rownames(ann_tab)
  ann_tab_pct$sample <- sample_name
  write.csv(ann_tab_pct, file = paste0(output_file_path,sample_name,"_proportions_",ident_of_interest1,"_",ident_of_interest2,".csv"))
}

###Function to generate boxplots per cell type of interest 
#save as pdf or svg
boxplot_cell_prop <- function(
  df,
  cell_type,
  color_oi,
  output_file_path,
  save_as
) {
  p <- ggplot(df, aes(x=X, y=df[[cell_type]], fill=X)) + theme_classic() +
    geom_boxplot(fill = color_oi,outlier.shape = NA) + theme(axis.text = element_text(size = 30))   +
    geom_jitter(color="black", size=5, alpha=0.9) +theme(axis.title= element_text(size = 25)) + 
    ggtitle(paste0(cell_type," proportion - hmFISH")) + xlab("Condition") + 
    ylab("Cell type proportion") +
    theme(plot.title = element_text(size = 25, face = "bold")) + theme(legend.position="none")
  p + ggsave(paste0(output_file_path,cell_type,save_as),width = 7, height = 10)
}


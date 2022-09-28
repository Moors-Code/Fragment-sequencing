###Function to test Cell type abundances between Mets_distance for annotation column of interest 
#following workflow of: http://bioconductor.org/books/3.15/OSCA.multisample/differential-abundance.html#performing-the-da-analysis 
DA_analysis_cell_type_abundance_Mets_distance <- function(
  seurat_object,
  seurat_object_anno, #f.e. seurat_object$annotation.broad
  output_file_path,
  sample_name
) {
  ###load and prepare data 
  #add extra meta data column that combines cell type and group ID of interest 
  seurat_object$sample_group_oi <- paste0(seurat_object$orig.ident,"_",seurat_object$Mets_distance)
  #convert to SCE object 
  sce <- as.SingleCellExperiment(seurat_object)

  cluster.counts <- table(seurat_object_anno,seurat_object$sample_group_oi)
  cluster.counts <- as(cluster.counts, "matrix")
  #add extra information 
  extra.info <- colData(sce)[match(colnames(cluster.counts),sce$sample_group_oi),]
  y.ab <- DGEList(cluster.counts, samples = extra.info)
  
  y.ab$samples$Mets_distance <- factor(y.ab$samples$Mets_distance)
  y.ab$samples$Mets_distance <- ordered(y.ab$samples$Mets_distance)

  ###Defining the model matrix
  #account for variability in samples with '~ factor(orig.ident) + ...'
  mdl <- model.matrix(~factor(orig.ident) + factor(Mets_distance),y.ab$samples)
  
  ###Follow the standard edgeR workflow
  y <- estimateDisp(y.ab, mdl)
  fit <- glmQLFit(y, mdl, robust=TRUE)
  res <- glmQLFTest(fit)
  
  ###save results for plotting 
  top <- topTags(res,n=nrow(y))$table
  top$Gene <- rownames(top)
  write.csv(top,paste0(output_file_path,"Mets_distance","_",sample_name, "_cell_type_abundance_Sphere_seq.csv"))
}

###Function to create table with cell type proportions per sample  
create_table_cell_type_prop <- function(
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
  write.csv(ann_tab_pct, file = paste0(output_file_path,sample_name,"_proportions_",ident_of_interest1,"_",ident_of_interest2,".csv"))
}

###Function to make boxplot of cell type abundance per sample of cell type of interest 
#you can save as pdf or svg 
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
    ggtitle(paste0(cell_type," proportion - sphere-seq")) + xlab("Condition") + 
    ylab("Cell type proportion") + 
    theme(plot.title = element_text(size = 25, face = "bold")) + theme(legend.position="none")
  p + ggsave(paste0(output_file_path,cell_type,save_as),width = 7, height = 10)
}


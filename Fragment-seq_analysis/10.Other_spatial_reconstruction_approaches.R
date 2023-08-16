######################### Part10: Other spatial reconstruction approach #################
#this code applies an alternative spatial reconstruction approach of fragments 

########## Prepare environment ##########
###Setting the working directory 
setwd("~/")

###Load packages and functions 
source("./functions_and_packages/1.Packages.R")

### load R object 
fragment_seq <- readRDS("./Sphere-seq_analysis/data_files_generated/LiverMerged_afterBC_anno_BS_5cells_zC_lobules_Mets_distance.Rda")
MC_data <- readRDS("./Highly_multiplexed_FISH_analysis/data_files_generated/Resolve_seurat_anno.rds")

### only take M5 sample (one metastatic liver)
Idents(fragment_seq) <- "orig.ident"
fragment_seq <- subset(fragment_seq, idents = "M5")

########## find metastatic proximal and distal marker genes in MC data ##########
#remove Fibroblasts, Hepatocytes, and Stellate cells --> not really present in fragment-seq data 
Idents(MC_data) <- "annotation"
MC_data <- subset(MC_data, idents = c("B","Kupffer","LECs","Metastasis","Monocytes","T"))
Idents(MC_data) <- "Mets_distance"
marker_genes <-  FindAllMarkers(object = MC_data, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
View(marker_genes %>% group_by(cluster) %>% top_n(n =10, wt = avg_log2FC))

### plot top 10 DEGs from proximal and distal 
top_10_genes <- (marker_genes %>% group_by(cluster) %>% top_n(n =10, wt = avg_log2FC))$gene

p <- DotPlot(MC_data, features = top_10_genes,dot.scale = 20, scale = FALSE) + 
  theme(legend.title = element_text(size = 20), legend.text = element_text(size = 20)) + 
  theme(title = element_text(size = 20))+ theme(axis.text = element_text(size = 20)) +
  ggtitle("Top5 MC") + theme(axis.text.x = element_text(angle = 90))
p + ggsave("dotPlot_top5_MC.svg", width = 12, height = 10)  

a <- AverageExpression(MC_data,group.by = "Mets_distance", return.seurat = FALSE, 
                  features = top_10_genes)

write.csv(a, "/Users/handler/Desktop/Sup11_a.csv")

########### use these top 20 genes to generate and plot Pseudobulks from fragment-seq data ###############
fragment_seq <- NormalizeData(fragment_seq, normalization.method = "LogNormalize",
                            scale.factor = 10000,
                            margin = 1, assay = "RNA")
fragment_seq_AE <- AverageExpression(fragment_seq,assays = "RNA",slot = "data", group.by = "sphere", return.seurat = TRUE, 
                                   features = top_10_genes)
fragment_seq_AE <- ScaleData(fragment_seq_AE, vars.to.regress = c("nFeature_RNA"), verbose = FALSE)
fragment_seq_AE <- FindVariableFeatures(object = fragment_seq_AE)
fragment_seq_AE <- RunPCA(object = fragment_seq_AE, features = top_10_genes, npcs = 10, verbose = FALSE)
print(ElbowPlot(fragment_seq_AE))
fragment_seq_AE <- FindNeighbors(object = fragment_seq_AE, dims = 1:10)
fragment_seq_AE <- FindClusters(fragment_seq_AE, resolution = 0.4, random.seed = 5, algorithm = 1)
fragment_seq_AE <- RunUMAP(fragment_seq_AE, dims = 1:10, seed.use = 5)
p <- DimPlot(fragment_seq_AE, reduction = "umap", label = TRUE, pt.size = 8, label.size = 10) 
p + ggsave("/Users/handler/Desktop/umap.svg",width = 12, height = 10)

### define marker genes from each cluster 
markers <- FindAllMarkers(fragment_seq_AE,only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
View(markers %>%top_n(n =20, wt = avg_log2FC))

#
p <- FeaturePlot(fragment_seq_AE, features = c("Plxnc1", "Pecam1","Spp1","Sema4d"), pt.size = 5) 
p + ggsave("/Users/handler/Desktop/umap_features.svg",width = 12, height = 10)

### subset fragment barcodes of cluster 2 and cluster 0/1 
sub <- subset(fragment_seq_AE, idents = c(1))
cluster2_bc <- as.data.frame(table(rownames(sub@meta.data)))$Var1
sub2 <- subset(fragment_seq_AE, idents = c(0))
cluster1_bc <- as.data.frame(table(rownames(sub2@meta.data)))$Var1

#ann new identity and analyse metastatic cells 
fragment_seq$Pseudobulk_cluster <- NA
fragment_seq@meta.data <- fragment_seq@meta.data %>%
  mutate(Pseudobulk_cluster = case_when(
    sphere %in% cluster2_bc ~ "cluster2",
    sphere %in% cluster1_bc ~ "cluster1",
    TRUE ~ NA_character_))

#DimPlot
p <- DimPlot(fragment_seq, label = TRUE,label.size = 5, group.by = "annotation.broad", pt.size = 0.5, 
             cols = c("#EF975B","#56595B","#EDB2D4","#91C6C4",
                      "#8E5229","#F46042","#E20FE8",
                      "#C491C6","#F4E740","#F90606","#19E80F",
                      "#AEAEAF" ,"#2C91E5", "#C1A08A","#2323E5"), split.by = "Pseudobulk_cluster") + 
  theme(legend.title = element_text(size = 22), legend.text = element_text(size = 22)) +
  theme(title = element_text(size = 25))+ theme(axis.text = element_text(size = 30)) +
  ggtitle("Broad annotation") 


### plot cell type proportions 
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

create_table_cell_type_prop_for_plot_one_sample <- function(
    df_cell_type_prop_per_sample, 
    sample_name,
    n_cell_types_dim
){
  df <- df_cell_type_prop_per_sample[df_cell_type_prop_per_sample$X %in% sample_name,]
  df <- t(df)
  df <- as.data.frame(df)
  colnames(df) <- NULL
  df <- df[n_cell_types_dim,]
  df$cell_types <- rownames(df)
  colnames(df) <- c("proportion","cell_types")
  rownames(df) <- NULL
  df$sample <- sample_name
  return(df)
}
create_table_cell_type_prop(fragment_seq, "Pseudobulk_cluster","annotation.broad","./","cell_types")
df <- read.csv("./cell_types_proportions_Pseudobulk_cluster_annotation.broad.csv", header = TRUE)

df$Total <- NULL

df1 <- create_table_cell_type_prop_for_plot_one_sample(df,"cluster1",c(2:16)) 
df2 <- create_table_cell_type_prop_for_plot_one_sample(df,"cluster2",c(2:16)) 

df <- rbind(df1,df2)
df <- rbind(df,df2)

df$proportion <- as.numeric(df$proportion)

p <- ggplot(data=df, aes(x=sample, y=proportion, fill = cell_types)) +
  geom_bar(stat="identity", position = "fill" ) + ggtitle("Cell type prop per sample") + 
  theme(axis.text = element_text(size = 5)) + theme(axis.title= element_text(size = 5)) +  
  theme(plot.title = element_text(size = 5, face = "bold")) + 
  theme(legend.title = element_text(size = 5), legend.text = element_text(size = 5))  + 
  xlab("Sample") + ylab("Cell type proportion") + theme(axis.text.x = element_text(angle = 90)) +theme_classic(base_size = 25)  + coord_flip() + 
  scale_fill_manual(values=  c("#EF975B","#56595B","#EDB2D4","#91C6C4",
                               "#8E5229","#F46042","#E20FE8",
                               "#C491C6","#F4E740","#F90606","#19E80F",
                               "#AEAEAF" ,"#2C91E5", "#C1A08A","#2323E5"))
p + ggsave("./cell_type_per_cond.svg", width = 12, height = 6)



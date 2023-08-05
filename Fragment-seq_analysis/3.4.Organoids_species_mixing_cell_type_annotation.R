########## Part 3.4: Organoids species mixing - cell type annotation ##########
#This part clusters and annotates cells as human or mouse

########## Prepare environment ##########
###Setting the working directory 
setwd("/mnt/khandler/R_projects/Fragment-sequencing/Fragment-seq_analysis/")

###Load packages and functions 
source("./functions_and_packages/1.Packages.R")

###Load data 
org_mix_ensemble <- readRDS("./data_files_generated/SpS_Organoid_mixing_ensembleIDs_BS_mito0.3_decontX.Rda")

########## Clustering ##########
org_mix_ensemble <- SCTransform(org_mix_ensemble, vars.to.regress = c("nCount_RNA"), verbose = FALSE, assay = "decontXcounts")
org_mix_ensemble <- RunPCA(object = org_mix_ensemble, features = VariableFeatures(object = org_mix_ensemble), npcs = 30, verbose = FALSE)
ElbowPlot(org_mix_ensemble)
org_mix_ensemble <- FindNeighbors(object = org_mix_ensemble, dims = 1:10, reduction = "pca")
org_mix_ensemble <- FindClusters(org_mix_ensemble, resolution = 0.05, random.seed = 5, algorithm = 4, graph.name = "SCT_nn")
org_mix_ensemble <- RunUMAP(org_mix_ensemble, dims = 1:10, seed.use = 5, reduction = "pca")
DimPlot(org_mix_ensemble, reduction = "umap", label = TRUE, pt.size = .5, group.by = "seurat_clusters")

########## Annotation ##########
dge_genes <- FindAllMarkers(org_mix_ensemble, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
View(dge_genes %>% group_by(cluster) %>% top_n(n =10, wt = avg_log2FC))

#cluster 1 = mouse, cluster 2 = human 
Idents(org_mix_ensemble) <- "seurat_clusters"
current.cluster.ids <- c(1:2)
new.cluster.ids <- c("mouse","human")
org_mix_ensemble$annotation <- plyr::mapvalues(x = org_mix_ensemble$seurat_clusters, from = current.cluster.ids, to = new.cluster.ids)

DimPlot(org_mix_ensemble, group.by = "annotation", label = TRUE)

########## save R object ##########
saveRDS(org_mix_ensemble,"./data_files_generated/SpS_Organoid_mixing_ensembleIDs_BS_mito0.3_decontX_annotated.Rda")

########## Part 2: Merging and batch effect correction  ##########
#This part merges both Visium samples and applies batch effect correction 

########## Prepare environment ##########
###Setting the working directory 
setwd("/mnt/khandler/R_projects/Sphere-sequencing/Visium_analysis/")

###Load packages and functions 
source("./functions_and_packages/1.Packages.R")

###load R object 
V1 <- readRDS("./data_files_generated/Visium_sample1.rds")
V2 <- readRDS("./data_files_generated/Visium_sample2.rds")

########## Merging of samples ##########
merged <- merge(V1, c(V2))
#apply 100 features per spot cutoff, you will loose many spots of the metastatic areas 
merged = merged[, merged$nFeature_Spatial > 100 ]

#normalization and scaling
merged <- SCTransform(merged, assay = "Spatial", verbose = TRUE, method = "poisson")

###Clustering 
merged <- RunPCA(merged, assay = "SCT", verbose = FALSE)
merged <- FindNeighbors(merged, reduction = "pca", dims = 1:30)
merged <- FindClusters(merged, verbose = FALSE)
merged <- RunUMAP(merged, reduction = "pca", dims = 1:30)
DimPlot(merged,group.by = "sample") 

########## Batch effect correction ##########
#create a list of the original data that we loaded to start with
st.list = list(V1,V2)
#run SCT on both datasets
st.list = lapply(st.list, SCTransform, assay = "Spatial", method = "poisson")

#need to set maxSize for PrepSCTIntegration to work
options(future.globals.maxSize = 2000 * 1024^2)  # set allowed size to 2K MiB

st.features = SelectIntegrationFeatures(st.list, nfeatures = 3000, verbose = FALSE)
st.list <- PrepSCTIntegration(object.list = st.list, anchor.features = st.features,
                              verbose = FALSE)

int.anchors <- FindIntegrationAnchors(object.list = st.list, normalization.method = "SCT",
                                      verbose = FALSE, anchor.features = st.features)
merged.integrated <- IntegrateData(anchorset = int.anchors, normalization.method = "SCT",
                                   verbose = FALSE)

merged.integrated <- RunPCA(merged.integrated, verbose = FALSE)
merged.integrated <- FindNeighbors(merged.integrated, dims = 1:10)
merged.integrated <- FindClusters(merged.integrated, verbose = FALSE,resolution =0.1)
merged.integrated <- RunUMAP(merged.integrated, dims = 1:10)

DimPlot(merged.integrated, reduction = "umap", group.by = c("sample")) 

########## save R object ##########
saveRDS(merged.integrated,"./data_files_generated/Visium_merged_batch_corrected.rds")




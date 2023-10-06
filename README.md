# Fragment-sequencing unveils local tissue microenvionments at single cell resolution  

This repository contains code from data analysis of fragment-sequencing split into 3 parts: 

- Fragment-sequencing analysis 
- Analysis of Highly multiplexed FISH (Molecular Cartography) data for validation of findings from fragment-sequencing
- Visium data analysis to compare fragment-sequencing


## Citation 

#### The code in this repository pertains to publication: 

### Fragment-sequencing unveils local tissue microenvionments at single cell resolution  

Kristina Handler, Karsten Bach, Costanza Borrelli, Xenia Ficht, Ilhan E. Acar, Andreas E. Moor


## Data
Data for the above publication can be found at GEO with the access number: GSE192742 and Resolve data at zenodo.org with the doi: 10.5281/zenodo.8413572


## Code description: 


### Fragment-seq_analysis:

#### Functions_and_packages: 

1.Packages.R: Code to load libraries of the required packages. 

2.Functions_preprocessing.R: Code contains functions for preprocessing steps: Demultiplexing of MULTI-seq libraries; Assignment of cells to their fragment-BCs (MULTI-seq BCs); Annotation of zUMI outputs with Gene IDs instead of Ensemble IDs; 

3.Functions_fragment_size_GFP_integration.R: Code contains functions for integrating fragment-seq with large fragment sorter outputs: calculating fragment-size and normalized GFP signals per fragment; plotting of fragment-size per fragment and sample; 

4.Functions_cells_per_fragment_cutoff.R: Code contains function for extracting fragments with an required amount of cells. 

5.Functions_lobule_layer.R: Code contains the following functions for zonated gene expression analysis: assignment of fragments to their lobule layer (L1-L10) and vein of origin (L1-L5 = central vein (CV); L5-L10 = portal vein (PV)) by calculating a zonation coordinate (ZC) per fragment; DGE analysis to find new zonation specific genes in a cell type of interest; plotting of gene expression in fragments across a scaled distance from CV to PV in Liver endothelial cells (LECs) and Kupffer cells (KCs);  

6.Functions_Cell_type_abundance.R: Code contains function for cell type abundance analysis: Differential abundance analysis between two conditions and plotting of cell type proportion of interest between two conditions.  

7.Functions_ligand_receptor_analysis.R: Code contains functions for L-R interaction analysis: Generation of input files for CellPhoneDB analysis; Generation of a matrix that contains interaction scores and p-values from CellPhoneDB comparing two conditions of interacting cell types of interest. 


#### Main analysis: 

##### 1.Pre-processing: 

1.1.Preprocessing_MULTIseq_demux.R: This part does demultiplexing of MULTI-seq FASTQ files to allocate single cells to their fragment of origin.  

1.2.Preprocessing_WTA_MULTIseqBC_integration.R: This part integrates fragment-BC with WTA Seurat objects. 


##### 2.Mouse metastatic liver fragment-seq anaysis: 

2.1.Liver_batch_effect_correction_and_clustering.R: This part merges all samples by applying bath effect correction; then, all samples undergo quality control, normalization, scaling and clustering. 

2.2.Liver_cell_type_annotation.R: This part annotates clustered Seurat object using marker genes from https://www.livercellatlas.org. 

2.3.Liver_fragment_size_integration.R: This part integrates fragment sizes calculated from biosorter outputs with merged Seurat object. 

2.4.Liver_cells_per_fragment_cutoff.R: This part applies a cutoff of a least 5 cells per fragment, all other fragments are removed from further analysis. 

2.5.1.Liver_lobule_layer_classification.R: This part allocates fragments into lobule layers between central and portal veins depending on landmark gene expression in LECs. 

2.5.2.Liver_zonated_gene_expression.R: This part investigates zonation specific gene expression in LECs and KCs. 

2.5.3.1.Liver_ligand_receptor_veins.R: This part first generates input files for CellPhoneDB from CV and PV areas and then compares interactions from CV and PV. 

2.5.3.2.CellPHoneDB_veins.sh: This part executes CellPhoneDB analysis of input files generated in 2.5.3.1.

2.6.1.Liver_metastatic_distance_classification.R: This part groups fragments into proximal and distal areas to metastatic sites. 

2.6.2.Liver_cell_type_abundance_mets_distance.R: This part compares cell type abundances of proximal and distal areas. 

2.6.3.1.Liver_ligand_receptor_mets_distance.R: This part first generates input files for CellPhoneDB from proximal and distal areas and then compares interactions from proximal and distal. 

2.6.3.2.CellPHoneDB_mets_distance.sh: This part executes CellPhoneDB analysis of input files generated in 2.6.3.1.


##### 3.CRC organoid mixing species fragment-seq analysis: 

3.1.Organoids_species_mixing_fragment_size_GFP_integration.R: This part integrates fragment size and GFP signal from biosorter data to data of sorted fragments. 

3.2.Organoids_species_mixing_QC.R: This part applies quality control measures and normalization.

3.3.Organoids_species_mixing_decontX.R: This part applies decontX to remove cell free RNA. 

3.4.Organoids_species_mixing_cell_type_annotation.R: This part clusters and annotates cells as human and mouse. 

3.5.Organoids_species_mixing_analysis_correctly_assigned_cells.R: This part investigates the fraction of correctly and wrongly assigned cell. 


##### Preliminary fragment-seq analysis of spleen and Crohn's disease biopsies: 

4.Spleen_preliminary_analysis.R: This part analyses preliminary spleen fragment-seq data. 

5.Crohn_preliminary_analysis.R: This part analyses preliminary Crohn’s disease fragment-seq data. 

##### Additional Fragment-seq analysis from CRC injected liver samples 

6.Fragment_size_cell_counts_bias.R: This part investigates if there are any biases introduced due to different fragment sizes and minimum cells per fragment cutoffs. 

6.1.CellPHoneDB_veins_bias_fragment_size_counts.sh: This part analyses L-R interactions from input data generated in 6.

7.Biosorter_efficiency_and_fragment_size_per_tissue.R: This part generated plots to show efficiency of the biosorter single fragment sorting and compares distribution of fragment-sizes between different tissues. 

8.Liver_SpS_comp_scRNAseq.R: This part compares fragment-seq to conventional scRNA-seq.

9.Liver_fragment_seq_vs_MC.R: This part compares cell types from fragment-seq to Molecular Cartography data. 

10.Other_spatial_reconstruction_approaches.R: This part applies an alternative reconstruction approach of fragments using marker genes from Molecular Cartography and pseudobulk clustering of fragments. 

### Highly_multiplexed_FISH_analysis (Molecular Cartography)

#### functions_and_packages 

1.Packages.R: Code to load libraries of the required packages. 

2.Functions_ImageJ_connection.R: Code contains functions to connect ImageJ readouts with Seurat analysis: extraction of x and y coordinates to match with spatial features from ImageJ; Generation of text files containing Cell IDs and annotation, this can be projected on images in ImageJ; 

3.Functions_DGE.R: Code contains functions for DGE analysis: between veins and between metastatic distances; Plotting of zonated gene expression in boxplots from KCs and LECs; 

4.Functions_cell_type_abundance.R: Code contains functions for comparing cell type abundance between two groups: differential abundance analysis between proximal and distal and plotting of cell type proportions from both groups in boxplots;  

5.Functions_cell_co_localisation.R: Code contains functions for colocalization analysis; 

#### Main analysis 

1.1.Cell_segmentation_Cellpose.sh: This part does cell segmentation of Molecular Cartography images using the python package Cellpose. 

1.2.Process_Segmentation.sh: This part generates count matrices of genes per single cell. 

1.3.Build_SCE_Object.R: This part produces a single cell experiment object of counts per cell. 

2.Preprocessing_Feature_integration.R: This part does preprocessing of SCE object and integrates x and y coordinates of manually drawn spatial areas (CV, PV, Metastasis).

3.Clustering_annotation.R: This part does clustering and annotation of Molecular Cartography data. 

4.1.Veins_zonated_genes_analysis.R: This part validates findings of zonated genes found from fragment-seq analysis. 

5.1.Mets_distance_DGE_analysis.R: This part does DGE analysis between monocytes from proximal and distal areas. 

5.2.Mets_distance_cell_type_abundance.R: This part validates findings of cell type abundance analysis of fragment-seq between proximal and distal metastatic areas. 

5.3.Mets_distance_cell_co_localisation.R: This part does cell colocalization analysis comapring proximal and distal areas. 



### Visium_analysis 

#### functions_and_packages

1.Packages.R: Code to load libraries of the required packages. 

#### Main analysis: 

1.Clustering_anno_nFeature_comp.R: This part does preprocesing of Visium data, annotation of spatial areas and comparison of number of gene features between different spatial areas with fragment-seq data. 

2.Merging_and_batch_effect_correction.R: This part merged both Visium samples and applies batch effect correction. 

3.Spot_deconvolution.R: This part does deconvolution of spots by integrating data from fragment-seq. 

4.Visium_public_data.R: This part uses public data from https://www.livercellatlas.org (Guilliams et al, 2022) to test newly found zonation specific genes in wild type and NAFLD mouse samples. 

#### System requirements 

R version 4.1.0 (2021-05-18)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 20.04.2 LTS

Matrix products: default
BLAS/LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.8.so

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8   
 [6] LC_MESSAGES=C              LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] RANN_2.6.1                  igraph_1.2.11               xbioc_0.1.19                AnnotationDbi_1.56.2        SCDC_0.0.0.9000            
 [6] Rtsne_0.15                  stringdist_0.9.8            ShortRead_1.52.0            GenomicAlignments_1.30.0    Rsamtools_2.10.0           
[11] Biostrings_2.62.0           XVector_0.32.0              BiocParallel_1.28.3         ggpattern_1.0.1             ggforce_0.3.3              
[16] RColorBrewer_1.1-3          singleCellTK_2.7.1          DelayedArray_0.18.0         plyr_1.8.6                  ggrepel_0.9.1              
[21] cowplot_1.1.1               celda_1.12.0                Matrix_1.3-4                edgeR_3.36.0                limma_3.50.1               
[26] SingleR_1.8.1               batchelor_1.10.0            scran_1.22.1                scuttle_1.4.0               SingleCellExperiment_1.16.0
[31] SummarizedExperiment_1.24.0 Biobase_2.52.0              GenomicRanges_1.44.0        GenomeInfoDb_1.30.1         IRanges_2.26.0             
[36] S4Vectors_0.30.0            BiocGenerics_0.40.0         MatrixGenerics_1.4.0        matrixStats_0.59.0          biomaRt_2.50.3             
[41] SeuratObject_4.0.2          Seurat_4.0.3                forcats_0.5.1               stringr_1.4.0               purrr_0.3.4                
[46] readr_2.1.2                 tidyr_1.2.0                 tibble_3.1.2                ggplot2_3.4.0               tidyverse_1.3.1            
[51] dplyr_1.0.7                 deMULTIplex_1.0.2          

loaded via a namespace (and not attached):
  [1] rappdirs_0.3.3             MCMCprecision_0.4.0        scattermore_0.8            R.methodsS3_1.8.2          pkgmaker_0.32.2.900       
  [6] bit64_4.0.5                irlba_2.3.5                R.utils_2.12.0             hwriter_1.3.2              data.table_1.14.2         
 [11] rpart_4.1-15               KEGGREST_1.34.0            RCurl_1.98-1.6             doParallel_1.0.17          generics_0.1.2            
 [16] ScaledMatrix_1.2.0         RSQLite_2.2.9              combinat_0.0-8             future_1.24.0              bit_4.0.4                 
 [21] tzdb_0.2.0                 spatstat.data_2.1-2        xml2_1.3.3                 lubridate_1.8.0            httpuv_1.6.5              
 [26] assertthat_0.2.1           hms_1.1.1                  promises_1.2.0.1           fansi_1.0.2                progress_1.2.2            
 [31] dbplyr_2.1.1               assertive.files_0.0-2      readxl_1.3.1               DBI_1.1.3                  htmlwidgets_1.5.4         
 [36] spatstat.geom_2.2-0        ellipsis_0.3.2             RSpectra_0.16-0            backports_1.4.1            deldir_1.0-6              
 [41] sparseMatrixStats_1.2.1    vctrs_0.5.1                here_1.0.1                 ROCR_1.0-11                abind_1.4-5               
 [46] cachem_1.0.5               RcppEigen_0.3.3.9.1        withr_2.5.0                GSVAdata_1.30.0            checkmate_2.0.0           
 [51] fastmatrix_0.4-12          sctransform_0.3.3          prettyunits_1.1.1          goftest_1.2-3              svglite_2.1.0             
 [56] cluster_2.1.2              lazyeval_0.2.2             crayon_1.5.0               labeling_0.4.2             pkgconfig_2.0.3           
 [61] tweenr_1.0.2               nlme_3.1-152               rlang_1.0.6                globals_0.14.0             lifecycle_1.0.3           
 [66] miniUI_0.1.1.1             registry_0.5-1             filelock_1.0.2             BiocFileCache_2.2.1        enrichR_3.0               
 [71] modelr_0.1.8               rsvd_1.0.5                 rprojroot_2.0.2            cellranger_1.1.0           polyclip_1.10-0           
 [76] lmtest_0.9-39              Rhdf5lib_1.16.0            zoo_1.8-9                  reprex_2.0.1               pheatmap_1.0.12           
 [81] ggridges_0.5.3             png_0.1-8                  viridisLite_0.4.1          rjson_0.2.21               bitops_1.0-7              
 [86] R.oo_1.25.0                KernSmooth_2.23-20         rhdf5filters_1.6.0         blob_1.2.2                 DelayedMatrixStats_1.12.3 
 [91] parallelly_1.30.0          jpeg_0.1-9                 gridGraphics_0.5-1         ggsignif_0.6.3             beachmat_2.10.0           
 [96] scales_1.2.1               memoise_2.0.1              magrittr_2.0.2             ica_1.0-2                  zlibbioc_1.40.0           
[101] compiler_4.1.0             dqrng_0.3.0                fitdistrplus_1.1-6         cli_3.1.0                  listenv_0.8.0             
[106] patchwork_1.1.1            pbapply_1.5-0              MASS_7.3-54                mgcv_1.8-35                tidyselect_1.1.2          
[111] stringi_1.6.2              BiocSingular_1.10.0        assertive.numbers_0.0-2    locfit_1.5-9.5             latticeExtra_0.6-29       
[116] grid_4.1.0                 tools_4.1.0                future.apply_1.8.1         parallel_4.1.0             rstudioapi_0.13           
[121] bluster_1.4.0              foreach_1.5.2              metapod_1.2.0              gridExtra_2.3              farver_2.1.1              
[126] assertive.types_0.0-3      DropletUtils_1.14.2        BiocManager_1.30.16        digest_0.6.28              shiny_1.7.1               
[131] Rcpp_1.0.7                 broom_0.7.12               later_1.3.0                RcppAnnoy_0.0.19           httr_1.4.2                
[136] L1pack_0.40                assertive.properties_0.0-5 colorspace_2.0-2           rvest_1.0.2                XML_3.99-0.9              
[141] fs_1.5.2                   tensor_1.5                 reticulate_1.24            splines_4.1.0              uwot_0.1.11               
[146] statmod_1.4.36             spatstat.utils_2.2-0       systemfonts_1.0.4          plotly_4.10.0              xtable_1.8-4              
[151] jsonlite_1.7.3             assertive.base_0.0-9       R6_2.5.1                   pillar_1.8.1               htmltools_0.5.2           
[156] mime_0.11                  nnls_1.4                   glue_1.5.0                 fastmap_1.1.0              BiocNeighbors_1.12.0      
[161] codetools_0.2-18           fishpond_2.0.1             utf8_1.2.2                 lattice_0.20-44            spatstat.sparse_2.1-0     
[166] ResidualMatrix_1.4.0       multipanelfigure_2.1.2     curl_4.3.2                 leiden_0.3.9               gtools_3.9.2              
[171] magick_2.7.3               survival_3.2-11            munsell_0.5.0              rhdf5_2.38.0               GenomeInfoDbData_1.2.7    
[176] iterators_1.0.14           HDF5Array_1.22.1           haven_2.4.3                reshape2_1.4.4             gtable_0.3.1              
[181] spatstat.core_2.2-0       

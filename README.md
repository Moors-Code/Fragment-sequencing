# Sphere-sequencing unveils local tissue microenvionments at single cell resolution  

This repository contains code from data analysis of sphere-sequencing split into 3 parts: 

- Sphere-sequencing analysis 
- Analysis of Highly multiplexed FISH (Molecular Cartography) data for validation of findings from sphere-sequencing
- Visium data analysis to compare sphere-sequencing


## Citation 

The code in this repository pertains to publication: 

Sphere-sequencing unveils local tissue microenvionments at single cell resolution  

Kristina Handler, Karsten Bach, Costanza Borrelli, Xenia Ficht, Ilhan E. Acar, Andreas E. Moor


## Data
Data for the above publication can be found at GEO with the access number: … and zenodo.org with the access number: 10.5281/zenodo.7197154.


## Code description: 


### Sphere-seq_analysis:

#### Functions_and_packages: 

1.Packages.R: Code to load libraries of the required packages. 

2.Functions_preprocessing.R: Code contains functions for preprocessing steps: Demultiplexing of MULTI-seq libraries; Assignment of cells to their sphere-BCs (MULTI-seq BCs); Annotation of zUMI outputs with Gene IDs instead of Ensemble IDs; 

3.Functions_sphere_size_GFP_integration.R: Code contains functions for integrating sphere-seq with large fragment sorter outputs: calculating sphere-size and normalized GFP signals per sphere; plotting of sphere-size per sphere and sample; 

4.Functions_cells_per_sphere_cutoff.R: Code contains function for extracting spheres with an required amount of cells. 

5.Functions_lobule_layer.R: Code contains the following functions for zonated gene expression analysis: assignment of spheres to their lobule layer (L1-L10) and vein of origin (L1-L5 = central vein (CV); L5-L10 = portal vein (PV)) by calculating a zonation coordinate (ZC) per sphere; DGE analysis to find new zonation specific genes in a cell type of interest; plotting of gene expression in spheres across a scaled distance from CV to PV in Liver endothelial cells (LECs) and Kupffer cells (KCs);  

6.Functions_Cell_type_abundance.R: Code contains function for cell type abundance analysis: Differential abundance analysis between two conditions and plotting of cell type proportion of interest between two conditions.  

7.Functions_ligand_receptor_analysis.R: Code contains functions for L-R interaction analysis: Generation of input files for CellPhoneDB analysis; Generation of a matrix that contains interaction scores and p-values from CellPhoneDB comparing two conditions of interacting cell types of interest. 


#### Main analysis: 

##### 1.Pre-processing: 

1.1.Preprocessing_MULTIseq_demux.R: This part does demultiplexing of MULTI-seq FASTQ files to allocate single cells to their sphere of origin.  

1.2.Preprocessing_WTA_MULTIseqBC_integration.R: This part integrates sphere-BC with WTA Seurat objects. 


##### 2.Mouse metastatic liver sphere-seq anaysis: 

2.1.Liver_batch_effect_correction_and_clustering.R: This part merges all samples by applying bath effect correction; then, all samples undergo quality control, normalization, scaling and clustering. 

2.2.Liver_cell_type_annotation.R: This part annotates clustered Seurat object using marker genes from https://www.livercellatlas.org. 

2.3.Liver_sphere_size_integration.R: This part integrates sphere sizes calculated from biosorter outputs with merged Seurat object. 

2.4.Liver_cells_per_sphere_cutoff.R: This part applies a cutoff of a least 5 cells per sphere, all other spheres are removed from further analysis. 

2.5.1.Liver_lobule_layer_classification.R: This part allocates spheres into lobule layers between central and portal veins depending on landmark gene expression in LECs. 

2.5.2.Liver_zonated_gene_expression.R: This part investigates zonation specific gene expression in LECs and KCs. 

2.5.3.1.Liver_ligand_receptor_veins.R: This part first generates input files for CellPhoneDB from CV and PV areas and then compares interactions from CV and PV. 

2.5.3.2.CellPHoneDB_veins.sh: This part executes CellPhoneDB analysis of input files generated in 2.5.3.1.

2.6.1.Liver_metastatic_distance_classification.R: This part groups spheres into proximal and distal areas to metastatic sites. 

2.6.2.Liver_cell_type_abundance_mets_distance.R: This part compares cell type abundances of proximal and distal areas. 

2.6.3.1.Liver_ligand_receptor_mets_distance.R: This part first generates input files for CellPhoneDB from proximal and distal areas and then compares interactions from proximal and distal. 

2.6.3.2.CellPHoneDB_mets_distance.sh: This part executes CellPhoneDB analysis of input files generated in 2.6.3.1.


##### 3.CRC organoid mixing species sphere-seq analysis: 

3.1.Organoids_species_mixing_sphere_size_GFP_integration.R: This part integrates sphere size and GFP signal from biosorter data to data of sorted spheres. 

3.2.Organoids_species_mixing_QC.R: This part applies quality control measures and normalization.

3.3.Organoids_species_mixing_decontX.R: This part applies decontX to remove cell free RNA. 

3.4.Organoids_species_mixing_cell_type_annotation.R: This part clusters and annotates cells as human and mouse. 

3.5.Organoids_species_mixing_analysis_correctly_assigned_cells.R: This part investigates the fraction of correctly and wrongly assigned cell. 


##### Preliminary sphere-seq analysis of spleen and Crohn's disease biopsies: 

4.Spleen_preliminary_analysis.R: This part analyses preliminary spleen sphere-seq data. 

5.Crohn_preliminary_analysis.R: This part analyses preliminary Crohn’s disease sphere-seq data. 



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

4.1.Veins_zonated_genes_analysis.R: This part validates findings of zonated genes found from sphere-seq analysis. 

5.1.Mets_distance_DGE_analysis.R: This part does DGE analysis between monocytes from proximal and distal areas. 

5.2.Mets_distance_cell_type_abundance.R: This part validates findings of cell type abundance analysis of sphere-seq between proximal and distal metastatic areas. 

5.3.Mets_distance_cell_co_localisation.R: This part does cell colocalization analysis comapring proximal and distal areas. 



### Visium_analysis 

#### functions_and_packages

1.Packages.R: Code to load libraries of the required packages. 

#### Main analysis: 

1.Clustering_anno_nFeature_comp.R: This part does preprocesing of Visium data, annotation of spatial areas and comparison of number of gene features between different spatial areas with sphere-seq data. 

2.Merging_and_batch_effect_correction.R: This part merged both Visium samples and applies batch effect correction. 

3.Spot_deconvolution.R: This part does deconvolution of spots by integrating data from sphere-seq. 

4.Visium_public_data.R: This part uses public data from https://www.livercellatlas.org (Guilliams et al, 2022) to test newly found zonation specific genes in wild type and NAFLD mouse samples. 


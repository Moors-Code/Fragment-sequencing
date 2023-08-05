#!/bin/bash
#Activate CellPhoneDB  environment 
source /mnt/khandler/cell_phone_db/bin/activate
WD="/mnt/khandler/R_projects/Fragment-sequencing/Fragment-seq_analysis/CellPhoneDB/Mets_distance"
CELLPHONE_DIR="/mnt/khandler/R_projects/Fragment-sequencing/Fragment-seq_analysis/CellPhoneDB/Mets_distance"

#Run cellphonedb analyses
cellphonedb method statistical_analysis ${CELLPHONE_DIR}/proximal_meta.txt ${CELLPHONE_DIR}/proximal_count.txt --output-path=${WD} --project-name="proximal" --counts-data=ensembl

cellphonedb method statistical_analysis ${CELLPHONE_DIR}/distal_meta.txt ${CELLPHONE_DIR}/distal_count.txt --output-path=${WD} --project-name="distal" --counts-data=ensembl

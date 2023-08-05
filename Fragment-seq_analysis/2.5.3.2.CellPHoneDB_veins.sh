#!/bin/bash
#Activate CellPhoneDB  environment 
source /mnt/khandler/cell_phone_db/bin/activate
WD="/mnt/khandler/R_projects/Fragment-sequencing/Fragment-seq_analysis/CellPhoneDB/Vein"
CELLPHONE_DIR="/mnt/khandler/R_projects/Fragment-sequencing/Fragment-seq_analysis/CellPhoneDB/Vein"

#Run cellphonedb analyses
cellphonedb method statistical_analysis ${CELLPHONE_DIR}/cv_meta.txt ${CELLPHONE_DIR}/cv_count.txt --output-path=${WD} --project-name="cv" --counts-data=ensembl

cellphonedb method statistical_analysis ${CELLPHONE_DIR}/pv_meta.txt ${CELLPHONE_DIR}/pv_count.txt --output-path=${WD} --project-name="pv" --counts-data=ensembl

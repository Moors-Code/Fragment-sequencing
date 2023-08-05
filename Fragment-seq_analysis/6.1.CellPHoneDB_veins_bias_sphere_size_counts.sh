#!/bin/bash
##small
#Activate CellPhoneDB  environment 
source /mnt/khandler/cell_phone_db/bin/activate
WD="/mnt/khandler/R_projects/Sphere-sequencing/Sphere-seq_analysis/figures/6/CellPhoneDB/small"
CELLPHONE_DIR="/mnt/khandler/R_projects/Sphere-sequencing/Sphere-seq_analysis/figures/6/CellPhoneDB/small"

#Run cellphonedb analyses
cellphonedb method statistical_analysis ${CELLPHONE_DIR}/cv_meta.txt ${CELLPHONE_DIR}/cv_count.txt --output-path=${WD} --project-name="cv" --counts-data=ensembl
cellphonedb method statistical_analysis ${CELLPHONE_DIR}/pv_meta.txt ${CELLPHONE_DIR}/pv_count.txt --output-path=${WD} --project-name="pv" --counts-data=ensembl

##large
#Activate CellPhoneDB  environment 
source /mnt/khandler/cell_phone_db/bin/activate
WD="/mnt/khandler/R_projects/Sphere-sequencing/Sphere-seq_analysis/figures/6/CellPhoneDB/large"
CELLPHONE_DIR="/mnt/khandler/R_projects/Sphere-sequencing/Sphere-seq_analysis/figures/6/CellPhoneDB/large"

#Run cellphonedb analyses
cellphonedb method statistical_analysis ${CELLPHONE_DIR}/cv_meta.txt ${CELLPHONE_DIR}/cv_count.txt --output-path=${WD} --project-name="cv" --counts-data=ensembl
cellphonedb method statistical_analysis ${CELLPHONE_DIR}/pv_meta.txt ${CELLPHONE_DIR}/pv_count.txt --output-path=${WD} --project-name="pv" --counts-data=ensembl


##C20 
#Activate CellPhoneDB  environment 
source /mnt/khandler/cell_phone_db/bin/activate
WD="/mnt/khandler/R_projects/Sphere-sequencing/Sphere-seq_analysis/figures/6/CellPhoneDB/C20"
CELLPHONE_DIR="/mnt/khandler/R_projects/Sphere-sequencing/Sphere-seq_analysis/figures/6/CellPhoneDB/C20"

#Run cellphonedb analyses
cellphonedb method statistical_analysis ${CELLPHONE_DIR}/cv_meta.txt ${CELLPHONE_DIR}/cv_count.txt --output-path=${WD} --project-name="cv" --counts-data=ensembl
cellphonedb method statistical_analysis ${CELLPHONE_DIR}/pv_meta.txt ${CELLPHONE_DIR}/pv_count.txt --output-path=${WD} --project-name="pv" --counts-data=ensembl

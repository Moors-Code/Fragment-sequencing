#!/bin/bash
#Activate CellPhoneDB  environment 
source /mnt/khandler/cell_phone_db/bin/activate
WD="/.../CellPhoneDB/Mets_distance"
CELLPHONE_DIR="/.../CellPhoneDB/Mets_distance"

#Run cellphonedb analyses
cellphonedb method statistical_analysis ${CELLPHONE_DIR}/cellphonedb_meta_close.txt ${CELLPHONE_DIR}/cellphonedb_count_close.txt --output-path=${WD} --project-name="close" --counts-data=ensembl
cellphonedb method statistical_analysis ${CELLPHONE_DIR}/cellphonedb_meta_close_S4.txt ${CELLPHONE_DIR}/cellphonedb_count_close_S4.txt --output-path=${WD} --project-name="close_S4" --counts-data=ensembl
cellphonedb method statistical_analysis ${CELLPHONE_DIR}/cellphonedb_meta_close_S6.txt ${CELLPHONE_DIR}/cellphonedb_count_close_S6.txt --output-path=${WD} --project-name="close_S6" --counts-data=ensembl
cellphonedb method statistical_analysis ${CELLPHONE_DIR}/cellphonedb_meta_close_S7.txt ${CELLPHONE_DIR}/cellphonedb_count_close_S7.txt --output-path=${WD} --project-name="close_S7" --counts-data=ensembl

cellphonedb method statistical_analysis ${CELLPHONE_DIR}/cellphonedb_meta_far.txt ${CELLPHONE_DIR}/cellphonedb_count_far.txt --output-path=${WD} --project-name="far" --counts-data=ensembl
cellphonedb method statistical_analysis ${CELLPHONE_DIR}/cellphonedb_meta_far_S4.txt ${CELLPHONE_DIR}/cellphonedb_count_far_S4.txt --output-path=${WD} --project-name="far_S4" --counts-data=ensembl
cellphonedb method statistical_analysis ${CELLPHONE_DIR}/cellphonedb_meta_far_S6.txt ${CELLPHONE_DIR}/cellphonedb_count_far_S6.txt --output-path=${WD} --project-name="far_S6" --counts-data=ensembl
cellphonedb method statistical_analysis ${CELLPHONE_DIR}/cellphonedb_meta_far_S7.txt ${CELLPHONE_DIR}/cellphonedb_count_far_S7.txt --output-path=${WD} --project-name="far_S7" --counts-data=ensembl


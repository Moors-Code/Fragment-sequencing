#!/bin/bash
#Activate CellPhoneDB  environment 
source /mnt/khandler/cell_phone_db/bin/activate
WD="/.../CellPhoneDB/Mets_distance"
CELLPHONE_DIR="/.../CellPhoneDB/Mets_distance"

#Run cellphonedb analyses
cellphonedb method statistical_analysis ${CELLPHONE_DIR}/proximal_meta.txt ${CELLPHONE_DIR}/proximal_count.txt --output-path=${WD} --project-name="proximal" --counts-data=ensembl
cellphonedb method statistical_analysis ${CELLPHONE_DIR}/proximal_S4_meta.txt ${CELLPHONE_DIR}/proximal_S4_count.txt --output-path=${WD} --project-name="proximal_S4" --counts-data=ensembl
cellphonedb method statistical_analysis ${CELLPHONE_DIR}/proximal_S6_meta.txt ${CELLPHONE_DIR}/proximal_S6_count.txt --output-path=${WD} --project-name="proximal_S6" --counts-data=ensembl
cellphonedb method statistical_analysis ${CELLPHONE_DIR}/proximal_S7_meta.txt ${CELLPHONE_DIR}/proximal_S7_count.txt --output-path=${WD} --project-name="proximal_S7" --counts-data=ensembl

cellphonedb method statistical_analysis ${CELLPHONE_DIR}/distal_meta.txt ${CELLPHONE_DIR}/distal_count.txt --output-path=${WD} --project-name="distal" --counts-data=ensembl
cellphonedb method statistical_analysis ${CELLPHONE_DIR}/distal_S4_meta.txt ${CELLPHONE_DIR}/distal_S4_count.txt --output-path=${WD} --project-name="distal_S4" --counts-data=ensembl
cellphonedb method statistical_analysis ${CELLPHONE_DIR}/distal_S6_meta.txt ${CELLPHONE_DIR}/distal_S6_count.txt --output-path=${WD} --project-name="distal_S6" --counts-data=ensembl
cellphonedb method statistical_analysis ${CELLPHONE_DIR}/distal_S7_meta.txt ${CELLPHONE_DIR}/distal_S7_count.txt --output-path=${WD} --project-name="distal_S7" --counts-data=ensembl


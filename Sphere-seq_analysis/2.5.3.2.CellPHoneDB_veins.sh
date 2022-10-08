#!/bin/bash
#Activate CellPhoneDB  environment 
source /mnt/khandler/cell_phone_db/bin/activate
WD="/.../CellPhoneDB/Vein"
CELLPHONE_DIR="/.../CellPhoneDB/Vein"

#Run cellphonedb analyses
cellphonedb method statistical_analysis ${CELLPHONE_DIR}/cv_meta.txt ${CELLPHONE_DIR}/cv_count.txt --output-path=${WD} --project-name="cv" --counts-data=ensembl
cellphonedb method statistical_analysis ${CELLPHONE_DIR}/cv_S1_meta.txt ${CELLPHONE_DIR}/cv_S1_count.txt --output-path=${WD} --project-name="cv_S1" --counts-data=ensembl
cellphonedb method statistical_analysis ${CELLPHONE_DIR}/cv_S2_meta.txt ${CELLPHONE_DIR}/cv_S2_count.txt --output-path=${WD} --project-name="cv_S2" --counts-data=ensembl
cellphonedb method statistical_analysis ${CELLPHONE_DIR}/cv_S3_meta.txt ${CELLPHONE_DIR}/cv_S3_count.txt --output-path=${WD} --project-name="cv_S3" --counts-data=ensembl
cellphonedb method statistical_analysis ${CELLPHONE_DIR}/cv_S4_meta.txt ${CELLPHONE_DIR}/cv_S4_count.txt --output-path=${WD} --project-name="cv_S4" --counts-data=ensembl
cellphonedb method statistical_analysis ${CELLPHONE_DIR}/cv_S5_meta.txt ${CELLPHONE_DIR}/cv_S5_count.txt --output-path=${WD} --project-name="cv_S5" --counts-data=ensembl
cellphonedb method statistical_analysis ${CELLPHONE_DIR}/cv_S6_meta.txt ${CELLPHONE_DIR}/cv_S6_count.txt --output-path=${WD} --project-name="cv_S6" --counts-data=ensembl
cellphonedb method statistical_analysis ${CELLPHONE_DIR}/cv_S7_meta.txt ${CELLPHONE_DIR}/cv_S7_count.txt --output-path=${WD} --project-name="cv_S7" --counts-data=ensembl
cellphonedb method statistical_analysis ${CELLPHONE_DIR}/cv_S8_meta.txt ${CELLPHONE_DIR}/cv_S8_count.txt --output-path=${WD} --project-name="cv_S8" --counts-data=ensembl
cellphonedb method statistical_analysis ${CELLPHONE_DIR}/cv_S9_meta.txt ${CELLPHONE_DIR}/cv_S9_count.txt --output-path=${WD} --project-name="cv_S9" --counts-data=ensembl

cellphonedb method statistical_analysis ${CELLPHONE_DIR}/pv_meta.txt ${CELLPHONE_DIR}/pv_count.txt --output-path=${WD} --project-name="pv" --counts-data=ensembl
cellphonedb method statistical_analysis ${CELLPHONE_DIR}/pv_S1_meta.txt ${CELLPHONE_DIR}/pv_S1_count.txt --output-path=${WD} --project-name="pv_S1" --counts-data=ensembl
cellphonedb method statistical_analysis ${CELLPHONE_DIR}/pv_S2_meta.txt ${CELLPHONE_DIR}/pv_S2_count.txt --output-path=${WD} --project-name="pv_S2" --counts-data=ensembl
cellphonedb method statistical_analysis ${CELLPHONE_DIR}/pv_S3_meta.txt ${CELLPHONE_DIR}/pv_S3_count.txt --output-path=${WD} --project-name="pv_S3" --counts-data=ensembl
cellphonedb method statistical_analysis ${CELLPHONE_DIR}/pv_S4_meta.txt ${CELLPHONE_DIR}/pv_S4_count.txt --output-path=${WD} --project-name="pv_S4" --counts-data=ensembl
cellphonedb method statistical_analysis ${CELLPHONE_DIR}/pv_S5_meta.txt ${CELLPHONE_DIR}/pv_S5_count.txt --output-path=${WD} --project-name="pv_S5" --counts-data=ensembl
cellphonedb method statistical_analysis ${CELLPHONE_DIR}/pv_S6_meta.txt ${CELLPHONE_DIR}/pv_S6_count.txt --output-path=${WD} --project-name="pv_S6" --counts-data=ensembl
cellphonedb method statistical_analysis ${CELLPHONE_DIR}/pv_S7_meta.txt ${CELLPHONE_DIR}/pv_S7_count.txt --output-path=${WD} --project-name="pv_S7" --counts-data=ensembl
cellphonedb method statistical_analysis ${CELLPHONE_DIR}/pv_S8_meta.txt ${CELLPHONE_DIR}/pv_S8_count.txt --output-path=${WD} --project-name="pv_S8" --counts-data=ensembl
cellphonedb method statistical_analysis ${CELLPHONE_DIR}/pv_S9_meta.txt ${CELLPHONE_DIR}/pv_S9_count.txt --output-path=${WD} --project-name="pv_S9" --counts-data=ensembl


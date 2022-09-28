#!/bin/bash
#Activate CellPhoneDB  environment 
source /mnt/khandler/cell_phone_db/bin/activate
WD="/.../CellPhoneDB/Vein"
CELLPHONE_DIR="/.../CellPhoneDB/Vein"

#Run cellphonedb analyses
cellphonedb method statistical_analysis ${CELLPHONE_DIR}/cellphonedb_meta_cv.txt ${CELLPHONE_DIR}/cellphonedb_count_cv.txt --output-path=${WD} --project-name="cv" --counts-data=ensembl
cellphonedb method statistical_analysis ${CELLPHONE_DIR}/cellphonedb_meta_cv_S1.txt ${CELLPHONE_DIR}/cellphonedb_count_cv_S1.txt --output-path=${WD} --project-name="cv_S1" --counts-data=ensembl
cellphonedb method statistical_analysis ${CELLPHONE_DIR}/cellphonedb_meta_cv_S2.txt ${CELLPHONE_DIR}/cellphonedb_count_cv_S2.txt --output-path=${WD} --project-name="cv_S2" --counts-data=ensembl
cellphonedb method statistical_analysis ${CELLPHONE_DIR}/cellphonedb_meta_cv_S3.txt ${CELLPHONE_DIR}/cellphonedb_count_cv_S3.txt --output-path=${WD} --project-name="cv_S3" --counts-data=ensembl
cellphonedb method statistical_analysis ${CELLPHONE_DIR}/cellphonedb_meta_cv_S4.txt ${CELLPHONE_DIR}/cellphonedb_count_cv_S4.txt --output-path=${WD} --project-name="cv_S4" --counts-data=ensembl
cellphonedb method statistical_analysis ${CELLPHONE_DIR}/cellphonedb_meta_cv_S5.txt ${CELLPHONE_DIR}/cellphonedb_count_cv_S5.txt --output-path=${WD} --project-name="cv_S5" --counts-data=ensembl
cellphonedb method statistical_analysis ${CELLPHONE_DIR}/cellphonedb_meta_cv_S6.txt ${CELLPHONE_DIR}/cellphonedb_count_cv_S6.txt --output-path=${WD} --project-name="cv_S6" --counts-data=ensembl
cellphonedb method statistical_analysis ${CELLPHONE_DIR}/cellphonedb_meta_cv_S7.txt ${CELLPHONE_DIR}/cellphonedb_count_cv_S7.txt --output-path=${WD} --project-name="cv_S7" --counts-data=ensembl
cellphonedb method statistical_analysis ${CELLPHONE_DIR}/cellphonedb_meta_cv_S8.txt ${CELLPHONE_DIR}/cellphonedb_count_cv_S8.txt --output-path=${WD} --project-name="cv_S8" --counts-data=ensembl
cellphonedb method statistical_analysis ${CELLPHONE_DIR}/cellphonedb_meta_cv_S9.txt ${CELLPHONE_DIR}/cellphonedb_count_cv_S9.txt --output-path=${WD} --project-name="cv_S9" --counts-data=ensembl

cellphonedb method statistical_analysis ${CELLPHONE_DIR}/cellphonedb_meta_pv.txt ${CELLPHONE_DIR}/cellphonedb_count_pv.txt --output-path=${WD} --project-name="pv" --counts-data=ensembl
cellphonedb method statistical_analysis ${CELLPHONE_DIR}/cellphonedb_meta_pv_S1.txt ${CELLPHONE_DIR}/cellphonedb_count_pv_S1.txt --output-path=${WD} --project-name="pv_S1" --counts-data=ensembl
cellphonedb method statistical_analysis ${CELLPHONE_DIR}/cellphonedb_meta_pv_S2.txt ${CELLPHONE_DIR}/cellphonedb_count_pv_S2.txt --output-path=${WD} --project-name="pv_S2" --counts-data=ensembl
cellphonedb method statistical_analysis ${CELLPHONE_DIR}/cellphonedb_meta_pv_S3.txt ${CELLPHONE_DIR}/cellphonedb_count_pv_S3.txt --output-path=${WD} --project-name="pv_S3" --counts-data=ensembl
cellphonedb method statistical_analysis ${CELLPHONE_DIR}/cellphonedb_meta_pv_S4.txt ${CELLPHONE_DIR}/cellphonedb_count_pv_S4.txt --output-path=${WD} --project-name="pv_S4" --counts-data=ensembl
cellphonedb method statistical_analysis ${CELLPHONE_DIR}/cellphonedb_meta_pv_S5.txt ${CELLPHONE_DIR}/cellphonedb_count_pv_S5.txt --output-path=${WD} --project-name="pv_S5" --counts-data=ensembl
cellphonedb method statistical_analysis ${CELLPHONE_DIR}/cellphonedb_meta_pv_S6.txt ${CELLPHONE_DIR}/cellphonedb_count_pv_S6.txt --output-path=${WD} --project-name="pv_S6" --counts-data=ensembl
cellphonedb method statistical_analysis ${CELLPHONE_DIR}/cellphonedb_meta_pv_S7.txt ${CELLPHONE_DIR}/cellphonedb_count_pv_S7.txt --output-path=${WD} --project-name="pv_S7" --counts-data=ensembl
cellphonedb method statistical_analysis ${CELLPHONE_DIR}/cellphonedb_meta_pv_S8.txt ${CELLPHONE_DIR}/cellphonedb_count_pv_S8.txt --output-path=${WD} --project-name="pv_S8" --counts-data=ensembl
cellphonedb method statistical_analysis ${CELLPHONE_DIR}/cellphonedb_meta_pv_S9.txt ${CELLPHONE_DIR}/cellphonedb_count_pv_S9.txt --output-path=${WD} --project-name="pv_S9" --counts-data=ensembl


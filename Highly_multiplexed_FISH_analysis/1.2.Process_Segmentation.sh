#!/bin/bash

python ./functions/0_process_cellpose.py 'data/resolve_out/32830-Slide1_B1-2_DAPI_seg.npy' --counts 'data/resolve_out/32830-Slide1_B1-2_results.txt' --sample 'Slide1_B1-2'
python ./functions/0_process_cellpose.py 'data/resolve_out/32830-Slide1_A1-1_DAPI_seg.npy' --counts 'data/resolve_out/32830-Slide1_A1-1_results.txt' --sample 'Slide1_A1-1'
python ./functions/0_process_cellpose.py 'data/resolve_out/32830-Slide1_A2-1_DAPI_seg.npy' --counts 'data/resolve_out/32830-Slide1_A2-1_results.txt' --sample 'Slide1_A2-1'
python ./functions/0_process_cellpose.py 'data/resolve_out/32830-Slide1_A2-2_DAPI_seg.npy' --counts 'data/resolve_out/32830-Slide1_A2-2_results.txt' --sample 'Slide1_A2-2'
python ./functions/0_process_cellpose.py 'data/resolve_out/32830-Slide1_B1-1_DAPI_seg.npy' --counts 'data/resolve_out/32830-Slide1_B1-1_results.txt' --sample 'Slide1_B1-1'
python ./functions/0_process_cellpose.py 'data/resolve_out/32830-Slide1_B2-1_DAPI_seg.npy' --counts 'data/resolve_out/32830-Slide1_B2-1_results.txt' --sample 'Slide1_B2-1'

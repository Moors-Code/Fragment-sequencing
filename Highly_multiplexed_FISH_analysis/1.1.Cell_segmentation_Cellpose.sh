#!/bin/bash
DIR="data/resolve_out/"
python -m cellpose --dir $DIR --pretrained_model nuclei --save_png --save_outlines --verbose --flow_threshold 0.5 --cellprob_threshold -0.2 --fast_mode

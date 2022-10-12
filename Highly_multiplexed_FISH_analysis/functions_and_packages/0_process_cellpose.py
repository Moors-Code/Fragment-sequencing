# coding: utf-8

# # Cellpose output analysis

import skimage as sk
import numpy as np
import matplotlib.pyplot as plt
import scipy
import cpfunctions
import pandas as pd
import argparse
from cellpose import plot, utils

## Get command line input
parser = argparse.ArgumentParser()
parser.add_argument("file",
                    help="absolute path to segmentation file")
parser.add_argument("-c", "--counts", action="store",
                    help="absolute path to count file")
parser.add_argument("-s", "--sample", action="store",
                    help="sample name")

args = parser.parse_args()

file_input = args.file
count_input = args.counts

# Read in segmentation
dat = np.load(file_input,
             allow_pickle=True).item()
nuclei = dat['masks']

print(f"Cellpose identified {nuclei.max()} nuclei in sample {args.sample}")

# Read in counts
counts = pd.read_csv(count_input,sep='\t', header=None)

# Expand nuclei
nuclei_expanded = sk.segmentation.expand_labels(nuclei,distance=10)

count_matrix = cpfunctions.count_molecules(masks = nuclei_expanded, spots = counts)

print(f"Count matrix computed for {args.sample}")

centroid = cpfunctions.get_seg_info(nuclei_expanded,dat['img'])

print(f"Centroid matrix computed for {args.sample}")

centroid_non_expanded = cpfunctions.get_seg_info(nuclei ,dat['img'])

print(f"Centroid matrix computed for non_expanded nuclei {args.sample}")
# Saving ROIs

print(f"Saving ROIs for {args.sample}")

from tifffile import TiffFile, imwrite
from roifile import ImagejRoi
coords = centroid['Coords'].tolist()
# We have to swap x and y in this case
labels = centroid['Cell'].tolist()
roi_file = "../data_out/nuclei_model/" + args.sample + "_ROIs.zip"
n = len(coords)
for i in range(0,n):
    roi = ImagejRoi.frompoints(points=coords[i],name=labels[i])
    roi.tofile(roi_file)


# Save output
count_file = "../data_out/nuclei_model/" + args.sample + "_counts.csv"
count_matrix.to_csv(count_file)

centroid_file = "../data_out/nuclei_model/" + args.sample + "_centroids.csv"
centroid.to_csv(centroid_file)

centroid_2_file = "../data_out/nuclei_model/" + args.sample + "_centroids_non_expanded.csv"
centroid.to_csv(centroid_2_file)

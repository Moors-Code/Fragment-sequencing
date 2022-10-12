import numpy as np
import scipy
import pandas as pd
import skimage as sk

## Function to count molecules within nuclei
def count_molecules(masks,spots):
    """
    Function to count resolve spots that fall into segmented nuclei
    Requires numpy and pandas
    Arguments:
        masks: output from cellpose 'masks' in the npy file
        spots: XX_results.txt file from resolve with coordinates and gene names
    Returns:
        m: count matrix with spots in rows and cells in columns
    """
    spots = spots.rename(columns = { 0 : 'x', 1: 'y', 2 : 'Tile', 3 : 'Gene', 4: 'tbd'})
    masks = masks.transpose()
    ## Assigning transcript position to cell identity
    coords = np.stack([spots['x'],spots['y']],axis=0)
    it = np.nditer(coords,flags=['external_loop'],order='F')
    cell_id = []
    for x in it:
        cell_id += [masks[x[0],x[1]]]
    cell_id = np.array(cell_id)
    # Now we can create a count matrix
    # I extract them as array assuming that's faster
    n = masks.max() + 1
    genes = spots['Gene']
    genes_uq = genes.unique()
    gncounts = []
    for i in range(1,n):
        if np.any(cell_id==i):
            tmp = genes[cell_id==i].value_counts()
        else:
            tmp = [0 for _ in range(genes_uq.size)]
            tmp = pd.Series(tmp)
            tmp.index = genes_uq
        gncounts += [tmp.rename("Cell_" + str(i))]
    count_matrix = pd.concat(gncounts,axis=1).fillna(0)
    return count_matrix

# Function to retrieve the position of the centroid for each cell
def get_seg_info(masks,image):
    """
    Function to compute centroid and area of segmented objects
    Requires numpy, scipy and pandas
    Arguments:
        masks: output from cellpose 'masks' in the npy file
        image: the DAPI image inside the numpy object from cell pose to extract image features
    Returns:
        seg_info pandas data frame with Cells as rows and various columns representing image features (inclduing coordinates and centroid (X,Y))
    """
## First the centroid of the cell using sparse matrix as this is way faster than np.where on the dense matrix!
#    masks_sparse = scipy.sparse.csr_matrix(masks)
#    nonzeros = scipy.sparse.find(masks_sparse)
    # Note that the masks are usually y*x hence 0=y
#    yindex = nonzeros[0]
#    xindex = nonzeros[1]
#    values = nonzeros[2]
#    x = []
#    y = []
#    area = []
#    n = masks.max()+1
#    for nucleus in range(1,n):
#        x += [xindex[values==nucleus].mean().round(0)]
#        y += [yindex[values==nucleus].mean().round(0)]
#        area += [len(xindex[values==nucleus])]
#    cellnames = ["Cell_" + str(i) for i in np.arange(1,n)]
#    seg_info = pd.DataFrame({"Cell" : cellnames,
#                             "X" : x,
#                             "Y" : y,
#                             "Area" : area})
#    return seg_info
# Update: Changed to using scikit's measure function that's speedier
    properties = sk.measure.regionprops(label_image=masks,intensity_image=image)
    stats = {
        'Area': [p.area for p in properties],
        'X': [p.centroid[1] for p in properties],
        'Y': [p.centroid[0] for p in properties],
        'Mean_Int': [p.mean_intensity for p in properties],
        'Max_Int': [p.max_intensity for p in properties],
        'Min_Int': [p.min_intensity for p in properties],
        'Eccentricity': [p.eccentricity for p in properties],
        'Perimeter': [p.perimeter for p in properties],
        'Cell': ["Cell_" + str(p.label) for p in properties],
        'AspectRatio': [p.major_axis_length / p.minor_axis_length for p in properties],
        'Major_axis': [p.major_axis_length for p in properties],
        'Minor_axis': [p.minor_axis_length for p in properties],
        'Euler_Number': [p.euler_number for p in properties],
        'extent': [p.euler_number for p in properties],
        'orientation': [p.orientation for p in properties],
        'solidity': [p.solidity for p in properties],
        'Coords': [p.coords[:, [1, 0]] for p in properties] # Note how I'm swapping the columns here so that 0 is X and 1 is Y
    }
    seg_info = pd.DataFrame(stats)
    return seg_info

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import linregress

########
# Set constants
########
PREFIX = "membrane"
LONG_PREFIX = f"Total_{PREFIX}"
PRECISION = 2
LIMX = 15
LIMY = 1e-4

###################
# Functions
###################

def fit_decay(x, y, LIMX=LIMX, LIMY=LIMY):
    """Fit a linear decay on log(y) vs x for x >= LIMX and y >= LIMY."""
    mask = (x >= LIMX) & (y >= LIMY)
    x_fit = x[mask]
    y_fit = y[mask]
    if len(x_fit) < 2:
        # Not enough points to fit
        return None
    slope, intercept, _, _, _ = linregress(x_fit, np.log(y_fit))
    return slope, intercept

def block_averaging(vect):
    """Divide vector into 3 blocks, compute decay for each block."""
    n = len(vect)
    block_edges = [0, n//3, 2*n//3, n]
    decays = []
    for i in range(3):
        subvect = vect[block_edges[i]:block_edges[i+1]]
        counts, bins = np.histogram(subvect, bins=np.arange(0, np.max(vect)+1)+0.5)
        x = (bins[:-1] + bins[1:])/2
        y = counts / counts.sum()
        fit = fit_decay(x, y)
        if fit is not None:
            slope, _ = fit
            decays.append(abs(1/slope))
        else:
            decays.append(np.nan)
    return decays

################################
# Packdef constant calculation #
################################

# Initialize results
packdef_constants = pd.DataFrame(
    0, 
    index=["PackDef_cst_global", "PackDef_cst_block1", "PackDef_cst_block2",
           "PackDef_cst_block3", "PackDef_cst_all_blocks", "error_all_blocks"],
    columns=["deep", "shallow", "all"]
)

# Open a text file to save results
with open(f"Res_{PREFIX}.txt", "w") as f:

    for DEFECT in ["deep", "shallow", "all"]:

        # Load data
        filename = f"{LONG_PREFIX}_{DEFECT}.txt"
        packdef_data = pd.read_csv(filename, sep="\s+", header=None, usecols=[1]).values.flatten()

        # Histogram and global fit
        counts, bins = np.histogram(packdef_data, bins=np.arange(0, np.max(packdef_data)+1)+0.5)
        x = (bins[:-1] + bins[1:])/2
        y = counts / counts.sum()
        fit = fit_decay(x, y)
        if fit is not None:
            slope, intercept = fit
            global_inv_decay = abs(1/slope)
        else:
            global_inv_decay = np.nan

        # Block averaging
        FITS_3blocks = block_averaging(packdef_data)
        inv_decay_block = np.nanmean(FITS_3blocks)
        error_inv_decay_block = np.nanstd(FITS_3blocks)

        # Save to dataframe
        packdef_constants.loc["PackDef_cst_global", DEFECT] = global_inv_decay
        packdef_constants.loc["PackDef_cst_block1", DEFECT] = FITS_3blocks[0]
        packdef_constants.loc["PackDef_cst_block2", DEFECT] = FITS_3blocks[1]
        packdef_constants.loc["PackDef_cst_block3", DEFECT] = FITS_3blocks[2]
        packdef_constants.loc["PackDef_cst_all_blocks", DEFECT] = inv_decay_block
        packdef_constants.loc["error_all_blocks", DEFECT] = error_inv_decay_block

        # Print results
        f.write(f"Results on {PREFIX} for {DEFECT} defects\n")
        f.write(f"Total number of defects = {len(packdef_data)}\n")
        f.write(f"Using all data: {round(global_inv_decay, PRECISION)} A^2\n")
        for i, val in enumerate(FITS_3blocks, 1):
            f.write(f"Using block {i}: {round(val, PRECISION)} A^2\n")
        f.write(f"Mean +/- sd on 3 blocks: {round(inv_decay_block, PRECISION)} +/- {round(error_inv_decay_block, PRECISION)} A^2\n\n")

################################
# Plotting
################################

# Colors
colors = ["black", "gray", "darkgray", "lightgray", "red"]

# Barplot of all constants
plt.figure(figsize=(8,6))
M = packdef_constants.iloc[0:5].to_numpy()
bp = plt.bar(np.arange(M.shape[1])-0.2, M[0], width=0.15, color=colors[0], label='Global')
plt.bar(np.arange(M.shape[1])-0.05, M[1], width=0.15, color=colors[1], label='Block1')
plt.bar(np.arange(M.shape[1])+0.05, M[2], width=0.15, color=colors[2], label='Block2')
plt.bar(np.arange(M.shape[1])+0.2, M[3], width=0.15, color=colors[3], label='Block3')
plt.bar(np.arange(M.shape[1])+0.35, M[4], width=0.15, color=colors[4], label='Mean')
plt.ylabel("Defect size constant (A^2)")
plt.xticks(np.arange(M.shape[1]), ["deep", "shallow", "all"])
plt.legend()
plt.title("Packing defect constants")
plt.show()

# Barplot with block-averaging only + errors
yvals = packdef_constants.loc["PackDef_cst_all_blocks"].to_numpy()
err = packdef_constants.loc["error_all_blocks"].to_numpy()
plt.figure(figsize=(6,6))
plt.bar(np.arange(len(yvals)), yvals, yerr=err, capsize=5, color=["blue", "green", "orange"])
plt.xticks(np.arange(len(yvals)), ["deep", "shallow", "all"])
plt.ylabel("Defect size constant (A^2)")
plt.title("Packing defect constants computed by block averaging")
plt.show()

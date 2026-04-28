#!/usr/bin/env python3
import sys
import numpy as np
from pathlib import Path
import pandas as pd
import os
import MDAnalysis as mda
import circle_fit as cf
import quadrilateral_mesh
import util
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from scipy.ndimage import gaussian_filter
import ctleelab_plothelper.plothelpers as ph

# pixi inputs
gro_file = sys.argv[1]
xtc_file = sys.argv[2]
interval = int(sys.argv[3])
system = sys.argv[4]
restriction = int(sys.argv[5])
analysis_folder = sys.argv[6]


os.makedirs(analysis_folder, exist_ok=True)

u = mda.Universe(gro_file, xtc_file)

neighborhood_cutoff = 10
mesh_size = 5

x_range_adj = (
    np.floor_divide(np.min(u.atoms.positions[:, 0] + neighborhood_cutoff), mesh_size)
    * mesh_size,
    np.floor_divide(np.max(u.atoms.positions[:, 0] - neighborhood_cutoff), mesh_size)
    * mesh_size,
)
y_range_adj = (
    np.floor_divide(np.min(u.atoms.positions[:, 1]) + neighborhood_cutoff, mesh_size)
    * mesh_size,
    np.floor_divide(np.max(u.atoms.positions[:, 1]) - neighborhood_cutoff, mesh_size)
    * mesh_size,
)

# # Create bin edges
x_edges = np.arange(x_range_adj[0], x_range_adj[1], mesh_size)
y_edges = np.arange(y_range_adj[0], y_range_adj[1], mesh_size)

x_centers = 0.5 * (x_edges[:-1] + x_edges[1:])
y_centers = 0.5 * (y_edges[:-1] + y_edges[1:])
X_c, Y_c = np.meshgrid(x_centers, y_centers, indexing="ij")

X, Y = np.mgrid[
    x_range_adj[0] : x_range_adj[1] : mesh_size,
    y_range_adj[0] : y_range_adj[1] : mesh_size,
]

num_frames = u.trajectory.n_frames // interval + 1
nx, ny = X.shape
n_quad = (nx) * (ny)
quad_shape = (nx, ny)

# initial calculation of mesh so the number of PO4 beads can be determined for the numpy array
Z, mesh_faces, mesh_vertices, mesh_normals, mesh_centroids, leaflets = (
    quadrilateral_mesh.interpolate_cubic_surface(u, "name PO4", X, Y)
)

upper = leaflets["upper"]
po4 = upper.select_atoms("name PO4").positions
n_po4 = len(po4)

# create empty numpy to store the curvature values
# each index will store the curvature value
curvature_per_frame_upper = np.empty((num_frames, n_quad))
PO4_positions_per_frame_upper = np.empty((num_frames, n_po4, 3))
Z_position_per_frame_upper = np.empty((num_frames, n_quad))
radii_curvature_per_frame_upper = np.empty((num_frames, 2))
po4_position_per_frame = np.empty((num_frames, n_po4))

index_of_frame = 0
for ts in u.trajectory[::interval]:
    # handle breaking of membrane where leaflet's are merged
    try:
        Z, mesh_faces, mesh_vertices, mesh_normals, mesh_centroids, leaflets = (
            quadrilateral_mesh.interpolate_cubic_surface(u, "name PO4", X, Y)
        )
    except (UnboundLocalError, ValueError):
        radii_curvature_per_frame_upper[index_of_frame:] = [np.nan, np.nan]
        print(f"Skipping frame {index_of_frame}: leaflet detection failed")
        # fill the rest of the numpy array with nan
        break

    # apply an gaussian filter to deal with artifically high curvature at the edges of the mesh
    Z_smooth = gaussian_filter(Z["upper"], sigma=3)
    Z_position_per_frame_upper[index_of_frame] = Z["upper"].ravel()
    mean_curvature = util.mean_curvature(Z_smooth / 10000, mesh_size / 10000)
    curvature_per_frame_upper[index_of_frame] = mean_curvature.ravel()

    po4 = upper.select_atoms("name PO4").positions
    PO4_positions_per_frame_upper[index_of_frame] = po4

    # po4 is (N,3) array of positions
    mean_z = np.mean(po4[:, 2])  # average membrane plane
    # compute deviation from mean
    deviation_z = po4[:, 2] - mean_z
    # find whether the tip is positive or negative
    if abs(np.max(deviation_z)) > abs(np.min(deviation_z)):
        # tip is upward
        direction = "up"
        sorted_indices = np.argsort(po4[:, 2])  # sort by z ascending
    else:
        # tip is downward
        direction = "down"
        sorted_indices = np.argsort(-po4[:, 2])  # sort by z descending

    # coordinates
    x_positions = po4[:, 0]
    z_positions = po4[:, 2]

    # identify peak of the buckle
    peak_idx = np.argmax(z_positions)
    peak_z_value = z_positions[peak_idx]
    peak_lipid_x = x_positions[peak_idx]

    # define a region around the peak
    z_threshold = peak_z_value - 40.0  # Å = 4 nm below the peak
    x_threshold_min = peak_lipid_x - restriction  # Å = 4 nm on each side
    x_threshold_max = peak_lipid_x + restriction

    # masks for the region around the peak
    peak_lipids_mask_z = z_positions >= z_threshold
    peak_lipid_mask_x = (x_positions >= x_threshold_min) & (
        x_positions <= x_threshold_max
    )

    # combined mask of x and z dimension
    total_peak_lipids_mask = peak_lipids_mask_z & peak_lipid_mask_x

    # selected lipids only from the peak
    peak_lipids = po4[total_peak_lipids_mask]

    # require enough points to determine the curvature
    if len(peak_lipids) < 10:
        radii_curvature_per_frame_upper[index_of_frame] = [np.nan, np.nan]
        index_of_frame += 1
        continue

    xy = peak_lipids[:, [0, 2]]

    # perfect fit circle
    xc, yc, r, s = cf.prattSVD(xy)
    fig, ax = ph.fixed_size_subplots(1, 1, subwidth=3, subheight=1.1)
    theta_fit = np.linspace(-np.pi, np.pi, 180)

    # shift coordinates
    x_shifted = xy[:, 0] - xc
    y_shifted = xy[:, 1] - yc

    x_fit_shifted = r * np.cos(theta_fit)
    y_fit_shifted = r * np.sin(theta_fit)

    #calculate curvature and radius of the fitted circle
    curvature = float(f"{1/(r/10000):.4g}")  # curvature in um^-1
    radius = float(f"{r/10:.4g}")  # nm radius

    #plot fitted circle and the PO4 beads 
    ax.plot(x_fit_shifted / 10, y_fit_shifted / 10, "b-", label="fitted circle", lw=2)
    ax.scatter(x_shifted / 10, y_shifted / 10, c="#e6cd49", label="data")

    # mark center of circle at (0,0)
    ax.plot([0], [0], mew=1, c="#666666")
    ax.set_xlabel("x (nm)")
    ax.set_ylabel("y (nm)")
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlim(-10, 10)
    ax.set_ylim(0, 10)
    ax.grid(False)
    ax.legend(loc="best", labelspacing=0.1)
    ax.grid()
    ax.set_title(
        f"R: {radius:.4g} nm      κ: {curvature:.4g} " + r"$\mathrm{\mu m^{-1}}$",
        fontsize=8,
    )
    fig.savefig(f"{analysis_folder}/circle-fit-frame{index_of_frame}_{restriction}.png")

    radii_curvature_per_frame_upper[index_of_frame] = [radius, curvature]
    index_of_frame += 1

# return the numpy array for curvature, 
np.save(
    f"{analysis_folder}/{system}-mean-upper-curvature.npy", curvature_per_frame_upper
)

np.save(f"{analysis_folder}/{system}-PO4-position.npy", PO4_positions_per_frame_upper)

np.save(f"{analysis_folder}/{system}-Z-position-upper.npy", Z_position_per_frame_upper)

np.save(
    f"{analysis_folder}/{system}-circle-perfect-fit-upper.npy",
    radii_curvature_per_frame_upper,
)


# compute the mean curvature or z value through each frame
mean_curvature_over_frames = np.nanmean(curvature_per_frame_upper, axis=0)
mean_Z_over_frames = np.nanmean(Z_position_per_frame_upper, axis=0)

# shape the data
mean_curvature_over_frames = mean_curvature_over_frames.reshape(quad_shape)
mean_Z_over_frames = mean_Z_over_frames.reshape(quad_shape)

# calculate the mean down the y-dim
mean_curvature_over_frames_xdim = np.nanmean(mean_curvature_over_frames, axis=1)


##-----graphs-----
X_centered = X - np.nanmean(X)
X_centered_only = X_centered[:, 0]

# create the mean curvature figure
# fig, ax = ph.fixed_size_subplots(1, 1, subwidth=4, subheight=1, rmargin_scale=1.5)
# sc = ax.pcolor(
#     X_centered,
#     Y,
#     mean_curvature_over_frames,
#     cmap="viridis",
#     shading="auto",
#     vmin=0,
#     vmax=250,
# )
# ph.add_fixed_colorbar(
#     sc, ax=ax, pad=0.05, aspect=30, label="Mean curvature 1/um ", extend="both"
# )

# ax.xaxis.set_major_locator(MultipleLocator(25))
# ax.set_xlabel("X (A)")
# ax.set_ylabel("Y (A)")
# fig.savefig(f"{analysis_folder}/{system}-mean-upper-curvature.png")
# plt.close(fig)


# # create figure for mean z in each grid
# fig, ax = ph.fixed_size_subplots(1, 1, subwidth=4, subheight=1, rmargin_scale=1.5)
# sc = ax.pcolor(X_centered, Y, mean_Z_over_frames, cmap="viridis", shading="auto")

# ph.add_fixed_colorbar(
#     sc, ax=ax, pad=0.05, aspect=30, label="height (A) ", extend="both"
# )

# ax.set_aspect("equal", adjustable="datalim")
# ax.xaxis.set_major_locator(MultipleLocator(10))
# ax.set_xlabel("X (A)")
# ax.set_ylabel("Y (A)")
# fig.savefig(f"{analysis_folder}/{system}-mean-Zcoords.png")
# plt.close(fig)


# compress height into the x-dimension
mean_Z_over_frames_xdim = np.mean(mean_Z_over_frames, axis=1)
fig, ax = ph.fixed_size_subplots(1, 1, subwidth=4, subheight=1, rmargin_scale=1.5)
sc = ax.scatter(
    X_centered_only,
    mean_Z_over_frames_xdim,
)
ax.set_xlabel("X (A)")
ax.set_ylabel("Mean height (A)")
ax.set_xlim(min(X_centered_only), max(X_centered_only))
fig.savefig(f"{analysis_folder}/{system}-mean-upper-Z-height-along-xdim.png")
plt.close(fig)

# calculate the curvature from the avg height function
mean_curvature_from_height = util.arc_length_curvature(
    mean_Z_over_frames_xdim / 10000, mesh_size / 10000
)
fig, ax = ph.fixed_size_subplots(1, 1, subwidth=4, subheight=1, rmargin_scale=1.5)
sc = ax.scatter(X_centered_only, mean_curvature_from_height)
ax.xaxis.set_major_locator(MultipleLocator(10))
ax.set_xlabel("X (A)")
ax.set_ylabel("Mean curvature 1/um")
fig.savefig(f"{analysis_folder}/{system}-mean-upper-curvature-from-avgheight.png")
plt.close(fig)


print("restriction:", restriction)
print("gro_file:", gro_file)
print("mean radius nm:", np.nanmean(radii_curvature_per_frame_upper[:, 0]))
print("std radius:", np.nanstd(radii_curvature_per_frame_upper[:, 0]))
print("mean curvature um:", np.nanmean(radii_curvature_per_frame_upper[:, 1]))
print("std of curvature:", np.nanstd((radii_curvature_per_frame_upper[:, 1])))
print("number frames:", len(radii_curvature_per_frame_upper[:, 1]))

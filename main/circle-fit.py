#!/usr/bin/env python3

import argparse
import sys
import numpy as np
from pathlib import Path

import MDAnalysis as mda

import circle_fit as cf
import quadrilateral_mesh
import util

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.ticker import MultipleLocator

import ctleelab_plothelper.plothelpers as ph

u = mda.Universe(sys.argv[1])

neighborhood_cutoff = 10
mesh_size = 20

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

Z, mesh_faces, mesh_vertices, mesh_normals, mesh_centroids, leaflets = (
    quadrilateral_mesh.interpolate_cubic_surface(u, "name PO4", X, Y)
)

mean_curvature = util.mean_curvature(Z["upper"] / 10000, mesh_size / 10000)
fig, ax = ph.fixed_size_subplots(1, 1, subwidth=4, subheight=1, rmargin_scale=1.5)
sc = ax.pcolor(
    X,
    Y,
    mean_curvature,
    cmap="viridis",
    shading="auto",
    vmin=0,
    vmax=150,
)
ph.add_fixed_colorbar(
    sc, ax=ax, pad=0.05, aspect=30, label="Mean curvature 1/um ", extend="both"
)

ax.set_aspect("equal", adjustable="datalim")
ax.xaxis.set_major_locator(MultipleLocator(25))
ax.set_xlabel("X (A)")
ax.set_ylabel("Y (A)")
fig.savefig(f"upper-curvature-tmp.png")

plt.close(fig)


upper = leaflets["upper"]

po4 = upper.select_atoms("name PO4").positions

threshold = 20

tip = po4[np.where(po4[:, 2] > np.max(po4[:, 2]) - threshold)]
# print(po4.shape, tip.shape, np.max(po4[:, 2]) - threshold)
xy = tip[:, [0, 2]]

print(xy.shape)

xc, yc, r, s = cf.prattSVD(xy)

print(xc, yc, r, s, f"{1/(r/10000)} um^-1")

x = xy[:, 0]
y = xy[:, 1]

fig, ax = ph.fixed_size_subplots(1, 1, subwidth=3, subheight=3)
theta_fit = np.linspace(-np.pi, np.pi, 180)
x_fit = xc + r * np.cos(theta_fit)
y_fit = yc + r * np.sin(theta_fit)
ax.plot(x_fit, y_fit, "b-", label="fitted circle", lw=2)
ax.plot([xc], [yc], "bD", mec="y", mew=1)
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_aspect("equal", adjustable="datalim")

# plot data
ax.scatter(x, y, c="red", label="data")
ax.legend(loc="best", labelspacing=0.1)
ax.grid()
ax.set_title("Fit Circle")

# plt.show()
fig.savefig(f"circle-fit-tmp.png")

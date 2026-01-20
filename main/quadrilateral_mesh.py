# libraries for multiprocessing running this code directly
# from pathlib import Path
# from tqdm.contrib.concurrent import process_map

# Cailyns outline for code
# goal create a mesh from the PO4 beads in each leaflets
# determine the min dimension of the system. Since the system is NVT there should be no change in the box size
# determine the size of the bins
# both in x and y direction
# divide the min dimension by the bin size to get the number of bins
# MembraneSpectralAnalysis

# size of bins 1A or .1nm
# step = .5

# path = Path("/home/casakurai/scratch/packmem/")

from scipy.interpolate import griddata
import numpy as np
import MDAnalysis as mda
import leaflet_finder

import pdb


def get_interpolated_z_surface(coordinates, X, Y, ag=None):
    """
    Derive surface from distribution of z coordinates in grid.

    Parameters
    ----------
    coordinates : numpy.ndarray
        Coordinates of AtomGroup. Numpy array of shape=(n_atoms, 3).
    X : numpy.ndarray
    Y : numpy.ndarray

    Returns
    -------
    z_surface: np.ndarray
        Surface derived from set of coordinates in grid of `x_range, y_range` dimensions.
        Returns Numpy array of floats of shape (`n_x_bins`, `n_y_bins`)

    """

    if ag is None:
        raise TypeError("Missing keyword argument 'ag' is required")

    _, _, z_coords = coordinates.T
    xy = np.delete(coordinates, -1, 1)

    box_x, box_y, _, _, _, _ = ag.dimensions
    xy_wrapped = np.vstack(
        (
            xy,
            xy + [box_x, 0],
            xy - [box_x, 0],
            xy + [0, box_y],
            xy - [0, box_y],
            xy + [box_x, box_y],
            xy - [box_x, box_y],
            xy + [box_x, -box_y],
            xy - [box_x, -box_y],
        )
    )

    z_wrapped = np.hstack(
        (
            z_coords,
            z_coords,
            z_coords,
            z_coords,
            z_coords,
            z_coords,
            z_coords,
            z_coords,
            z_coords,
        )
    )

    return griddata(xy_wrapped, z_wrapped, (X, Y), method="cubic")


def interpolate_cubic_surface(ag, selection, X, Y, radius=3.6199):
    leaflets, _ = leaflet_finder.determine_leaflets(ag, selection)

    Z = {}
    Z["upper"] = (
        get_interpolated_z_surface(
            leaflets["upper"].select_atoms(selection).positions,
            X,
            Y,
            ag=leaflets["upper"],
        )
        + radius
    )

    Z["lower"] = (
        get_interpolated_z_surface(
            leaflets["lower"].select_atoms(selection).positions,
            X,
            Y,
            ag=leaflets["lower"],
        )
        - radius
    )

    nx, ny = X.shape
    n_quad = (nx - 1) * (ny - 1)

    # Vertices
    vertices = {
        "upper": np.vstack((X.ravel(), Y.ravel(), Z["upper"].ravel())).T,
        "lower": np.vstack((X.ravel(), Y.ravel(), Z["lower"].ravel())).T,
    }

    def ij2index(i, j, ny):
        return j + i * ny

    faces = np.empty((n_quad, 4), dtype=int)
    quad_index = 0
    for i in range(nx - 1):
        for j in range(ny - 1):
            p0 = ij2index(i, j, ny)
            p1 = ij2index(i + 1, j, ny)
            p2 = ij2index(i + 1, j + 1, ny)
            p3 = ij2index(i, j + 1, ny)

            faces[quad_index] = [p0, p1, p2, p3]
            quad_index += 1

    normals = {}
    normals["upper"] = (
        np.linalg.cross(
            vertices["upper"][faces[:, 1]] - vertices["upper"][faces[:, 0]],
            vertices["upper"][faces[:, 2]] - vertices["upper"][faces[:, 0]],
        )
        + np.linalg.cross(
            vertices["upper"][faces[:, 2]] - vertices["upper"][faces[:, 0]],
            vertices["upper"][faces[:, 3]] - vertices["upper"][faces[:, 0]],
        )
    ) / 2

    normals["lower"] = (
        np.linalg.cross(
            vertices["lower"][faces[:, 2]] - vertices["lower"][faces[:, 0]],
            vertices["lower"][faces[:, 1]] - vertices["lower"][faces[:, 0]],
        )
        + np.linalg.cross(
            vertices["lower"][faces[:, 3]] - vertices["lower"][faces[:, 0]],
            vertices["lower"][faces[:, 2]] - vertices["lower"][faces[:, 0]],
        )
    ) / 2

    normals["upper"] /= np.linalg.norm(normals["upper"], axis=1)[:, np.newaxis]
    normals["lower"] /= np.linalg.norm(normals["lower"], axis=1)[:, np.newaxis]

    centroids = {
        "upper": np.mean(vertices["upper"][faces], axis=1),
        "lower": np.mean(vertices["lower"][faces], axis=1),
    }

    return Z, faces, vertices, normals, centroids, leaflets

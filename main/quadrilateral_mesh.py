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


def interpolate_cubic_surface(ag, selection, step: float = 1, radius=3.6199):
    x_range = (
        np.floor_divide(np.min(ag.atoms.positions[:, 0]), step) * step,
        np.floor_divide(np.max(ag.atoms.positions[:, 0]), step) * step,
    )
    y_range = (
        np.floor_divide(np.min(ag.atoms.positions[:, 1]), step) * step,
        np.floor_divide(np.max(ag.atoms.positions[:, 1]), step) * step,
    )

    # # Create bin edges
    x_edges = np.arange(x_range[0], x_range[1], step)
    y_edges = np.arange(y_range[0], y_range[1], step)

    # # Bin centers for plotting
    x_centers = 0.5 * (x_edges[:-1] + x_edges[1:])
    y_centers = 0.5 * (y_edges[:-1] + y_edges[1:])

    X, Y = np.meshgrid(x_centers, y_centers, indexing="ij")

    P, Q = np.mgrid[
        x_range[0] : x_range[1] : step,
        y_range[0] : y_range[1] : step,
    ]

    leaflets, _ = leaflet_finder.determine_leaflets(ag, selection)

    z_surface = {}
    z_surface["upper"] = (
        get_interpolated_z_surface(
            leaflets["upper"].select_atoms(selection).positions,
            P,
            Q,
            ag=leaflets["upper"],
        )
        + radius
    )

    z_surface["lower"] = (
        get_interpolated_z_surface(
            leaflets["lower"].select_atoms(selection).positions,
            P,
            Q,
            ag=leaflets["lower"],
        )
        - radius
    )

    return (P, Q, z_surface, leaflets, X, Y)


def get_flat_surface(ag, selection, step: float = 1, radius=3.6199):
    x_range = (
        np.floor_divide(np.min(ag.atoms.positions[:, 0]), step) * step,
        np.floor_divide(np.max(ag.atoms.positions[:, 0]), step) * step,
    )
    y_range = (
        np.floor_divide(np.min(ag.atoms.positions[:, 1]), step) * step,
        np.floor_divide(np.max(ag.atoms.positions[:, 1]), step) * step,
    )

    # # Create bin edges
    x_edges = np.arange(x_range[0], x_range[1], step)
    y_edges = np.arange(y_range[0], y_range[1], step)

    # # Bin centers for plotting
    x_centers = 0.5 * (x_edges[:-1] + x_edges[1:])
    y_centers = 0.5 * (y_edges[:-1] + y_edges[1:])

    X, Y = np.meshgrid(x_centers, y_centers, indexing="ij")

    P, Q = np.mgrid[
        x_range[0] : x_range[1] : step,
        y_range[0] : y_range[1] : step,
    ]

    leaflets, _ = leaflet_finder.determine_leaflets(ag, selection)

    flat_z_upper = np.mean(leaflets["upper"].positions[:, 2])
    flat_z_lower = np.mean(leaflets["lower"].positions[:, 2])

    z_surface = {}
    z_surface["upper"] = np.full((P.shape), flat_z_upper + radius)
    z_surface["lower"] = np.full((P.shape), flat_z_lower - radius)

    return (P, Q, z_surface, leaflets, X, Y)

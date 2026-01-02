import pickle
from pathlib import Path
import numpy as np
import MDAnalysis as mda
import csv

import BasicFunctions as bf
import cProfile


# loading data from mc and formating data for analysis
def format_data(X, Y, Zupper, Zlower):
    # remove the extra dimension from Z_upper
    # (X,Y) should be vertices of the quadrilateral

    Zupper_2d = Zupper.squeeze()  # removes the dimension representing the Z column
    Zlower_2d = Zlower.squeeze()  # removes the dimension representing the Z column
    Y_np = np.array(Y)
    X_np = np.array(X)

    return X_np, Y_np, Zupper_2d, Zlower_2d


# determines which residues belong to the upper vs lower leaflet
def residue_leaflet_determine(u, selection):
    import main.leaflet_finder as leaflet_finder
    import MDAnalysis as mda

    leaflet, residue_leaflet = leaflet_finder.determine_leaflets(u, selection)

    return leaflet, residue_leaflet


def vectorized_normal_vector_triangles(X, Y, Zupper, Zlower, normalize=True):
    """
    Takes the parallelogram grids that create the membrane mesh, and splits each one to create two triangle grids.
    Computes the triangle grid centroid, vector of triangular grid, coordinates of each vertex of triangle, edges of each triangle.

    Args:
        X (np.array): x-coordinate of the centroid of the parallelogram grids 
        Y (np.array): y-coordinate of the centroid of the parallelogram grids 
        Zupper (np.array): z-coordinate of the centroid of the parallelogram grids for the upper leaflet 
        Zlower (np.array): z-coordinate of the centroid of the parallelogram grids for the upper lower

    Returns: 
        centroid_triangles_upper (np.array): centroid coordinates for each triangle in the upper leaflet 
        centroid_triangles_lower (np.array): centroid coordinates for each triangle in the lower leaflet 
        normal_vector_triangles_upper (np.array): normal vector coordinates for each triangle in the upper leaflet 
        normal_vector_triangles_lower (np.array): normal vector coordinates for each triangle in the lower leaflet 
        coords_triangles_upper (np.array): vertex coordinates for each triangle in the upper leaflet 
        coords_triangles_lower (np.array): vertex coordinates for each triangle in the lower leaflet 
        edges_triangles_upper (np.array): edge coordinates for each triangle in the upper leaflet 
        edges_triangles_lower (np.array): edge coordinates for each triangle in the lower leaflet 


    """
    import pdb

    nx, ny = X.shape
    n_tri = 2 * (nx - 1) * (ny - 1)  # two triangles per quad per leaflet

    # Vertices
    Vupper = np.vstack((X.ravel(), Y.ravel(), Zupper.ravel())).T
    Vlower = np.vstack((X.ravel(), Y.ravel(), Zlower.ravel())).T

    def ij2index(i, j, ny):
        return j + i * ny

    # Triangle indices
    tri_indices = np.empty((n_tri, 3), dtype=int)

    triangle_index = 0
    for i in range(nx - 1):
        for j in range(ny - 1):
            p0 = ij2index(i, j, ny)
            p1 = ij2index(i + 1, j, ny)
            p2 = ij2index(i + 1, j + 1, ny)
            p3 = ij2index(i, j + 1, ny)

            tri_indices[triangle_index] = [p0, p1, p2]
            tri_indices[triangle_index + 1] = [p0, p2, p3]
            triangle_index += 2

    # Calculate normal vectors
    normal_vector_triangles_upper = np.linalg.cross(
        Vupper[tri_indices[:, 1]] - Vupper[tri_indices[:, 0]],
        Vupper[tri_indices[:, 2]] - Vupper[tri_indices[:, 0]],
    )
    normal_vector_triangles_lower = np.linalg.cross(
        Vlower[tri_indices[:, 1]] - Vlower[tri_indices[:, 0]],
        Vlower[tri_indices[:, 2]] - Vlower[tri_indices[:, 0]],
    )
    if normalize:
        normal_vector_triangles_upper = np.divide(
            normal_vector_triangles_upper,
            np.linalg.norm(normal_vector_triangles_upper, axis=1).reshape(-1, 1),
        )
        normal_vector_triangles_lower = np.divide(
            normal_vector_triangles_lower,
            np.linalg.norm(normal_vector_triangles_lower, axis=1).reshape(-1, 1),
        )

    centroid_triangles_upper = np.mean(Vupper[tri_indices], axis=1)
    centroid_triangles_lower = np.mean(Vlower[tri_indices], axis=1)

    # Pull coordinates of triangle vertices
    coords_triangles_upper = Vupper[tri_indices]
    coords_triangles_lower = Vlower[tri_indices]

    # Edge differences
    edges_triangles_upper = np.hstack(
        (
            Vupper[tri_indices[:, 1]] - Vupper[tri_indices[:, 0]],
            Vupper[tri_indices[:, 2]] - Vupper[tri_indices[:, 1]],
            Vupper[tri_indices[:, 0]] - Vupper[tri_indices[:, 2]],
        )
    ).reshape(-1, 3, 3)
    edges_triangles_lower = np.hstack(
        (
            Vlower[tri_indices[:, 1]] - Vlower[tri_indices[:, 0]],
            Vlower[tri_indices[:, 2]] - Vlower[tri_indices[:, 1]],
            Vlower[tri_indices[:, 0]] - Vlower[tri_indices[:, 2]],
        )
    ).reshape(-1, 3, 3)

    return (
        centroid_triangles_upper,
        centroid_triangles_lower,
        normal_vector_triangles_upper,
        normal_vector_triangles_lower,
        coords_triangles_upper,
        coords_triangles_lower,
        edges_triangles_upper,
        edges_triangles_lower,
    )


# determine which grid gets assigned to which lipid
# residue_leaflet_assignment will tell if the lipid belongs to the uppper or lower leaflet
# uses euclidean distance between gl2 coordinates and centroid of the quadrilateral triangle to determine which grid is the closet to the lipid
# output: grid mapping = gl2 coords with the coordinates of the closest triangle grid to it (should reset with each new residue)
# output: closest_mapping = residue:triangular grid with the shortest ED
def closest_grid(
    gl2_coords,
    empty_res,
    residue_leaflet_assignment,
    centroid_triangles_upper,
    centroid_triangles_lower,
):

    # create a copy of the empty_res this is where to closet grid will be stored
    grid_mapping = empty_res.copy()
    closest_mapping = empty_res.copy()

    # loop through the gl2_coords
    # if the atom is a GL2 then do the following
    for residue in gl2_coords:
        residue_gl2_coords = gl2_coords[residue]
        gl2_ED = {}  ##triangle number: all EDs

        if residue in residue_leaflet_assignment["upper"]:
            for triangle, centroid_coords in centroid_triangles_upper.items():
                centroid_coords = centroid_triangles_upper[triangle]
                difference = np.array(centroid_coords) - np.array(
                    residue_gl2_coords
                )  ##must convert to numpy array to do the subtraction
                distance = np.sqrt(np.sum(difference**2))
                gl2_ED[triangle] = distance
            grid_mapping[residue] = gl2_ED
            min_triangle = min(gl2_ED, key=gl2_ED.get)
            closest_mapping[residue] = min_triangle
        if residue in residue_leaflet_assignment["lower"]:
            for triangle, centroid_coords in centroid_triangles_lower.items():
                centroid_coords = centroid_triangles_upper[triangle]
                difference = np.array(centroid_coords) - np.array(residue_gl2_coords)
                distance = np.sqrt(np.sum(difference**2))
                gl2_ED[triangle] = distance
            grid_mapping[residue] = gl2_ED
            min_triangle = min(gl2_ED, key=gl2_ED.get)
            closest_mapping[residue] = min_triangle

        ##determining the lowest ED and the triangle associated with it
        ##gets the triangle with the smallest ED to the gl2 bead

        ##pair the min_triangle to it's gl2 bead

        # if residue in residue_leaflet_assignment["lower"]:
        #     difference = centroid_triangles_lower - coords
        #     distances = np.sqrt(np.sum(difference**2, axis = 1))
        #     ED_lower = {i: d for i, d in enumerate(distances)}
        #     min_key = min(ED_lower, key = ED_lower.get)
        #     grid_mapping[residue] = min_key

    return closest_mapping

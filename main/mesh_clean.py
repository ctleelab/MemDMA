import pickle
from pathlib import Path
import numpy as np
import MDAnalysis as mda
import csv

import cProfile


# determines which residues belong to the upper vs lower leaflet
def residue_leaflet_determine(u, selection):
    import main.leaflet_finder as leaflet_finder
    import MDAnalysis as mda

    leaflet, residue_leaflet = leaflet_finder.determine_leaflets(u, selection)

    return leaflet, residue_leaflet


def vectorized_normal_vector_triangles(X, Y, z_surface, normalize=True):
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
        normal_upper (np.array): normal vector coordinates for each triangle in the upper leaflet
        normal_lower (np.array): normal vector coordinates for each triangle in the lower leaflet
        vertices_upper (np.array): vertex coordinates for each triangle in the upper leaflet
        vertices_lower (np.array): vertex coordinates for each triangle in the lower leaflet
        edges_triangles_upper (np.array): edge coordinates for each triangle in the upper leaflet
        edges_triangles_lower (np.array): edge coordinates for each triangle in the lower leaflet


    """
    nx, ny = X.shape
    n_tri = 2 * (nx - 1) * (ny - 1)  # two triangles per quad per leaflet

    # Vertices
    vertices_upper = np.vstack((X.ravel(), Y.ravel(), z_surface["upper"].ravel())).T
    vertices_lower = np.vstack((X.ravel(), Y.ravel(), z_surface["lower"].ravel())).T

    def ij2index(i, j, ny):
        return j + i * ny

    # Triangle indices
    faces = np.empty((n_tri, 3), dtype=int)

    triangle_index = 0
    for i in range(nx - 1):
        for j in range(ny - 1):
            p0 = ij2index(i, j, ny)
            p1 = ij2index(i + 1, j, ny)
            p2 = ij2index(i + 1, j + 1, ny)
            p3 = ij2index(i, j + 1, ny)

            faces[triangle_index] = [p0, p1, p2]
            faces[triangle_index + 1] = [p0, p2, p3]
            triangle_index += 2

    # Calculate normal vectors
    normal_vector_triangles_upper = np.linalg.cross(
        vertices_upper[faces[:, 1]] - vertices_upper[faces[:, 0]],
        vertices_upper[faces[:, 2]] - vertices_upper[faces[:, 0]],
    )
    normal_vector_triangles_lower = np.linalg.cross(
        vertices_lower[faces[:, 1]] - vertices_lower[faces[:, 0]],
        vertices_lower[faces[:, 2]] - vertices_lower[faces[:, 0]],
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

    centroid_triangles_upper = np.mean(vertices_upper[faces], axis=1)
    centroid_triangles_lower = np.mean(vertices_lower[faces], axis=1)

    vertices = {"upper": vertices_upper, "lower": vertices_lower}
    centroids = {"upper": centroid_triangles_upper, "lower": centroid_triangles_lower}
    normals = {
        "upper": normal_vector_triangles_upper,
        "lower": normal_vector_triangles_lower,
    }

    return (faces, vertices, centroids, normals)


def vectorized_normal_vector_quads(X, Y, z_surface, normalize=True):
    nx, ny = X.shape
    n_tri = 2 * (nx - 1) * (ny - 1)  # two triangles per quad per leaflet

    # Vertices
    vertices_upper = np.vstack((X.ravel(), Y.ravel(), z_surface["upper"].ravel())).T
    vertices_lower = np.vstack((X.ravel(), Y.ravel(), z_surface["lower"].ravel())).T

    def ij2index(i, j, ny):
        return j + i * ny

    # Triangle indices
    faces = np.empty((n_tri, 3), dtype=int)

    triangle_index = 0
    for i in range(nx - 1):
        for j in range(ny - 1):
            p0 = ij2index(i, j, ny)
            p1 = ij2index(i + 1, j, ny)
            p2 = ij2index(i + 1, j + 1, ny)
            p3 = ij2index(i, j + 1, ny)

            faces[triangle_index] = [p0, p1, p2]
            faces[triangle_index + 1] = [p0, p2, p3]
            triangle_index += 2

    # Calculate normal vectors
    normal_vector_triangles_upper = np.linalg.cross(
        vertices_upper[faces[:, 1]] - vertices_upper[faces[:, 0]],
        vertices_upper[faces[:, 2]] - vertices_upper[faces[:, 0]],
    )
    normal_vector_triangles_lower = np.linalg.cross(
        vertices_lower[faces[:, 1]] - vertices_lower[faces[:, 0]],
        vertices_lower[faces[:, 2]] - vertices_lower[faces[:, 0]],
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

    centroid_triangles_upper = np.mean(vertices_upper[faces], axis=1)
    centroid_triangles_lower = np.mean(vertices_lower[faces], axis=1)

    vertices = {"upper": vertices_upper, "lower": vertices_lower}
    centroids = {"upper": centroid_triangles_upper, "lower": centroid_triangles_lower}
    normals = {
        "upper": normal_vector_triangles_upper,
        "lower": normal_vector_triangles_lower,
    }

    return (faces, vertices, centroids, normals)

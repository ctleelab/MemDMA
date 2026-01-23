from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree

list_types_atoms = [
    "W",
    "NA",
    "CL",
    "NC3",
    "PO4",
    "GL1",
    "GL2",
    "C1A",
    "D2A",
    "C3A",
    "C4A",
    "C1B",
    "D2B",
    "C3B",
    "C4B",
    "C2A",
    "C2B",
]
list_types_atoms_lipid = [
    "NC3",
    "PO4",
    "GL1",
    "GL2",
    "C1A",
    "D2A",
    "C3A",
    "C4A",
    "C1B",
    "D2B",
    "C3B",
    "C4B",
    "C2A",
    "C2B",
]
SIZE = 0.5
lipid = "name NC3 PO4 GL1 GL2 C1A D2A C3A C4A C1B D2B C3B C4B C2A C2B"
membrane_sel = "resname DPPC DOPC DOPA DOPS DPPA DPPS"

# setting paths
script_dir = Path(__file__).resolve().parent
analysis_file_path = Path(script_dir.parent / "analysis")


def parse_radius_aliphatic_dict(filename: str) -> dict:
    property_dict = {}
    for line in open(filename, "r"):
        line = line.split("#")[0].strip()  # Remove comments and whitespace
        if line:
            _tmp = line.split()
            key = f"{_tmp[0]} {_tmp[1]}"
            property_dict[key] = (float(_tmp[2]), _tmp[3])
    return property_dict


def cross(v1, v2):
    return np.array(
        [
            v1[1] * v2[2] - v1[2] * v2[1],
            v1[2] * v2[0] - v1[0] * v2[2],
            v1[0] * v2[1] - v1[1] * v2[0],
        ]
    )


def norm(v):
    return np.sqrt(np.dot(v, v))


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
    n_quad = (nx - 1) * (ny - 1)

    # Vertices
    vertices_upper = np.vstack((X.ravel(), Y.ravel(), z_surface["upper"].ravel())).T
    vertices_lower = np.vstack((X.ravel(), Y.ravel(), z_surface["lower"].ravel())).T

    def ij2index(i, j, ny):
        return j + i * ny

    # Triangle indices
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

    #  faces[triangle_index] = [p0, p1, p2]
    # faces[triangle_index + 1] = [p0, p2, p3]

    # Calculate normal vectors
    normal_vector_triangles_upper_a = np.linalg.cross(
        vertices_upper[faces[:, 1]] - vertices_upper[faces[:, 0]],
        vertices_upper[faces[:, 2]] - vertices_upper[faces[:, 0]],
    )
    normal_vector_triangles_upper_b = np.linalg.cross(
        vertices_upper[faces[:, 2]] - vertices_upper[faces[:, 0]],
        vertices_upper[faces[:, 3]] - vertices_upper[faces[:,]],
    )

    normal_upper = (
        normal_vector_triangles_upper_a + normal_vector_triangles_upper_b
    ) / 2

    print(normal_upper.shape, n_quad)

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

    return faces, vertices, centroids, normals



def mean_curvature(Z, h):
    """
    Calculates mean curvature from Z cloud points.


    Parameters
    ----------
    Z: np.ndarray.
        Multidimensional array of shape (n,n).
    h: float.
        Regular grid separation

    Returns
    -------
    H : np.ndarray.
        The result of mean curvature of Z. Returns multidimensional
        array object with values of mean curvature of shape `(n, n)`.

    """

    Zx, Zy = np.gradient(Z, h)
    Zxx, Zxy = np.gradient(Zx, h)
    _, Zyy = np.gradient(Zy, h)

    H = (1 + Zx**2) * Zyy + (1 + Zy**2) * Zxx - 2 * Zx * Zy * Zxy
    H = -H / (2 * (1 + Zx**2 + Zy**2) ** (1.5))

    return H

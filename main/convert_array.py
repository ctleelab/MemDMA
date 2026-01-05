import numpy as np
import util
import csv


def create_triangle_edge_and_coord_arrays(
    list_of_neighborhood_grids,
    edges_triangles,
    coords_vertex_triangles,
    norm_vector_triangles,
    centroid_coords_triangles: dict,
):
    """
    This function creates a dictionary for all neighborhood triangle edges, centroids, normal vectors for each res id

    Args:
        list_of_neighborhood_grids (dict): A dictionary with residue number (key): (neighborhood_list of the grids numbers)
        edges_triangles (dict): A dictionary containing edge vectors for triangles as a list. {triangle_number: [edge1, edge2, edge3]}
        coords_vertex_triangles (dict): A dictionary containing the coordinates of all vertexes of the triangles {14959: array([[113.63737488,  21.55670166, 282.9483378 ],[114.63858032,  21.55670166, 282.60440123],[114.63858032,  22.55933762, 282.61294111]])}
        norm_vector_triangles (dict): A dictionary with triangle number and the coordinates of the normal vector (28104: array([ 0.29990481, -0.31881054,  0.89912009])
        centroid_coords_triangles (dict): A dictionary with triangle number and the coordinates of the centroid of the triangle

    Returns:
        lipid_ids_all_edges (dict): {res_id: array of all neighborhood triangles edges}
        lipid_ids_all_coords_vertex_array (dict): {res_id: array of all neighborhood centroids triangles}
        lipid_ids_all_norms (dict): {res_id: array of all neighborhood triangles normal vectors}
        lipid_ids_all_centroid_coords (dict): {res_id: array of all neighborhood triangles centroid coordinates}
        lipid_ids_all_grid_neighborhood (dict): {res_id: [list of neighborhood grids]}

    """

    lipid_ids_all_edges_array = {}
    lipid_ids_all_coords_vertex_array = {}
    lipid_ids_all_norms_array = {}
    lipid_ids_all_centroid_coords_array = {}
    lipid_ids_all_grid_neighborhood = {}

    for lipid_id, grid_id_list in list_of_neighborhood_grids.items():

        all_edges = []
        all_vertex_coords = []
        all_norms = []
        all_centroids_coords = []
        all_neighborhood_grids = []

        for grid_id in grid_id_list:

            # Extract 3 edge vectors → (3,3) array
            edges = np.array(edges_triangles[grid_id])

            # Extract triangle vertex coordinates → already (3,3)
            coords = coords_vertex_triangles[grid_id]

            # Extract normal vector of the triangle
            norms = norm_vector_triangles[grid_id]

            # Extract the centroid coordinates of the triangles
            centroids_coords = centroid_coords_triangles[grid_id]

            all_edges.append(edges)
            all_vertex_coords.append(coords)
            all_norms.append(norms)
            all_centroids_coords.append(centroids_coords)
            all_neighborhood_grids.append(grid_id)

        # Convert lists to numpy arrays
        all_edges_array = np.array(all_edges)  # shape (N_triangles, 3, 3)
        all_coords_array = np.array(all_vertex_coords)  # shape (N_triangles, 3, 3)
        all_norms_array = np.array(all_norms)  # shape (N_triangles, 3)
        all_centroids_coord_array = np.array(
            all_centroids_coords
        )  # shape (N_triangles, 3)

        # store with a specific residue id number
        lipid_ids_all_edges_array[lipid_id] = all_edges_array
        lipid_ids_all_coords_vertex_array[lipid_id] = all_coords_array
        lipid_ids_all_norms_array[lipid_id] = all_norms_array
        lipid_ids_all_centroid_coords_array[lipid_id] = all_centroids_coord_array
        lipid_ids_all_grid_neighborhood[lipid_id] = all_neighborhood_grids
    return (
        lipid_ids_all_edges_array,
        lipid_ids_all_coords_vertex_array,
        lipid_ids_all_norms_array,
        lipid_ids_all_centroid_coords_array,
        lipid_ids_all_grid_neighborhood,
    )

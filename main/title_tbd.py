import numpy as np
import csv
import util

def create_triangle_edge_and_coord_arrays(list_of_neighborhood_grids, edges_triangles, coords_triangles, norm_vector_triangles, centroid_coords_triangles:dict):
    '''
    This function create a numpy array for edges and coords of triangles for each lipid's neighborhood grids

    Args:
        list_of_neighborhood_grids (dict): A dictionary with lipid (key): (neighborhood_list of the grid_idex)
        edges_triangles (dict): A dictionary containing edge vectors for triangles as a list. {triangle_number: [edge1, edge2, edge3]}
        coords_triangles (dict): A dictionary containing the coordinates of all vertexes of the triangles {14959: array([[113.63737488,  21.55670166, 282.9483378 ],[114.63858032,  21.55670166, 282.60440123],[114.63858032,  22.55933762, 282.61294111]])}
        norm_vector_triangles (dict): A dictionary with triangle number and the coordinates of the normal vector (28104: array([ 0.29990481, -0.31881054,  0.89912009])
        centroid_coords_triangles (dict): A dictionary with triangle number and the coordinates of the centroid of the triangle
    
    Returns:
        lipid_ids_all_edges_array (dict): {res_id: array of all neighborhood triangles edges}
        lipid_ids_all_coords_array (dict): {res_id: array of all neighborhood vertex of the triangles}
        lipid_ids_all_norms_array (dict): {res_id: array of all neighborhood triangles normal vectors}
        lipid_ids_all_centroid_coords_array (dict): {res_id: array of all neighborhood triangles centroid coordinates}
        lipid_ids_all_grid_neighborhood (dict): {res_id: [list of neighborhood grids]}

    '''


    lipid_ids_all_edges_array = {}
    lipid_ids_all_coords_array = {}
    lipid_ids_all_norms_array = {}
    lipid_ids_all_centroid_coords_array = {}
    lipid_ids_all_grid_neighborhood = {}
    lipid_ids_all_centroid_coords_array_dict = {}

    #looping through the lipid_ids and pulling out the neighborhood of grids
    for lipid_id, grid_id_list in list_of_neighborhood_grids.items():
        sorted_grid_ids = sorted(grid_id_list) #by adding this sorting function the order of the grid list matches the index of the numpy array

        all_edges = []
        all_coords = []
        all_norms = []
        all_centroids_coords = []
        all_neighborhood_grids = []
        grid_id_dict_and_centroids = {}

        for grid_id in sorted_grid_ids:

            # Extract 3 edge vectors → (3,3) array
            edges = np.array(edges_triangles[grid_id])

            # Extract triangle vertex coordinates → already (3,3)
            coords = coords_triangles[grid_id]             

            # Extract normal vector of the triangle
            norms = norm_vector_triangles[grid_id] 

            # Extract the centroid coordinates of the triangles
            centroids_coords = centroid_coords_triangles[grid_id]

            

            all_edges.append(edges)

            #create a dict of grid id: centroids
            grid_id_dict_and_centroids[grid_id] = centroids_coords

            all_coords.append(coords)
            all_norms.append(norms)
            all_centroids_coords.append(centroids_coords)
            all_neighborhood_grids.append(grid_id)



        # Convert lists to numpy arrays
        all_edges_array = np.array(all_edges)     # shape (N_triangles, 3, 3)
        all_coords_array = np.array(all_coords)   # shape (N_triangles, 3, 3)
        all_norms_array = np.array(all_norms)     # shape (N_triangles, 3)
        all_centroids_coord_array = np.array(all_centroids_coords)  # shape (N_triangles, 3)

        lipid_ids_all_edges_array[lipid_id] = all_edges_array
        lipid_ids_all_coords_array[lipid_id] = all_coords_array
        lipid_ids_all_norms_array[lipid_id] = all_norms_array
        lipid_ids_all_centroid_coords_array[lipid_id] = all_centroids_coord_array
        lipid_ids_all_grid_neighborhood[lipid_id] = all_neighborhood_grids
        lipid_ids_all_centroid_coords_array_dict[lipid_id] = grid_id_dict_and_centroids




    return  lipid_ids_all_edges_array, lipid_ids_all_coords_array, lipid_ids_all_norms_array, lipid_ids_all_centroid_coords_array, lipid_ids_all_grid_neighborhood
import pickle
import numpy as np
import csv
import util
import neighborhood_grid_graph as ngg
import BasicFunctions as bfg
def sum_atom_mapping_upper_np_array(atom_features:dict, radius:dict,aliphatic:dict, 
                                    lipid_ids_all_grid_edges:dict, lipid_ids_all_grid_vertices_coords:dict, lipid_ids_all_norms_upper: dict, lipid_ids_all_centroid_coords_upper:dict, lipid_grid_assignment:dict, 
                                    atom_neighborhood_grid_id:dict, grid_num:list, frame:int, lipid_type:str):
    '''
    For each atom calculate the mapping value for all grids in a neighborhood at once using numpy arrays 

    Args:
        atom_features (dict): A dictionary containing features of atoms, including their coordinates and type.output ex: {atm#: GL2', 'lower', array([18.222,  4.37 , 14.85 ]), 'DPPC GL2', resnumber }
        radius (dict): A dictionary mapping atom names to their respective radii.
        aliphatic (dict): A dictionary mapping atom names to their type (e.g., polar or non-polar).
        
        lipid_ids_all_grid_edges (dict): {res_id: array of all neighborhood triangles edges represented as a single vector}
        lipid_ids_all_grid_vertices_coords (dict): {res_id: array of all neighborhood vertex of the triangles}
        lipid_ids_all_norms_upper (dict): {res_id: array of all neighborhood triangles normal vectors}
        lipid_ids_all_centroid_coords_upper(dict): A dictionary {res_id: array of all neighborhood triangles centroid coordinates}
        lipid_grid_assignment (dict): A dictionary with lipid residue number (key): (neighborhood_list of the grid_idex)
        
        atom_neighborhood_grid_id (dict): A dictionary containing the res_number the atom belongs to, coords of atom, list of close by grids for each atom number
        grid_num (list): a list of all the grids numbers in either the upper or lower leaflet
        frame (int): the frame number of the trajectory file
        lipid_type: tells you the type of lipids in the simulation

    Returns:
        dict: A dictionary mapping triangle IDs to their respective atom assignment values.
    '''
    #first loop over each atom
    #pull out the residue number it belongs to 
    triangle_atom_assignment_upper = {}

    for atom_num, data in atom_neighborhood_grid_id.items():
        res_num = data[0] #pulls out res_num
        grid_ids = np.asarray(lipid_grid_assignment[res_num], dtype=int)

        grid_to_row = {gid: i for i, gid in enumerate(grid_ids)}


        atom_num = int(
            atom_num
        )  # changes atom_num to a int so I can apply the if statement
       
        if atom_num > 0: #change if want to troubleshoot so only do a set number of atoms
            atom_num = str(
                atom_num
            )  # changes atom_num to a string so it matches atom_features keys


            #-----determine various characteristics of the atom-----
            lipid_atm_name = atom_features[atom_num][3] #determine the resdiue number of the atom
            atom_coords = atom_features[atom_num][2] #gets the coordianates of the atom
            atm_radius = [radius[lipid_atm_name]] #gets the radius of the atom 
            atm_type = aliphatic[lipid_atm_name] #gets the atom type (polar vs non-polar)

            #deals with an atom number if it does not exist in the atom features dictionary
            if atom_num not in atom_features:
                print(f"Warning: atom_num {atom_num} missing in atom_features")
                # fails.append(int(atom_num))
                break

            #Handeling of only atoms that are in the upper leaflet
            if (
                atom_features[atom_num][1] == "lower"
            ):  ##skips atoms that belong to the lower leaflet
                continue

            else:
                
                #pulls the array dependent on the lipid residue number the atom belongs to
                neighborhood_edges_array = lipid_ids_all_grid_edges[res_num]
                neighborhood_grid_vertices_coords_array = lipid_ids_all_grid_vertices_coords[res_num]
                neighborhood_norms_array = lipid_ids_all_norms_upper[res_num]
                neighborhood_centroid_coords_array = lipid_ids_all_centroid_coords_upper[res_num]

                #-----Project the atom onto the same plane as the triangle-----

                # Calculate vector from triangle centroids to atom
                w_vectors =  atom_coords - neighborhood_centroid_coords_array

                # Calculate the normal direction of w_vector by computing the dot product the w_vector to the triangle normal vector
                #sum across rows (sum of x,y,z components) to get the scalar projection of w onto the normal vector 
                dot_products = np.sum(w_vectors * neighborhood_norms_array, axis=1) 
                
                #Computes the normal component of the w_vector
                #dot_products[:, None] changes from row of number to column of numbers to do the multiplication
                #(3,1) * (3,3)
                n_contribution = dot_products[:, None] * neighborhood_norms_array

                #Take the atom and move down the normal component of w to get the projection onto the plane of the grid 
                atom_coords_projections = atom_coords - n_contribution

                #graph of projected atom coordinates, grid centroid, and the vertices of the triangle 
                #ngg.projected_atom_grid(atom_coords_projections, neighborhood_centroid_coords_array,neighborhood_grid_vertices_coords_array, atom_coords, neighborhood_norms_array, atom_num, atm_radius) 
                
                
                #-----Determine if the atom lies within the triangle-----

                #graph of determining if the atom lies within the triangle
                #ngg.atom_projection_grid_assignment(atom_coords_projections, neighborhood_grid_vertices_coords_array, vertex_to_atom) 

                #Computes the vector from each vertex of the traingle to the atom projection
                vector_from_vetex_to_atom_projection = atom_coords_projections[:, None, :] - neighborhood_grid_vertices_coords_array

                #Calculate the euclidian distance from each vertex of the triangle to the atom projection
                atom_coords_projections_transpose = atom_coords_projections[:, None, :] # reshape to (N,1,3) (number of projections, 1, x/y/z)
                euclidian_distance_vertex_to_atom = bfg.Edist3D(atom_coords_projections_transpose, neighborhood_grid_vertices_coords_array)

                #dot product between each edge of the triangle and the vector from the vertex to the atom projection
                #shape (N,3) (number of triangles, 3 vertices)
                dot_product_edge_vertex_to_atom_projection = np.sum(neighborhood_edges_array * vector_from_vetex_to_atom_projection, axis=2) #shape (N,3) number of triangles, 3 edges

                #compute the magnitude of each edge vector of the triangle
                edge_magnitude = np.linalg.norm(neighborhood_edges_array, axis=2)

                for i, grid_id in enumerate(grid_ids):
                    #1) first determine if atom project is within the triangle by looking at the vertices
                    #2) if not, look at the special cases

                    # Convert grid ID to row index
                    row = grid_to_row[grid_id]
                    if (grid_id not in triangle_atom_assignment_upper):  # only initialize the grid if the grid doesn't exist yet
                        triangle_atom_assignment_upper[grid_id] = 0

                    #-----deal with the vertices cases first-----
                    atom_maps_to_triangle = False #first set the atom not mapping to the triangle, which will hold true if none of the conditions are met

                    #gets all distances from atom to all vertices of a given triangle
                    shortest_distance_to_atom_all_vertex = euclidian_distance_vertex_to_atom[row] 
                    
                    # Get the minimum distance from atom to all vertices of the triangle
                    min_distance_to_vertices = np.min(shortest_distance_to_atom_all_vertex)   

                    #handles the case where the atom projection is within radius to a vertex of the triangle
                    if min_distance_to_vertices <= atm_radius:
                        atom_maps_to_triangle = True
                    
                    #-----deals with cases where projected atom is closer to edge of triangle rather than vertices of triangle----

                    else:
                        #pulls out the dot product values for the vertexs for a given triangle
                        tri_vertex = dot_product_edge_vertex_to_atom_projection[row]

                        for tri_vertex_indx, dotproduct in enumerate(tri_vertex):
                            if dotproduct > 0:

                                #pull out the magnitude of the edge vector
                                edge_magnitude_tri_vertex_indx = edge_magnitude[row][tri_vertex_indx]

                                #atom should be within the length of the edge vector
                                if dotproduct <= edge_magnitude_tri_vertex_indx:

                                    #-----compute the area of the parallelogram and divide by the base to get the hieght-----

                                    #compute the area of the parallelogram formed by the edge vector and the atom projection vector 
                                    area_parallelogram = np.linalg.norm(np.cross(vector_from_vetex_to_atom_projection[row][tri_vertex_indx], neighborhood_edges_array[row][tri_vertex_indx]), axis=0)
                                    
                                    #computes the magintude of the edge vector, which is the base of the paralleogram
                                    base_parallelogram = np.linalg.norm(neighborhood_edges_array[row][tri_vertex_indx]) 
                                    
                                    height_parallelogram = area_parallelogram / base_parallelogram

                                    if height_parallelogram <= atm_radius:
                                        atom_maps_to_triangle = True

                    if atom_maps_to_triangle:
                        if (atm_type == "n"):  
                            triangle_atom_assignment_upper[grid_id] += 1.0
        
                        else:
                            triangle_atom_assignment_upper[grid_id] += 0.001


    for grid in grid_num:
        grid_int = int(grid)  # normalize key type
        if grid_int not in triangle_atom_assignment_upper:
            triangle_atom_assignment_upper[grid_int] = 0.0

    triangle_atom_assignment_upper_sorted = dict(
        sorted(triangle_atom_assignment_upper.items())
    )
    
    with open(
        util.trajectory_path
        / f"test_{lipid_type}_triangle_atom_assignment_upper_strain2_frame{frame}.pickle",
        "wb",
    ) as f:
        pickle.dump(triangle_atom_assignment_upper_sorted, f)

    # Create a CSV file to store the results
    with open(
        util.trajectory_path
        / f"test_{lipid_type}_triangle_atom_assignment_upper_strain2_frame{frame}.csv",
        "w",
        newline="",
    ) as f:
        writer = csv.writer(f)
        writer.writerow(["triangle_id", "assignment"])  # optional header
        for key, value in triangle_atom_assignment_upper_sorted.items():
            writer.writerow([key, value])

    # Count number of grids with assignment < 1
    num_grids_less_than_1 = sum(1 for value in triangle_atom_assignment_upper_sorted.values() if value < 1)

    print(f"Number of grids with assignment < 1, atom_num{atom_num} atom_type{atm_type}: {num_grids_less_than_1}")

    return triangle_atom_assignment_upper_sorted

import pickle
import numpy as np
import csv
import util
import neighborhood_grid_graph as ngg
import BasicFunctions as bfg

def sum_atom_mapping_upper_np_array(atom_features:dict, radius:dict,aliphatic:dict, 
                                    lipid_ids_all_grid_edges:dict, lipid_ids_all_grid_vertices_coords:dict, lipid_ids_all_norms_upper: dict, lipid_ids_all_centroid_coords_upper:dict, lipid_grid_assignment:dict,
                                    atom_neighborhood_grid_id:dict, grid_num:list, frame:int):
    '''
    For each atom calculate the mapping value for all grids in a neighborhood at once using numpy arrays 

    Args:
        atom_features (dict): A dictionary containing features of atoms, including their coordinates and type.output ex: {atm#: GL2', 'lower', array([18.222,  4.37 , 14.85 ]), 'DOPC GL2', resnumber }
        radius (dict): A dictionary mapping atom names to their respective radii.
        aliphatic (dict): A dictionary mapping atom names to their type (e.g., polar or non-polar).
        
        lipid_ids_all_grid_edges (dict): {res_id: array of all neighborhood triangles edges}
        lipid_ids_all_grid_vertices_coords (dict): {res_id: array of all neighborhood vertex of the triangles}
        lipid_ids_all_norms_upper (dict): {res_id: array of all neighborhood triangles normal vectors}
        lipid_ids_all_centroid_coords_upper(dict): A dictionary {res_id: array of all neighborhood triangles centroid coordinates}
        lipid_grid_assignment (dict): A dictionary with lipid residue number (key): (neighborhood_list of the grid_idex)

        lipid_ids_all_edges_dict (dict): {res_id: array of all neighborhood triangles edges {grid id: edges}}
        
        atom_neighborhood_grid_id (dict): A dictionary containing the res_number the atom belongs to, coords of atom, list of close by grids for each atom number
        grid_num (list): a list of all the grids numbers in either the upper or lower leaflet
        frame (int): the frame number of the trajectory file

    Returns:
        dict: A dictionary mapping triangle IDs to their respective atom assignment values.
    '''
    #first loop over each atom
    #pull out the residue number it belongs to 
    triangle_atom_assignment_upper = {}

    for atom_num, data in atom_neighborhood_grid_id.items():
        res_num = data[0] #will need this to search the array
        grid_ids_unsorted = np.array(lipid_grid_assignment[res_num], dtype=int) #list of grids in the neighborhood for each resnum

        order = np.argsort(grid_ids_unsorted)
        grid_ids_sorted = grid_ids_unsorted[order]

        grid_to_row = {gid: i for i, gid in enumerate(grid_ids_sorted)} 

        atom_num = int(
            atom_num
        )  # changes atom_num to a int so I can apply the if statement
       
        if atom_num > 0: #change if want to troubleshoot so only do a set number of atoms
            atom_num = str(
                atom_num
            )  # changes atom_num to a string so it matches atom_features keys


            # determine various characteristics of the atom 
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
                #pulls the array dependent on the lipid residue number the atom belongs to
                neighborhood_edges_array = lipid_ids_all_grid_edges[res_num][order]
                neighborhood_grid_vertices_coords_array = lipid_ids_all_grid_vertices_coords[res_num][order]
                neighborhood_norms_array = lipid_ids_all_norms_upper[res_num][order]
                neighborhood_centroid_coords_array = lipid_ids_all_centroid_coords_upper[res_num][order]

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

                #Issue the edge of the triangles are infinite and not restricted to the length of the triangle edge
                #create a new variable that states hit or not, then evaluate it at the end?
                
                #Determine the distance from atom projection to each vertex of the triangle
                #If the distance shortest than radius of the atom, it is considered a hit 
                #Input: atom_coords_projections (N,3), neighborhood_grid_vertices_coords_array (N,3,3) (number of triangles, 3 vertices, x/y/z)


                #graph of determining if the atom lies within the triangle
                ngg.atom_projection_grid_assignment(atom_coords_projections, neighborhood_grid_vertices_coords_array, vertex_to_atom) 

                atom_coords_projections_transpose = atom_coords_projections[:, None, :] # reshape to (N,1,3) (number of projections, 1, x/y/z)

                euclidian_distance_vertex_to_atom = bfg.Edist3D(atom_coords_projections_transpose, neighborhood_grid_vertices_coords_array)

                #output: array (N,3,3)(triangle number, triangle vertice number, x/y/z)
                #each of the three points of the triangle get subtracted from the atom_coords_projections
                #Defines the vector from the vertex of the triangle to the coordiantes of the projected atom 
                vertex_to_atom = atom_coords_projections[:, None, :] - neighborhood_grid_vertices_coords_array 

                #each column is the area of the parallelogram that is formed by the vector of the atom and each edge of the triangle
                area_parallelogram = np.linalg.norm(np.cross(vertex_to_atom, neighborhood_edges_array, axis=2), axis=2)
                
                #computes the magnitude of the vector for each edge of a given triangle
                base_parallelogram = np.linalg.norm(neighborhood_edges_array, axis=2)  # base of parallelogram
                
                #shortest distance from atom to each edge of triangle
                shortest_distance_to_atom = np.divide(area_parallelogram, base_parallelogram)   
                

                min_distance_any_grid = np.min(shortest_distance_to_atom)

                if min_distance_any_grid > .34796:
                    print(f"Terminated because atom#{atom_num} did not map to any grid!")
                    break
                #print("max shortest distance to atom:", np.max(shortest_distance_to_atom))
    
                #now loop over each triangle in the array if shortest distance < radius assign to triangle then execute the following 
                #loop over grid_ids 
                #each row in the array is a triangle grid in the grids_ids

                for i, grid_id in enumerate(grid_ids_sorted):
                    # Convert grid ID to row index
                    row = grid_to_row[grid_id]
                    if (
                        grid_id not in triangle_atom_assignment_upper
                    ):  # only initialize the grid if the grid doesn't exist yet
                        triangle_atom_assignment_upper[grid_id] = 0

                    #get shortest distance for a given triangle
                    shortest_distance_to_atom_all_edges = shortest_distance_to_atom[row] #gives me the shortest distance to each edge
                    
                    # Get the minimum distance from atom to all edges of the triangle
                    min_distance_to_edges = np.min(shortest_distance_to_atom_all_edges)
                    
                    #print("min distance to edges:", min_distance_to_edges, "atom radius:", atm_radius)

                    #if any of the shortest distances are smaller than atom radius complete the following
                    if min_distance_to_edges <= atm_radius:
                        if (
                            atm_type == "n"
                        ):  
                            triangle_atom_assignment_upper[grid_id] += 1.0
        
                        else:
                            triangle_atom_assignment_upper[grid_id
                            ] += 0.001

    for grid in grid_num:
        grid_int = int(grid)  # normalize key type
        if grid_int not in triangle_atom_assignment_upper:
            print(f"Adding missing grid: {grid_int}")
            triangle_atom_assignment_upper[grid_int] = 0.0

    triangle_atom_assignment_upper_sorted = dict(
        sorted(triangle_atom_assignment_upper.items())
    )
    
    with open(
        util.trajectory_path
        / f"test_DPPC_triangle_atom_assignment_upper_strain2_frame{frame}.pickle",
        "wb",
    ) as f:
        pickle.dump(triangle_atom_assignment_upper_sorted, f)

    # use csv file to check code is working
    with open(
        util.trajectory_path
        / f"test_DPPC_triangle_atom_assignment_upper_strain2_frame{frame}.csv",
        "w",
        newline="",
    ) as f:
        writer = csv.writer(f)
        writer.writerow(["triangle_id", "assignment"])  # optional header
        for key, value in triangle_atom_assignment_upper_sorted.items():
            writer.writerow([key, value])

    return triangle_atom_assignment_upper_sorted

import pickle
import numpy as np
import csv
import util
import neighborhood_grid_graph as ngg 

def sum_atom_mapping_upper_np_array(atom_features:dict, radius:dict,aliphatic:dict, 
                                    lipid_ids_all_grid_edges:dict, lipid_ids_all_grid_coords:dict, lipid_ids_all_norms_upper: dict, lipid_ids_all_centroid_coords_upper:dict, lipid_grid_assignment:dict,
                                    atom_neighborhood_grid_id:dict, grid_num:list, frame:int):
    '''
    For each atom calculate the neighborhood assignment all at once using numpy arrays 

    Args:
        atom_features (dict): A dictionary containing features of atoms, including their coordinates and type.output ex: {atm#: GL2', 'lower', array([18.222,  4.37 , 14.85 ]), 'DOPC GL2', resnumber }
        radius (dict): A dictionary mapping atom names to their respective radii.
        aliphatic (dict): A dictionary mapping atom names to their type (e.g., polar or non-polar).
        
        lipid_ids_all_grid_edges (dict): {res_id: array of all neighborhood triangles edges}
        lipid_ids_all_grid_coords (dict): {res_id: array of all neighborhood vertex of the triangles}
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
       
        if atom_num < 3: #change if want to troubleshoot 
            atom_num = str(
                atom_num
            )  # changes atom_num to a string so it matches atom_features keys

            ##checks to make sure the grid ids match in the array and the grid_ids
            #print("atom num:", atom_num)
            #print("grid ids:", grid_ids)
            #print("neighborhood:",neighborhood_edges_dict.keys())

            # determine various characteristics of the atom 
            lipid_atm_name = atom_features[atom_num][3] #determine the resdiue number of the atom
            #print(lipid_atm_name)
            atom_coords = atom_features[atom_num][2] #gets the coordianates of the atom
            print("atom number:", atom_num)
            print("atom cords:", atom_coords)
            #print("atom coords:", atom_coords)
            atm_radius = [radius[lipid_atm_name]] #gets the radius of the atom 
            print("atom radius:", atm_radius)
            atm_type = aliphatic[lipid_atm_name] #gets the atom type (polar vs non-polar)

            #deals with an atom number that does not exist in the atom features dictionary
            if atom_num not in atom_features:
                print(f"Warning: atom_num {atom_num} missing in atom_features")
                # fails.append(int(atom_num))
                continue

            #Handeling of only atoms that are in the upper leaflet
            if (
                atom_features[atom_num][1] == "lower"
            ):  ##skips atoms that belong to the lower leaflet
                continue

            else:
                
                #pulls the array dependent on the lipid residue number the atom belongs to
                #pulls the array dependent on the lipid residue number the atom belongs to
                neighborhood_edges_array = lipid_ids_all_grid_edges[res_num][order]
                neighborhood_coords_array = lipid_ids_all_grid_coords[res_num][order]
                neighborhood_norms_array = lipid_ids_all_norms_upper[res_num][order]
                neighborhood_centroid_coords_array = lipid_ids_all_centroid_coords_upper[res_num][order]

                #print("neighborhood_centroid_coords_array:", neighborhood_centroid_coords_array)

                ####Project the atom onto the same plane as the triangle
                ####Something is wrong here

                # Calculate vector from triangle centroids to atom
                w_vectors =  atom_coords - neighborhood_centroid_coords_array


                # Calculate the normal direction by computing the dot product of the vector from the triangle centroid to the atom and the normalized triangle normal
                #sum across rows (sum of x,y,z components)
                #scalar projection of w onto the normal vector 
                dot_products = np.sum(w_vectors * neighborhood_norms_array, axis=1) #computes the dot product from each w_vector and the normals of the corresponding grid (output is a scalar)
                
                #neighborhood_norms_array (N, 3)
                #dot_products is a scalar 
                #dot_products[:, None] changes from row of number to column of numbers to do the multiplication
                #(3,1) * (3,3)
                #gives you the normal component of w 
                n_contribution = dot_products[:, None] * neighborhood_norms_array

                #each row in the array represents the projection of the atom onto a different triangle plane
                #get the vector parallel to the plane 
                #describes: How far do I move from the atom to reach its projection onto the plane

                #-----either I can move the atom down or move the centroid over, both should give the same results-----
                #check this 
                #more intuitive to take atom and move down the normal component of w

                #takes the original coordinate of atom and moves it by atom_coords_projections
                atom_coords_projections = atom_coords - n_contribution

                #graph of projected atom coordinates, grid centroid, and the vertices of the triangle 
                ngg.projected_atom_grid(atom_coords_projections, neighborhood_centroid_coords_array,neighborhood_coords_array, atom_coords, neighborhood_norms_array, atom_num) 

                
                
                
                
                #-----Determine if the atom lies within the triangle-----

                #output: array (N,3,3)(triangle number, point number, x/y/z)
                vertex_to_atom = atom_coords_projections[:, None, :] - neighborhood_coords_array #each of the three points of the triangle get subtracted from the atom_coords_projections

            
                #each row is a different triangle
                #each column is the area of the parallelogram that is formed by the vector of the atom and each edge of the triangle
                #issue: neighborhood edges are off?
                area_parallelogram = np.linalg.norm(np.cross(vertex_to_atom, neighborhood_edges_array, axis=2), axis=2) 
                
                #computes the magnitude of the vector for each edge of a given triangle
                base_parallelogram = np.linalg.norm(neighborhood_edges_array, axis=2)  # base of parallelogram
                
                #shortest distance from atom to each edge of triangle
                ##issue: some atoms are not mapping anywhere even though they should be
                shortest_distance_to_atom = np.divide(area_parallelogram, base_parallelogram)   
                

                #something is wrong here because the shortest distance to should at least be smaller than the atom radius for some grid
                min_distance_any_grid = np.min(shortest_distance_to_atom)
                #print("min shortest distance to atom:", min_distance_any_grid)
                if min_distance_any_grid > 3.4796:
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
                    min_dist_per_atom = np.min(shortest_distance_to_atom_all_edges) #gives me the minimum distance for each triangle, taking into account all edges of the triangle


                    #if any of the shortest distances are smaller than atom radius complete the following
                    if min_dist_per_atom <= atm_radius:
                        if (
                            atm_type == "n"
                        ):  # add an if statement where if the atom is polar exit and move onto a new triangle
                            triangle_atom_assignment_upper[grid_id] += 1.0
        
                        else:
                            triangle_atom_assignment_upper[grid_id
                            ] += 0.001

    for grid in grid_num:
        grid_int = int(grid)  # normalize key type
        if grid_int not in triangle_atom_assignment_upper:
            # print(f"Adding missing grid: {grid_int}")
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

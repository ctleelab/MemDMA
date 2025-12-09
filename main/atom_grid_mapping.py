import pickle
import numpy as np
import csv
import util
import neighborhood_grid_graph as ngg
import BasicFunctions as bfg

def sum_atom_mapping_upper_np_array(atom_features:dict, radius:dict,aliphatic:dict, 
                                    lipid_ids_all_grid_edges:dict, lipid_ids_all_grid_vertices_coords:dict, lipid_ids_all_norms_upper: dict, lipid_ids_all_centroid_coords_upper:dict, lipid_grid_assignment:dict, lipid_ids_all_grid_edges_2coords:dict, 
                                    atom_neighborhood_grid_id:dict, grid_num:list, frame:int):
    '''
    For each atom calculate the mapping value for all grids in a neighborhood at once using numpy arrays 

    Args:
        atom_features (dict): A dictionary containing features of atoms, including their coordinates and type.output ex: {atm#: GL2', 'lower', array([18.222,  4.37 , 14.85 ]), 'DOPC GL2', resnumber }
        radius (dict): A dictionary mapping atom names to their respective radii.
        aliphatic (dict): A dictionary mapping atom names to their type (e.g., polar or non-polar).
        
        lipid_ids_all_grid_edges (dict): {res_id: array of all neighborhood triangles edges represented as a single vector}
        lipid_ids_all_grid_vertices_coords (dict): {res_id: array of all neighborhood vertex of the triangles}
        lipid_ids_all_norms_upper (dict): {res_id: array of all neighborhood triangles normal vectors}
        lipid_ids_all_centroid_coords_upper(dict): A dictionary {res_id: array of all neighborhood triangles centroid coordinates}
        lipid_grid_assignment (dict): A dictionary with lipid residue number (key): (neighborhood_list of the grid_idex)
        ipid_ids_all_grid_edges_2coords (dict): {res_id: array of all neighborhood triangles edges represented as two vertex coordinates}

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
       
        if atom_num < 3: #change if want to troubleshoot so only do a set number of atoms
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

                #Determine the distance from atom projection to each vertex of the triangle
                #If the distance shortest than radius of the atom, it is considered a hit 
                #Input: atom_coords_projections (N,3), neighborhood_grid_vertices_coords_array (N,3,3) (number of triangles, 3 vertices, x/y/z)


                #graph of determining if the atom lies within the triangle
                #ngg.atom_projection_grid_assignment(atom_coords_projections, neighborhood_grid_vertices_coords_array, vertex_to_atom) 


                # -----Dealing with special cases where you cannot use the triangle vertex to provide atom assignment-----

                #Determine what type of special case you're working with 
                #1)dot product between edge of triangle and vector from vertex to atom projection is negative
                #2)dot product between edge of triangle and vector from vertex to atom is positive, and less than the magnitude of edge vector
                #3)dot product between edge of triangle and vector from vertex to atom is positive, and greater than the magnitue of edge vector 

                #Computes the vector from each vertex of the traingle to the atom projection
                vector_from_vetex_to_atom_projection = atom_coords_projections[:, None, :] - neighborhood_grid_vertices_coords_array


                #Calculate the euclidian distance from each vertex of the triangle to the atom projection
                atom_coords_projections_transpose = atom_coords_projections[:, None, :] # reshape to (N,1,3) (number of projections, 1, x/y/z)
                euclidian_distance_vertex_to_atom = bfg.Edist3D(atom_coords_projections_transpose, neighborhood_grid_vertices_coords_array)

                #dot product between the edge vector and the vector from vertex to atom projection
                #shape (N,3) (number of triangles, 3 vertices)
                #output: the dot product between each edge of the triangle and the vector from the vertex to the atom projection
                dot_product_edge_vertex_to_atom_projection = np.sum(neighborhood_edges_array * vector_from_vetex_to_atom_projection, axis=2) #shape (N,3) number of triangles, 3 edges


                #compute the magnitude of each edge vector of the triangle
                #this can most likely be precomputed to save time
                edge_magnitude = np.linalg.norm(neighborhood_edges_array, axis=2)

                for i, grid_id in enumerate(grid_ids_sorted):
                    # Convert grid ID to row index
                    row = grid_to_row[grid_id]
                    if (
                        grid_id not in triangle_atom_assignment_upper
                    ):  # only initialize the grid if the grid doesn't exist yet
                        triangle_atom_assignment_upper[grid_id] = 0

                    #-----deal with the vertices cases first-----
                    atom_maps_to_triangle = False 

                    #gets all distances from atom to all vertices of a given triangle
                    shortest_distance_to_atom_all_vertex = euclidian_distance_vertex_to_atom[row] 
                    
                    # Get the minimum distance from atom to all vertices of the triangle
                    min_distance_to_vertices = np.min(shortest_distance_to_atom_all_vertex)   

                    if min_distance_to_vertices <= atm_radius:
                        print("atom maps to triangle")
                        atom_maps_to_triangle = True   
                    
                    if min_distance_to_vertices > atm_radius:
                        print("atom does not map to triangle by looking at vertices, checking special cases now")

                        # #need to evaluate the other special cases below

                        # tri_vertex = dot_product_edge_vertex_to_atom_projection[grid_id]

                        # for tri_vertex_indx, dotproduct in enumerate(tri_vertex):
                        #     if dotproduct > 0:

                        #         #compare to the magnitude of the edge vector
                        #         edge_magnitude_tri_vertex_indx = edge_magnitude[grid_id][tri_vertex_indx]

                        #         if dotproduct <= edge_magnitude_tri_vertex_indx:
                        #             #compute the area of the parallelogram and divide by the base to get the hieght 

                        #             #compute the area of the parallelogram formed by the edge vector and the atom projection vector 
                        #             area_parallelogram = np.linalg.norm(np.cross(vector_from_vetex_to_atom_projection[grid_id][tri_vertex_indx], neighborhood_edges_array[grid_id][tri_vertex_indx]), axis=0)
                                    
                        #             #computes the magintude of the edge vector, which is the base of the paralleogram
                        #             base_parallelogram = np.linalg.norm(neighborhood_edges_array[grid_id][tri_vertex_indx]) 
                                    
                        #             height_parallelogram = area_parallelogram / base_parallelogram

                        #             if height_parallelogram <= atm_radius:
                        #                 atom_maps_to_triangle = True

                        #         #deals with atoms that are beyond the edge in the positive direction 
                        #         if dotproduct > edge_magnitude_tri_vertex_indx:

                        #             minimum_distance_to_vertex_triangle = np.min(euclidian_distance_vertex_to_atom[triangle_vertices])

                        #             #hadling the assignment value based on atom type
                        #             if minimum_distance_to_vertex_triangle <= atm_radius:
                        #                 atom_maps_to_triangle = True


                        #     #deals with atoms that are beyond the edge in the negative direction
                        #     elif dotproduct < 0:
                        #         #print(f"negative dot product case", dotproduct)
                        #         minimum_distance_to_vertex_triangle = np.min(euclidian_distance_vertex_to_atom[triangle_vertices])
                        #         if minimum_distance_to_vertex_triangle <= atm_radius:

                        #             #hadling the assignment value based on atom type
                        #             if minimum_distance_to_vertex_triangle <= atm_radius:
                        #                 atom_maps_to_triangle = True
                    
                    if atom_maps_to_triangle:
                        triangle_atom_assignment_upper[grid_id] += (1.0 if atm_type=='n' else 0.001)

                # if atom_maps_to_triangle == False:
                #     print(f"Atom#{atom_num} does not map to any triangle in the neighborhood!")
                #     break                      

        

                #Construct the appropriate vectors to do the dot product calculation
                #Defining a vector from the vertex of the triangle to the atom projection


                # #output: array (N,3,3)(triangle number, triangle vertice number, x/y/z)
                # #each of the three points of the triangle get subtracted from the atom_coords_projections
                # #Defines the vector from the vertex of the triangle to the coordiantes of the projected atom 
                # vertex_to_atom = atom_coords_projections[:, None, :] - neighborhood_grid_vertices_coords_array 

                # #each column is the area of the parallelogram that is formed by the vector of the atom and each edge of the triangle
                # area_parallelogram = np.linalg.norm(np.cross(vertex_to_atom, neighborhood_edges_array, axis=2), axis=2)
                
                # #computes the magnitude of the vector for each edge of a given triangle
                # base_parallelogram = np.linalg.norm(neighborhood_edges_array, axis=2)  # base of parallelogram
                
                # #shortest distance from atom to each edge of triangle
                # shortest_distance_to_atom = np.divide(area_parallelogram, base_parallelogram)   
                

                # min_distance_any_grid = np.min(shortest_distance_to_atom)

                #testing the euclidian distance calculation



                # for i, grid_id in enumerate(grid_ids_sorted):
                #     # Convert grid ID to row index
                #     row = grid_to_row[grid_id]
                #     if (
                #         grid_id not in triangle_atom_assignment_upper
                #     ):  # only initialize the grid if the grid doesn't exist yet
                #         triangle_atom_assignment_upper[grid_id] = 0

                #     #get shortest distance for a given triangle
                #     shortest_distance_to_atom_all_vertex = euclidian_distance_vertex_to_atom[row] #gives me the shortest distance from the atom to a vertex of
                    
                #     # Get the minimum distance from atom to all edges of the triangle
                #     min_distance_to_edges = np.min(shortest_distance_to_atom_all_vertex)
                    
                #     #print("min distance to edges:", min_distance_to_edges, "atom radius:", atm_radius)

                #     #if any of the shortest distances are smaller than atom radius complete the following
                #     if min_distance_to_edges <= atm_radius:
                #         if (
                #             atm_type == "n"
                #         ):  
                #             triangle_atom_assignment_upper[grid_id] += 1.0
        
                #         else:
                #             triangle_atom_assignment_upper[grid_id
                #             ] += 0.001

    for grid in grid_num:
        grid_int = int(grid)  # normalize key type
        if grid_int not in triangle_atom_assignment_upper:
            #print(f"Adding missing grid: {grid_int}")
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

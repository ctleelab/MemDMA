import pickle 
import numpy as np
import csv
import util


########################################################
##now determine what atoms lie inside in each triangle##
########################################################



# with open("neighborhood_projections_upper_test.pickle", "rb") as f:
#     neighborhood_projections_upper = pickle.load(f)

# with open("neighborhood_projections_atom_id_test.pickle", "rb") as f:
#     neighborhood_grid_id_upper = pickle.load(f)


triangle_atom_assignment_upper  = {} 
triangle_atom_assignment_lower  = {} 

####################
##trouble-shooting##  
####################
#success = [] 
#fails = [] 
#print(max(neighborhood_grid_id_upper))


#function is for the upper leaflet
#inputs: (atom_features, radius, aliphatic, edges_triangles_upper)
#output {grid #: mapping value}
#neighborhood_grid_id: {atom number: list of grid number}
def sum_atom_mapping_upper(atom_features: dict, radius: dict, aliphatic:dict, edges_triangles:dict, atom_neighborhood_grid_id:dict, coords_triangles:dict, grid_num:list, grid_coords_np, norm_vector_triangles:dict, frame: int) -> dict:
    """
    This function is mapping atoms to grids in the upper leaflet of a lipid bilayer.

    Args:
        atom_features (dict): A dictionary containing features of atoms, including their coordinates and type.output ex: {atm#: GL2', 'lower', array([18.222,  4.37 , 14.85 ]), 'DOPC GL2', resnumber }
        radius (dict): A dictionary mapping atom names to their respective radii.
        aliphatic (dict): A dictionary mapping atom names to their type (e.g., polar or non-polar).
        edges_triangles (dict): A dictionary containing edge vectors for triangles as a list. {triangle_number: [edge1, edge2, edge3]}
        atom_neighborhood_grid_id (dict): A dictionary containing the res_number the atom belongs to, coords of atom, list of close by grids for each atom number
        coords_triangles (dict): A dictionary containing the coordinates of all vertexes of the triangles {14959: array([[113.63737488,  21.55670166, 282.9483378 ],[114.63858032,  21.55670166, 282.60440123],[114.63858032,  22.55933762, 282.61294111]])}
        grid_num (list): a list of all the grids numbers in either the upper or lower leaflet
        grid_coords_np (array): An array where is index represents the grid number, and has the (x,y,z) coordinates of the grid in the designated leaflet
        norm_vector_triangles (dict): A dictionary with triangle number a nd the coordinates of the normal vector (28104: array([ 0.29990481, -0.31881054,  0.89912009])
        frame (int): the frame number of the trajectory file     

    Returns:
        dict: A dictionary mapping triangle IDs to their respective atom assignment values.
    """
    for atom_num, data in atom_neighborhood_grid_id.items(): 
        grid_ids = data[2] 
        #print("atomnum:", atom_num)
        atom_num = str(atom_num) #changes atom_num to a string so it matches atom_features keys
        if atom_num not in atom_features:
            print(f"Warning: atom_num {atom_num} missing in atom_features")
                #fails.append(int(atom_num))
            continue
        #print(atom_features[atom_num][1])
            # Determine the coordinates for the atom by taking atom_id and searching in the same triangle
            # in neighborhood_projections and pulling out its coordinates
        if atom_features[atom_num][1] == "lower": ##skips atoms that belong to the lower leaflet
            break
        else:
                #success.append(int(atom_num))
            for grid_id in grid_ids: #loops over the neighborhood of grids for each atom
                shortest_distance_to_atom = 1000 #setting arbritraury distance that much greater than the radius of all atoms 


                atom_coords = atom_features[atom_num][2]*10

                ###Determine the coords if atom projected onto the same plane as the triangle grid
                ##inputs: atom_num, grid_id, 

                #subtract 1 because need to index the numpy grid that starts at 0 and grid_id starts at 1
                grid_coord = grid_coords_np[grid_id-1] 
                w_vectors = grid_coord - atom_coords  #change the naming of the grid_coords to triangle_centroid_coords 

                # scalars = np.dot(w_vectors, normals[tri_indx+1])
                # n_contributions = np.outer(scalars, normals[tri_indx+1])

                # Calculate the normal direction by computing the dot product of the vector from the triangle centroid to the atom and the normalized triangle normal
                n_contributions = np.dot(w_vectors, norm_vector_triangles[grid_id]) * norm_vector_triangles[grid_id]

                # Project the atom onto the triangle plane by removing the normal contribution
                atom_coords_projections = w_vectors - n_contributions

                #determine the radius of the atom you're working with
                lipid_atm_name = atom_features[atom_num][3]
                atm_radius = radius[lipid_atm_name]

                #print("atom_radius:", atm_radius)
                atm_type = aliphatic[lipid_atm_name]
                #print(atm_type)

                triangle_idex_num = int(grid_id)

                ##import edges of triangle and import points of triangles
                edges_and_points = {
                    "edge1": (edges_triangles[grid_id][0], coords_triangles[grid_id][0]),
                    "edge2": (edges_triangles[grid_id][1], coords_triangles[grid_id][1]),
                    "edge3": (edges_triangles[grid_id][2], coords_triangles[grid_id][2])
                }

                ##determine the distance between each edge
                ##once the distance is <1 exit the loop of the code 
                
                if triangle_idex_num not in triangle_atom_assignment_upper: #only initialize the grid if the grid doesn't exist yet
                    triangle_atom_assignment_upper[triangle_idex_num] = 0

                ##the plane should be the one defined by the triangle 
                ##point_on_edge: is a vertex of the triangle
                for edge, (edge_vector, point_on_edge) in edges_and_points.items():

                    #triangle_atom_assignment_upper[triangle_idex_num] = []
                    if shortest_distance_to_atom > atm_radius:
                        #vector from point of vector to coords_atom_plane_of_triangle
                        
                        vector_to_atom =  atom_coords_projections-point_on_edge 
                        area_parallelogram = np.linalg.norm(np.cross(vector_to_atom, edge_vector)) 
                        base_parallelogram = np.linalg.norm(edge_vector) #length of parallelogram
                        shortest_distance_to_atom = np.divide(area_parallelogram,base_parallelogram)

                        #print("shortest_distance:", shortest_distance_to_atom)
                        #print("shortest:", shortest_distance_to_atom)

                        ###description of why 1.0 vs .001 (these will be the default values)
                        if shortest_distance_to_atom <= atm_radius:
                            if atm_type == "n":  #add an if statement where if the atom is polar exit and move onto a new triangle
                                triangle_atom_assignment_upper[triangle_idex_num] += 1.0 
                                # if triangle_idex_num not in triangle_atom_assignment_upper:
                                #     triangle_atom_assignment_upper[triangle_idex_num] = []
                                #triangle_atom_assignment_upper[triangle_idex_num].append((atom_id).strip()) 
                                #triangle_atom_assignment_upper[triangle_idex_num].append((atom_id).strip()) 
                                ##rather than appends the atom id, maybe just add 1 
                                #stop_triangle = True
                            else:
                                triangle_atom_assignment_upper[triangle_idex_num] += .001
                                # if triangle_idex_num not in triangle_atom_assignment_upper:
                                #     triangle_atom_assignment_upper[triangle_idex_num] = []
                                #triangle_atom_assignment_upper[triangle_idex_num].append((atom_id).strip()) 
                            break # once you hit then you do not need to keep checking other edges


    #check to see grid already present 
    for grid in grid_num:
        grid_int = int(grid)  # normalize key type
        if grid_int not in triangle_atom_assignment_upper:
            #print(f"Adding missing grid: {grid_int}")
            triangle_atom_assignment_upper[grid_int] = 0.0

    triangle_atom_assignment_upper_sorted = dict(sorted(triangle_atom_assignment_upper.items()))

    
    ####################
    ##trouble-shooting##  
    ####################
    #print(fails)
    #print("min fails:", min(fails))
    #print("max fails:", max(fails))
    #print(success)
    #print("min success:", min(success))
    #print("max success:", max(success))
    #print("length of atom features:", len(atom_features))
                    
    # Save as a pickle file
    with open(util.trajectory_path/f"DPPC_triangle_atom_assignment_upper_strain2_frame{frame}.pickle", "wb") as f:
        pickle.dump(triangle_atom_assignment_upper_sorted, f) 

    #use csv file to check code is working
    with open(util.trajectory_path/f"DPPC_triangle_atom_assignment_upper_strain2_frame{frame}.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["triangle_id", "assignment"])  # optional header
        for key, value in triangle_atom_assignment_upper_sorted.items():
            writer.writerow([key, value])
    
    return triangle_atom_assignment_upper_sorted


def sum_atom_mapping_lower(atom_features, radius, aliphatic, edges_triangles_lower, neighborhood_grid_id_lower, coords_triangles_lower, grid_num): 
    for atom_num, grid_ids in neighborhood_grid_id_lower.items():  
        #print("atomnum:", atom_num)
        stop_triangle = False #just emplaced so that I can troubleshoot the code
        atom_num = str(atom_num) #changes atom_num to a string so it matches atom_features keys
        #print("atom_num:", atom_num)
        if atom_num not in atom_features:
            print(f"Warning: atom_num {atom_num} missing in atom_features")
                #fails.append(int(atom_num))
            continue
        #print(atom_features[atom_num][1])
            # Determine the coordinates for the atom by taking atom_id and searching in the same triangle
            # in neighborhood_projections and pulling out its coordinates
        if atom_features[atom_num][1] == "upper": ##skips atoms that belong to the upper leaflet
            break
        else:
                #success.append(int(atom_num))
            for grid_id in grid_ids: #loops over the neighborhood of grids for each atom
                shortest_distance_to_atom = 1000 #setting arbitrary distance that much greater than the radius of all atoms 
                grid_id += 1 
                if stop_triangle:  #will move onto a new triangle if an polar atom was mapped
                    break


                #coordinates of atom in plane of triangle 
                # do the projection here, rather the indexing 
                atom_coords = matching.get_projection(atom_num, grid_id, neighborhood_projections_upper, neighborhood_grid_id_upper, atom_features)

                #atom_coords = atom_features[atom_num][2]*10

                #determine the radius of the atom you're working with
                lipid_atm_name = atom_features[atom_num][3]
                atm_radius = radius[lipid_atm_name]
                print("atom_radius:", atm_radius)
                atm_type = aliphatic[lipid_atm_name]
                #print(atm_type)

                triangle_idex_num = int(grid_id)

                ##import edges of triangle and import points of tringles
                edges_and_points = {
                    "edge1": (edges_triangles_lower[grid_id][0], coords_triangles_lower[grid_id][0]),
                    "edge2": (edges_triangles_lower[grid_id][1], coords_triangles_lower[grid_id][1]),
                    "edge3": (edges_triangles_lower[grid_id][2], coords_triangles_lower[grid_id][2])
                }

                ##determine the distance between each edge
                ##once the distance is <1 exit the loop of the code 
                
                if triangle_idex_num not in triangle_atom_assignment_lower: #only initialize the grid if the grid doesn't exist yet
                    triangle_atom_assignment_upper[triangle_idex_num] = 0

                ##the plane should be the one defined by the triangle 
                for edge, (edge_vector, point_on_edge) in edges_and_points.items():

                    #triangle_atom_assignment_upper[triangle_idex_num] = []
                    if shortest_distance_to_atom > atm_radius:
                        #vector from point of vector to coords_atom_plane_of_triangle
                        
                        vector_to_atom = atom_coords-point_on_edge
                        area_parallelogram = np.linalg.norm(np.cross(vector_to_atom, edge_vector)) 
                        base_parallelogram = np.linalg.norm(edge_vector)
                        shortest_distance_to_atom = np.divide(area_parallelogram,base_parallelogram)

                        print("shortest_distance:", shortest_distance_to_atom)
                        #print("shortest:", shortest_distance_to_atom)

                        if shortest_distance_to_atom <= atm_radius:
                            if atm_type == "n":  #add an if statement where if the atom is polar exit and move onto a new triangle
                                triangle_atom_assignment_lower[triangle_idex_num] += 1.0 
                                # if triangle_idex_num not in triangle_atom_assignment_upper:
                                #     triangle_atom_assignment_upper[triangle_idex_num] = []
                                #triangle_atom_assignment_upper[triangle_idex_num].append((atom_id).strip()) 
                                #triangle_atom_assignment_upper[triangle_idex_num].append((atom_id).strip()) 
                                ##rather than appends the atom id, maybe just add 1 
                                #stop_triangle = True
                                continue 
                            else:
                                triangle_atom_assignment_lower[triangle_idex_num] += .001
                                # if triangle_idex_num not in triangle_atom_assignment_upper:
                                #     triangle_atom_assignment_upper[triangle_idex_num] = []
                                #triangle_atom_assignment_upper[triangle_idex_num].append((atom_id).strip()) 
                                continue


    #check to see grid already present 
    for grid in grid_num:
        grid_int = int(grid)  # normalize key type
        if grid_int not in triangle_atom_assignment_lower:
            #print(f"Adding missing grid: {grid_int}")
            triangle_atom_assignment_lower[grid_int] = 0.0

    triangle_atom_assignment_lower_sorted = dict(sorted(triangle_atom_assignment_lower.items()))

    
    ####################
    ##trouble-shooting##  
    ####################
    #print(fails)
    #print("min fails:", min(fails))
    #print("max fails:", max(fails))
    #print(success)
    #print("min success:", min(success))
    #print("max success:", max(success))
    #print("length of atom features:", len(atom_features))
                    
    # Save as a pickle file
    with open("triangle_atom_assignment_lower_flat.pickle", "wb") as f:
        pickle.dump(triangle_atom_assignment_lower_sorted, f) 

    with open("triangle_atom_assignment_lower_flat.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["triangle_id", "assignment"])  # optional header
        for key, value in triangle_atom_assignment_lower_sorted.items():
            writer.writerow([key, value])
    
    return triangle_atom_assignment_lower_sorted
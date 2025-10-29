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
def sum_atom_mapping_upper(atom_features: dict, radius: dict, aliphatic:dict, edges_triangles_upper, neighborhood_grid_id_upper, coords_triangles_upper, grid_num, frame: int) -> dict:
    """
    This function is mapping atoms to grids in the upper leaflet of a lipid bilayer.

    Args:
        atom_features (dict): A dictionary containing features of atoms, including their coordinates and type.
        radius (dict): A dictionary mapping atom names to their respective radii.
        aliphatic (dict): A dictionary mapping atom names to their type (e.g., polar or non-polar).
        edges_triangles_upper: Data structure containing edge vectors for triangles in the upper leaflet.
        neighborhood_grid_id_upper: A dictionary that assigns a list of grid per atom  [neighborhood_grid_id: {atom number: list of grid number}]
        coords_triangles_upper: Data structure containing coordinates of triangle vertices in the upper leaflet.        

    Returns:
        dict: A dictionary mapping triangle IDs to their respective atom assignment values.
    """
    for atom_num, grid_ids in neighborhood_grid_id_upper.items():  
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
        if atom_features[atom_num][1] == "lower": ##skips atoms that belong to the lower leaflet
            break
        else:
                #success.append(int(atom_num))
            for grid_id in grid_ids: #loops over the neighborhood of grids for each atom
                shortest_distance_to_atom = 1000 #setting arbritraury distance that much greater than the radius of all atoms 
                grid_id += 1 
                if stop_triangle:  #will move onto a new triangle if an polar atom was mapped
                    break


                #coordinates of atom in plane of triangle 
                #atom_coords = matching.get_projection(atom_num, grid_id, neighborhood_projections_upper, neighborhood_grid_id_upper, atom_features)
                
                #multiply by 10 to covert into Angstrom
                #
                atom_coords = atom_features[atom_num][2]*10


                #determine the radius of the atom you're working with
                lipid_atm_name = atom_features[atom_num][3]
                atm_radius = radius[lipid_atm_name]

                #print("atom_radius:", atm_radius)
                atm_type = aliphatic[lipid_atm_name]
                #print(atm_type)

                triangle_idex_num = int(grid_id)

                ##import edges of triangle and import points of triangles
                edges_and_points = {
                    "edge1": (edges_triangles_upper[grid_id][0], coords_triangles_upper[grid_id][0]),
                    "edge2": (edges_triangles_upper[grid_id][1], coords_triangles_upper[grid_id][1]),
                    "edge3": (edges_triangles_upper[grid_id][2], coords_triangles_upper[grid_id][2])
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
                        
                        vector_to_atom = atom_coords-point_on_edge 
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
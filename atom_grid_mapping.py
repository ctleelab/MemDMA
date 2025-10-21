import pickle 
import numpy as np
import csv


########################################################
##now determine what atoms lie inside in each triangle##
########################################################



# with open("neighborhood_projections_upper_test.pickle", "rb") as f:
#     neighborhood_projections_upper = pickle.load(f)

# with open("neighborhood_projections_atom_id_test.pickle", "rb") as f:
#     neighborhood_grid_id_upper = pickle.load(f)


triangle_atom_assignment_upper  = {} 
triangle_atom_assignment_lower  = {} 
#success = [] #troubleshooting
#fails = [] #troublshooting




#trouble-shooting    
#print(max(neighborhood_grid_id_upper))



#function
#inputs: (atom_features, radius, aliphatic, edges_triangles_upper)


#output grid, mapping value

def sum_atom_mapping(atom_features, radius, aliphatic, edges_triangles_upper, neighborhood_grid_id_upper, coords_triangles_upper): 
    for atom_num, grid_ids in neighborhood_grid_id_upper.items(): 
        stop_triangle = False
        atom_num = int(atom_num)
        if atom_num > 0: #just emplaced so that I can troubleshoot the code
            atom_num = str(atom_num) #changes atom_num to a string so it matches atom_features keys
            if atom_num not in atom_features:
                #print(f"Warning: atom_num {atom_num} missing in atom_features")
                #fails.append(int(atom_num))
                continue
            # Determine the coordinates for the atom by taking atom_id and searching in the same triangle
            # in neighborhood_projections and pulling out its coordinates
            else:
                #success.append(int(atom_num))
                for grid_id in grid_ids: #loops over the neighborhood of grids for each atom
                    ##Issue: for some reason atoms are mising
                    if stop_triangle:   #will move onto a new triangle if an polar atom was mapped
                        break
                    if atom_features[atom_num][1] == "lower": ##skips atoms that belong to the lower leaflet
                        continue

                    #coordinates of atom in plane of triangle 
                    #atom_coords = matching.get_projection(atom_num, grid_id, neighborhood_projections_upper, neighborhood_grid_id_upper, atom_features)

                    atom_coords = atom_features[atom_num][2]

                    #determine the radius of the atom you're working with
                    lipid_atm_name = atom_features[atom_num][3]
                    atm_radius = radius[lipid_atm_name]
                    atm_type = aliphatic[lipid_atm_name]


                    triangle_idex_num = int(grid_id)

                    ##import edges of triangle and import points of tringles
                    edges_and_points = {
                        "edge1": (edges_triangles_upper[grid_id][0], coords_triangles_upper[grid_id][0]),
                        "edge2": (edges_triangles_upper[grid_id][1], coords_triangles_upper[grid_id][1]),
                        "edge3": (edges_triangles_upper[grid_id][2], coords_triangles_upper[grid_id][2])
                    }

                ##determine the distance between each edge
                ##once the distance is <1 exit the loop of the code 
                    shortest_distance_to_atom = 1000 #setting arbritraury distance that much greater than the radius of all atoms 

                    if triangle_idex_num not in triangle_atom_assignment_upper: #only initialize the grid if the grid doesn't exist yet
                        triangle_atom_assignment_upper[triangle_idex_num] = 0

                    ##the plane should be the one defined by the triangle...
                    for edge, (edge_vector, point_on_edge) in edges_and_points.items():
                        #triangle_atom_assignment_upper[triangle_idex_num] = []
                        if shortest_distance_to_atom > atm_radius:
                            #vector from point of vector to coords_atom_plane_of_triangle
                            
                            vector_to_atom = atom_coords- point_on_edge
                            area_parallelogram = np.linalg.norm(np.cross(vector_to_atom, edge_vector))  #this is giving me coordinates and not a single value....
                            base_parallelogram = np.linalg.norm(edge_vector)

                            shortest_distance_to_atom = np.divide(area_parallelogram,base_parallelogram)
                        if shortest_distance_to_atom <= atm_radius:
                            if atm_type == "n":  #add an if statement where if the atom is polar exit and move onto a new triangle
                                print("polar atom")
                                triangle_atom_assignment_upper[triangle_idex_num] += 1.0 
                                # if triangle_idex_num not in triangle_atom_assignment_upper:
                                #     triangle_atom_assignment_upper[triangle_idex_num] = []
                                #triangle_atom_assignment_upper[triangle_idex_num].append((atom_id).strip()) 
                                #triangle_atom_assignment_upper[triangle_idex_num].append((atom_id).strip()) 
                                ##rather than appends the atom id, maybe just add 1 
                                stop_triangle = True
                                break 
                            else:
                                triangle_atom_assignment_upper[triangle_idex_num] += .001
                                # if triangle_idex_num not in triangle_atom_assignment_upper:
                                #     triangle_atom_assignment_upper[triangle_idex_num] = []
                                #triangle_atom_assignment_upper[triangle_idex_num].append((atom_id).strip()) 
                                continue
    

    #troubleshooting
    #print(fails)
    #print("min fails:", min(fails))
    #print("max fails:", max(fails))
    #print(success)
    #print("min success:", min(success))
    #print("max success:", max(success))
    #print("length of atom features:", len(atom_features))
                    
    # Save as a pickle file
    with open("triangle_atom_assignment_upper_test7.pickle", "wb") as f:
        pickle.dump(triangle_atom_assignment_upper, f) 

    with open("triangle_atom_assignment_upper_test7.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["triangle_id", "assignment"])  # optional header
        for key, value in triangle_atom_assignment_upper.items():
            writer.writerow([key, value])
    
    return 




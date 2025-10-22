#build a list of Z-coords based on the normal vectors

#input: normal vectors to the centroids of the triangular mesh, grid_mapping(tells the grid closest to each lipid)
#output: normal vector to each atom

import numpy as np
import util as util
import re

def atm_normal_vector_assignment(pdblines, residue_leaflet, closest_mapping, norm_up, norm_l):
    normal_vector_atoms = {}

    #need to import the gro file 
    for i, line_atom in enumerate(pdblines):
        if i > 1 and i < (len(pdblines)-1):
            #z-coord is dependent on which lipid the atom belong to 
            #search function to find which atom the lipid belongs to 
            res_name_number = line_atom[0:12].strip()
            match = re.match(r"(\d+)([A-Za-z]+)", res_name_number)
            if match:
                res_number = int(match.group(1))
                #res_name = int(match.group(2))
            
            #list of atom types 
            list_types_atoms = util.list_types_atoms

            atm_name_number = str(line_atom[12:22].strip())

            #extracts atom name and number
            for j in list_types_atoms:
                if atm_name_number.startswith(j):
                    atm_name = j
                    atom_number = atm_name_number[len(j):]
                    break

            #first determine if atom is in the upper or lower leaflet depending on its residue number
            leaflet = find_leaflet(res_number, residue_leaflet)
                
            if leaflet == "upper": ##check the value of the dictionary 
                for k in closest_mapping:
                    if k == res_number: 
                        normal_vector_atom = norm_up[k]
                        normal_vector_atoms[atom_number]= normal_vector_atom
            
            if leaflet == "lower":
                for k in closest_mapping:
                    if k == res_number: 
                        normal_vector_atom = norm_l[k]
                        normal_vector_atoms[atom_number] = normal_vector_atom
    return normal_vector_atoms

#determines if a lipid is in the lower or upper leaflet depending on its residue number
def find_leaflet(res_number, residue_leaflet):
    for leaflet, residues in residue_leaflet.items():
        if res_number in residues:
            return leaflet
        

def atm_centroid_assignment(pdblines, residue_leaflet, closest_mapping, centroid_up, centroid_l):
    centroid_atoms = {}  ##creates a dictionary that tells which centroid will go with what atm

    #need to import the gro file 
    for i, line_atom in enumerate(pdblines):
        if i > 1 and i < (len(pdblines)-1):
            #z-coord is dependent on which lipid the atom belong to 
            #search function to find which atom the lipid belongs to 
            res_name_number = line_atom[0:12].strip()
            match = re.match(r"(\d+)([A-Za-z]+)", res_name_number)
            if match:
                res_number = int(match.group(1))
                #res_name = int(match.group(2))
            
            #list of atom types 
            list_types_atoms = util.list_types_atoms

            atm_name_number = str(line_atom[12:22].strip())

            #extracts atom name and number
            for j in list_types_atoms:
                if atm_name_number.startswith(j):
                    atm_name = j
                    atom_number = atm_name_number[len(j):]
                    break

            #first determine if atom is in the upper or lower leaflet depending on its residue number
            leaflet = find_leaflet(res_number, residue_leaflet)
                
            if leaflet == "upper": ##check the value of the dictionary 
                for k in closest_mapping:
                    if k == res_number: 
                        centroid_atom = centroid_up[k]
                        centroid_atoms[atom_number]= centroid_atom
            
            if leaflet == "lower":
                for k in closest_mapping:
                    if k == res_number: 
                        centroid_atom = centroid_l[k]
                        centroid_atoms[atom_number] = centroid_atom
    return centroid_atoms 



#creating the slices of the normal vector 
#how can I move 1A up the normal vector? 
#do I need to feed in the centroid? 
#input: feed in atom coordinates, need to feed in the grid assignment for the atom, length of how much you incriment along the normal
#output: lower_slices and upper_slices for each atom
def grid_slices(pdblines, residue_leaflet, centroid_atoms, normal_vector_atoms):

    upper_slices = {}
    lower_slices = {} ##{atom: centroids of grid coords (x,y,z)}
    #for each atom, center its normal vector onto its centroid
    #then slice up and down the centroid (for 1A it would be .1)
    #the slice would be 1A x 1A

    #need to import the gro file 
    for i, line_atom in enumerate(pdblines):
        if i > 1 and i < (len(pdblines)-1):
            #z-coord is dependent on which lipid the atom belong to 
            #search function to find which atom the lipid belongs to 
            res_name_number = line_atom[0:12].strip()
            match = re.match(r"(\d+)([A-Za-z]+)", res_name_number)
            if match:
                res_number = int(match.group(1))
                #res_name = int(match.group(2))
            
            #list of atom types 
            list_types_atoms = util.list_types_atoms

            atm_name_number = str(line_atom[12:22].strip())

            #extracts atom name and number
            for j in list_types_atoms:
                if atm_name_number.startswith(j):
                    atm_name = j
                    atom_number = atm_name_number[len(j):]
                    break

            #first determine if atom is in the upper or lower leaflet depending on its residue number
            leaflet = find_leaflet(res_number, residue_leaflet)
            if atm_name != "GL2":
                if leaflet == "upper":
                    #for each atom, center its normal vector onto its centroid
                    centroid = centroid_atoms[atom_number] ##grab atom's centroid
                    normal_vector = normal_vector_atoms[atom_number] ##grab atom's normal vector

                    #turn the normal_vector into a unit vector, so that the increment along the normal is 1A
                    unit_normal = normal_vector/ np.linalg.norm(normal_vector) #np.linalg.norm(normal) computes the magnitude of the vector 

                    #upper limit: slightly above the glycerol bead
                    #lower limit: unsure what this limit should be...maybe just the thickness of the leaflet? or maybe the length of a lipid? 
                    #each point will be the centroid of the triangle 
                    lower_limit = [0,0,0]             #bottom_portion_of_leaflet, this would have to be change to be a correct representation
                    upper_limit = centroid + [1,1,1] #[1,1,1] just represents some distance above the gl2, I will have to change this later
                    #start from gl2 coordinate of that residue 
                    #the centroid will be the first grid
                    lower_slice_coord_centroid = centroid + util.SIZE * unit_normal
                    
                    #I will need to add some offset so that it covers above the centroid (maybe to the PO4 bead?)
                    #centroid is the starting value, and will go to lower_limit decrease by 1A each time along the normal
                    centroids_along_normal_vector = create_list_descend(centroid, lower_limit, unit_normal *-1) #multiply by -1 because you want to descend

                    #need to append the centroids_along_normal_vector for each atom
                    upper_slices[atom_number] = centroids_along_normal_vector

                if leaflet == "lower":
                    continue

            else:
                continue
    print("completed grid_slicing")
    return upper_slices 

#########information about code##########
# create list (This is directly from PackMem)
## v1 starting point 
## v2 ending point
## unit normal: tells the magnitude and the direction of the increment 
## output centroids 
def create_list_descend(upper_limit, end_point, unit_normal, step_size = .1):
    centroids_along_normal_vector = []
    #create numpy array of the start_point so vector math can be complete
    upper_limit_np = np.array(upper_limit)
    end_point_np = np.array(end_point)
    i = upper_limit_np.copy()

    ## need to use a distance measurement 

    while np.dot(end_point_np - i, unit_normal) > 0: ##0: represents the point at which i is == end_point ##when dot product is negative it means you passed the end_point
        centroids_along_normal_vector.append(i)
        i = i + unit_normal * step_size

    return centroids_along_normal_vector

def atm_lower_limit_assignment(pdblines,leaflet_residues_new, atm_normal_vectors, atm_centroids, f_mid):
    atm_lower_limit = {}
    #goal
        ## to take each atom
        ## then take its vector associated with the atom
        ## then compute the intersection between the normal vector and the f_mid 
    
    #first set the normal vectors coming out from each centroid that is corresponds to 
    #then see where the two intersect
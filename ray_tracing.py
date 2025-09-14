import numpy as np
from scipy.spatial import cKDTree

def neighborhood_project_atoms_to_triangles(atom_features, centroids, normals, cutoff, leaflet):
    """ 
    Projects atoms in a given leaflet onto triangle planes.

    Create local neighborhoods by only storing atoms that are near the triangle centroid.
    
    atom_featuers: dict
        atm_number: atm_name, leaflet, atm_coords, lipid_atm_name

    centroids: dict
        triangle_number: centroid coordinates(x,y,z)
    
    normals: dict
        triangle_number: normal vector (x,y,z)


      """
    
    #extract and label features of the atom 
    atom_numbers = []
    atom_coords = []

    for atom_num, feat in atom_features.items(): #atom_name (key) feat(dictionary value)
        if feat[1] == leaflet: #checks if the atom is in the upper or lower leaflet
            atom_numbers.append(atom_num)
            atom_coords.append(feat[2])
        #turn into numpy arrays
        atom_numbers_np = np.array(atom_numbers)
        atom_coords_np = np.array(atom_coords)

    
    #KDTree for fast neighbor search 
    tree = cKDTree(atom_coords)

    
    centroids_only = []
    projections = []
    for tri_idx, centroid in centroids.items(): #stores the index of centroids as tri_idx and the value as centroid
        centroids_only.append(centroid) #stores the centroid coordinates 
        
        normal = normals[tri_idx]

        #find only atoms close to this triangle centroid 
        nearby_idx = tree.query_ball_point(centroid, cutoff) #centroid is the point wehere you want to check neighbors around 

        tri_projections = []

        for j in nearby_idx:
            atom_coord = atom_coords[j] #stores the atom coords for each j (atom_num that is near a triangle)

        
            #projection of atom onto the same plane of the triangle
            w_vector = np.subtract(centroid, atom_coord) #vector from centroid of triangle to atom coord 
            scalar = np.dot(w_vector,normal)
            n_contribution_in_w_vector = np.multiply(scalar, normal)
            vector_orthogonal_normal = np.subtract(w_vector, n_contribution_in_w_vector)

            tri_projections.append(vector_orthogonal_normal)
    
        projections.append(np.array(tri_projections))
    
    centroids_only_np = np.array(centroids_only, dtype = float)



    #checking to see if the atoms and the centroids span the same approximate space
    print("Atom x min/max:", np.min(atom_coords_np[:,0]), np.max(atom_coords_np[:,0]))
    print("Centroid x min/max:", np.min(centroids_only_np[:,0]), np.max(centroids_only_np[:,0]))
    
    return projections









#determines the coords for each atom in a specific triangle 
import numpy as np 

def get_projection(atom_num, grid_id, neighborhood_projections, neighborhood_grid_id, atom_features):
    """
    Returns the (x,y,z) projection for a specific atom in a specific triangle.
    """
    #take atom coords 
    atom_coords = atom_features[atom_num][2]

    
    # Find the index of the atom number in the array
    # np.where returns a tuple of arrays; [0][0] grabs the first match
    atom_idx = np.where(atom_ids == grid_id)[0][0]
    
    # Use the index to get the projected coordinates
    atom_coords_on_plane = neighborhood_projections[atom_num][atom_idx]
    
    return atom_coords_on_plane
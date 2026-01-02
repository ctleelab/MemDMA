import numpy as np
import neighborhood_grids as ng


# get the index of the smallest value (by using the function np.argmin)
# the output of this is the same for all keys
# grid_dict: grid with it's neighborhood
# i
def grid_assignment(
    gl2_coords: dict,
    centroid_triangles,
    tree_upper,
    grid_keys_upper: list,
    centroids_upper,
    cutoff: int,
    atom_features: dict,
):
    """
    Assigns a grid to a GL2 bead then computes the grid neighborhood using cKDTree.

    Args:
        gl2_coords (dict): A dictionary containing GL2 bead coordinates. {GL2 bead number: (x,y,z)}
        centroid_triangles (np.ndarray): A matrix of all centroids of the grid, where each index is the grid number.
        tree_upper (cKDTree): A cKDTree object for efficient spatial queries.
        grid_keys_upper (list): A list of grid keys corresponding to the grid matrix.
        centroids_upper (np.ndarray): An array of centroid coordinates for the upper grid.
        cutoff (int): The cutoff distance for neighborhood determination.
        atom_features (dict):{atom: (atom_name), (leaflet_assignment), (atom_coords), (residue number)}
    Returns:
        gl2_grid_assignment (dict): A dictionary with GL2 number (key): (closet grid_index, grid_coords of closest grid_index, neighborhood_list of the grid_idex).
        lipid_grid_assignment (dict): A dictionary with lipid (key): (neighborhood_list of the grid_idex)
    """

    gl2_grid_assignment = {}
    lipid_grid_assignment = {}

    for key, gl2_points in gl2_coords.items():

        # compute the euclidean distance between each grid and GL2 bead
        # each row represent one grid point
        distance = np.linalg.norm((centroid_triangles - gl2_points), axis=1)
        # return the smallest value of the distance array index value, which will be the grid closest to the GL2 bead
        min_distance_index = np.argmin(distance)

        # takes the index of the closest grid to the GL2 bead and gets their centroid coordinates
        grid_coords = centroid_triangles[min_distance_index]

        # determine the list of neighborhood grids numbers for the GL2 bead, based on the closest grid centroid coordinates
        list_neighborhood = ng.neighborhood(
            cutoff, min_distance_index, tree_upper, grid_keys_upper, centroids_upper
        )

        # add the GL2 bead number to the new dictionary with the following values (index of closest grid, centroid of closest grid, list of neighborhood grid numbers)
        gl2_grid_assignment[key] = (min_distance_index, grid_coords, list_neighborhood)

        # determine the lipid_residue that the GL2 bead belongs to
        lipid_res = atom_features[key][4]

        # creates a dictionary {lipid_res: neighborhood grid list based on GL2 bead in residue}
        lipid_grid_assignment[lipid_res] = list_neighborhood

    return gl2_grid_assignment, lipid_grid_assignment


# gl2_neighborhood_grids: {gl2_atom_number: grid_number, (coords of grid), list of neighborhood grids}
# all_atom: {atom_number: resnumber, coords}
# gl2_beads_resnumber: {gl2_atom_number: resnumber}


def atom_grid_mapping(
    all_atom: dict, gl2_neighborhood_grids: dict, gl2_beads_resnumber: dict
):
    """
    Append a list of neighborhood grids to the all_atoms dictionary.
    Args:
        all_atom (dict): {atom number: res_number, coords}
        gl2_neighborhood_grids (dict): {A dictionary with GL2 atom number: (closet grid_index, grid_coords of closest grid_index, neighborhood_list of the grid_idex)}
        gl2_beads_resnumber (dict): {atom number of gl2 bead: res_number}
    Returns:
        all_atom (dict): {atom number: res_number, coords, list of neighborhood grid number}
    """

    # Precompute reverse lookup. Mapping residue number -> gl2_atom
    res_to_gl2_atom = {v: str(k).strip() for k, v in gl2_beads_resnumber.items()}

    for key1, v1 in all_atom.items():

        res_number = v1[0]  # pulls out residue number

        # Find corresponding GL2 atom for each residue number
        gl2_atom = res_to_gl2_atom.get(res_number)

        if (
            gl2_atom is None
        ):  # Error handeling of residues that are not present in GL2 mapping
            print("Error, residue not found in the Gl2 mapping.")
            break

        values_1 = gl2_neighborhood_grids.get(
            gl2_atom
        )  # pulls out the values of the gl2_neighborhood_grids dictionary

        neighborhood_list = values_1[
            2
        ]  # pulls out the neighborhood grid list of the gl2_neighborhood_grids dictionary

        # Append neighborhood list to atom entry and preserves the original data
        all_atom[key1] = (*v1, neighborhood_list)

    return all_atom

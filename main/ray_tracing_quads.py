# goal: to map each atom to the same plane as the grids
# take the atom coordinates, then project each atom onto the neighborhood of grids, then record it's point

# get atom number
# get atom coords
# get list of close by grids
# get center of close by grids
# project the atom onto the plane of the grids

import numpy as np


def project_block_of_triangles(
    quadrilateral_normal_upper, atom_grid_assignment_upper, upper_centroids
):
    block_projections = {}

    for atom_num, atom_info in atom_grid_assignment_upper.items():
        if (
            len(atom_info) < 3
        ):  # deals with dictionary entries that don't have three items
            # print(f"⚠️ Skipping atom {atom_num}: malformed entry {atom_info}")
            continue

        atom_coords_np = np.array(atom_info[1])  # pull out atom coordinates
        nearby_grids = atom_info[2]  # list of nearby grids

        for grid in nearby_grids:

            coords = np.array(
                upper_centroids[grid]
            )  # grabs the coordinate of the centroid of each grid

            # vector from atom to centroid
            w_vectors = coords - atom_coords_np

            # Projection of w_vectors onto the normal
            scalars = np.dot(w_vectors, quadrilateral_normal_upper[grid])  # shape (M,)
            n_contributions = np.outer(
                scalars, quadrilateral_normal_upper[grid]
            )  # shape (M, 3)

            # Remove normal component → get in-plane vectors
            projections = w_vectors - n_contributions  # shape (M, 3)

            block_projections[grid] = projections

    return block_projections

import numpy as np
from scipy.spatial import cKDTree
from concurrent.futures import ProcessPoolExecutor, as_completed

# from concurrent.futures import ThreadPoolExecutor, as_completed
from tqdm import tqdm
from multiprocessing import shared_memory
import pickle

# blocks should be the atoms
# shared memory should be the grid coordinates


# ###No longer need####
# # Extract atoms coords in the given leaflet
def extract_atoms(atom_features, leaflet):
    atom_num = [
        atom_num for atom_num, feat in atom_features.items() if feat[1] == leaflet
    ]
    atom_coords = [feat[2] for feat in atom_features.values() if feat[1] == leaflet]
    atom_coords_np = np.array(atom_coords)
    return atom_num, atom_coords_np


# Extract the grid coords in the given leaflet, used to create the shared memory
# input: centroid_triangle dictionary {grid_num: (x,y,z)}
def extract_grid_coords(centroid_triangles):
    grid_num = [grid_num for grid_num in centroid_triangles.keys()]
    grid_coords = list(centroid_triangles.values())
    grid_coords_np = np.array(grid_coords)
    return grid_num, grid_coords_np


# Creates a shared memory of the grid coordinates array
def create_shared_atom_coords(arr, name=None):
    shm = shared_memory.SharedMemory(create=True, size=arr.nbytes, name=name)
    shared_arr = np.ndarray(
        arr.shape, dtype=arr.dtype, buffer=shm.buf
    )  # turn into numpy array   ##do I need this?
    shared_arr[:] = arr[:]  # copy data into shared memory
    meta = {
        "name": shm.name,
        "shape": arr.shape,
        "dtype": arr.dtype,
    }  # small dictionary containing information about the memory block
    return shm, shared_arr, meta


# Function that actually executes the projection of the atoms on the same plane as the triangle
# will want to import also the neighborhood
# Need to change this function: where only loop through the neighborhood list


def _project_block_of_atoms(
    block_atoms, normals, grid_num, coords_meta, atom_grid_neighborhood_assignment
):
    block_projections = {}
    block_atom_numbers = {}

    shm = shared_memory.SharedMemory(name=coords_meta["name"])
    grid_coords_np = np.ndarray(
        coords_meta["shape"], dtype=coords_meta["dtype"], buffer=shm.buf
    )

    for atom in block_atoms:
        values = atom_grid_neighborhood_assignment.get(atom)
        if not values or len(values) != 3:
            print(f"⚠️ Skipping malformed or missing entry for atom {atom}: {values}")
            continue

        atom_resid = values[0]
        atom_coords = values[1]
        nearby_idx_list = values[2]

        if atom not in block_projections:
            block_projections[atom] = []

        if atom not in block_atom_numbers:
            block_atom_numbers[atom] = []

        for tri_indx in nearby_idx_list:
            # print(f"grid_coords_np shape: {grid_coords_np.shape}")
            tri_indx -= 1
            # if tri_indx not in normals:
            #     print(f"KeyError: tri_indx {tri_indx} not found in normals")
            #     continue
            # if tri_indx >= len(grid_coords_np):
            #     print(f"Invalid index: {tri_indx}, max allowed: {len(grid_coords_np) - 1}")
            #     continue

            grid_coord = grid_coords_np[tri_indx]
            w_vectors = (
                grid_coord - atom_coords
            )  # change the naming of the grid_coords to triangle_centroid_coords

            # scalars = np.dot(w_vectors, normals[tri_indx+1])
            # n_contributions = np.outer(scalars, normals[tri_indx+1])

            # Calculate the normal direction by computing the dot product of the vector from the triangle centroid to the atom and the normalized triangle normal
            n_contributions = (
                np.dot(w_vectors, normals[tri_indx + 1]) * normals[tri_indx + 1]
            )

            # Project the atom onto the triangle plane by removing the normal contribution
            projections = w_vectors - n_contributions

            block_projections[atom].append(projections.tolist())
            block_atom_numbers[atom].append(int(tri_indx))  # use atom ID directly

    shm.close()
    return block_projections, block_atom_numbers


# max_workers tells you how many threads are being used at once
def neighborhood_project_atoms_to_triangles(
    atom_num,
    grid_num,
    coords_np,
    normals,
    atom_grid_neighborhood_assignment,
    max_workers=24,
    block_size=1000,
):
    """
    Projects atoms in a given leaflet onto triangle planes, in parallel.
    """
    # Load shared memory
    shm_coords, shared_coords, coords_meta = create_shared_atom_coords(coords_np)

    # Split atoms into blocks
    atom_indices = list(atom_num)
    blocks = [
        atom_indices[i : i + block_size]
        for i in range(0, len(atom_indices), block_size)
    ]

    # Dictionaries to store results
    projections_dict = {}
    atom_numbers_dict = {}

    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = {
            executor.submit(
                _project_block_of_atoms,
                block,
                normals,
                grid_num,
                coords_meta,
                atom_grid_neighborhood_assignment,
            ): block
            for block in blocks
        }

        for future in tqdm(
            as_completed(futures), total=len(futures), desc="Processing triangle blocks"
        ):
            block_projections, block_atom_numbers = future.result()

            # Merge block results into global dict
            for atom in block_projections:
                if atom not in projections_dict:
                    projections_dict[atom] = []
                    atom_numbers_dict[atom] = []
                projections_dict[atom].extend(block_projections[atom])
                atom_numbers_dict[atom].extend(block_atom_numbers[atom])

    shm_coords.close()
    shm_coords.unlink()

    # with open("neighborhood_projections_upper_test.pickle", "wb") as f:
    #     pickle.dump(projections_dict, f)

    # with open("neighborhood_projections_upper_atom_id_test.pickle", "wb") as f:
    #     pickle.dump(atom_numbers_dict, f)

    return projections_dict, atom_numbers_dict

    # # Keep triangles sorted
    # sorted_projections = [projections_list[tri_idx] for tri_idx in sorted(centroids)]  # give coordinates of the atoms
    # sorted_atom_numbers = [atom_numbers_list[tri_idx] for tri_idx in sorted(centroids)]  # gives the atom number


# #max_workers tells you how many threads are being used at once
# ##uses the CKDTree to find nearby atoms
# def neighborhood_project_atoms_to_triangles_tree(atom_num, atom_coords_np, centroids, normals, tree, cutoff, max_workers=24, block_size=1000):
#     """
#     Projects atoms in a given leaflet onto triangle planes.
#     """
#     #load in shared arrays
#     shm_coords, shared_coords, coords_meta = create_shared_atom_coords(atom_coords_np)

#     # Parallelization over blocks of triangles
#     num_triangles = len(centroids)
#     projections_list = [None] * num_triangles
#     atom_numbers_list = [None] * num_triangles

#     # Split triangle indices into blocks
#     triangle_indices = list(centroids.keys())
#     blocks = [triangle_indices[i:i + block_size] for i in range(0, len(triangle_indices), block_size)]

#     with ProcessPoolExecutor(max_workers=max_workers) as executor:
#         futures = {
#             executor.submit(
#                 _project_block_of_triangles,
#                 block,
#                 centroids,
#                 normals,
#                 atom_num,
#                 coords_meta,
#                 tree,
#                 cutoff
#             ): block
#             for block in blocks
#         }
#         for future in tqdm((as_completed(futures)), total=len(futures), desc="Processing triangle blocks"):
#             block_projections, block_atom_numbers = future.result()
#             for tri_idx in block_projections:
#                 projections_list[tri_idx] = block_projections[tri_idx]
#                 atom_numbers_list[tri_idx] = block_atom_numbers[tri_idx]

#     shm_coords.close(); shm_coords.unlink()

#     # # Keep triangles sorted
#     # sorted_projections = [projections_list[tri_idx] for tri_idx in sorted(centroids)]  # give coordinates of the atoms
#     # sorted_atom_numbers = [atom_numbers_list[tri_idx] for tri_idx in sorted(centroids)]  # gives the atom number
#     return projections_list, atom_numbers_list

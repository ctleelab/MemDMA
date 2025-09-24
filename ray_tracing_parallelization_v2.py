import numpy as np
from scipy.spatial import cKDTree
#from concurrent.futures import ProcessPoolExecutor, as_completed
from concurrent.futures import ThreadPoolExecutor, as_completed
from tqdm import tqdm


# Worker function, which will be executed per triangle
# This determines the normal vector for each triangle

#builds a KDtree of atom positions 
def tree_cKDTree(atom_features, leaflet):
    # Extract atoms coords in the given leaflet
    atom_num = [atom_num for atom_num, feat in atom_features.items() if feat[1] == leaflet]
    atom_coords = [feat[2] for feat in atom_features.values() if feat[1] == leaflet]
    atom_coords_np = np.array(atom_coords)
    tree = cKDTree(atom_coords_np)
    return tree, atom_num, atom_coords_np

def _project_single_triangle(tri_idx, centroid, normal, atom_num, atom_coords_np, tree, cutoff):
    nearby_idx = tree.query_ball_point(centroid, cutoff)
    tri_projections = []
    tri_atom_numbers = []
    for j in nearby_idx:
        atom_coord = atom_coords_np[j]
        atom_id = atom_num[j]

        # projection of atom onto the same plane of the triangle
        w_vector = centroid - atom_coord
        scalar = np.dot(w_vector, normal)
        n_contribution = scalar * normal
        vector_orthogonal_normal = w_vector - n_contribution

        tri_projections.append(vector_orthogonal_normal)
        tri_atom_numbers.append(atom_id)

    return tri_idx, np.array(tri_projections), np.array(tri_atom_numbers)

#max_workers tells you how many threads are being used at once 
def neighborhood_project_atoms_to_triangles(atom_num, atom_coords_np, centroids, normals, tree, cutoff, max_workers=24):
    """
    Projects atoms in a given leaflet onto triangle planes, in parallel.
    """
    # Parallelization over triangles
    projections_dict = {}
    atoms_dict = {}
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {
            executor.submit(
                _project_single_triangle,
                tri_idx,
                centroids[tri_idx],
                normals[tri_idx],
                   # pass full array
                atom_num,
                atom_coords_np, 
                tree,
                cutoff
            ): tri_idx
            for tri_idx in centroids
        }
        for future in tqdm((as_completed(futures)), total = len(futures), desc="Processing triangles") :
            tri_idx, tri_proj, tri_atoms = future.result()
            projections_dict[tri_idx] = tri_proj
            atoms_dict[tri_idx] = tri_atoms

    # Keep triangles sorted
    #sorted sorts the keys of the dictionary, thereby sorting by triangle number [0] = triangle 1 
    sorted_projections = [projections_dict[tri_idx] for tri_idx in sorted(centroids)] #give coordinates of the atoms
    sorted_atom_numbers = [atoms_dict[tri_idx] for tri_idx in sorted(centroids)] #gives the atom number 
    return sorted_projections, sorted_atom_numbers

    ## Optional: convert centroids to numpy array for checks
    # centroids_only_np = np.array([centroids[tri_idx] for tri_idx in sorted(centroids)], dtype=float)
    # print("Atom x min/max:", np.min(atom_coords_np[:,0]), np.max(atom_coords_np[:,0]))
    # print("Centroid x min/max:", np.min(centroids_only_np[:,0]), np.max(centroids_only_np[:,0]))


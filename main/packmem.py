#!/usr/bin/env python3

# general libraries
import argparse
import sys
import numpy as np
from pathlib import Path

import MDAnalysis as mda
from MDAnalysis.transformations.translate import center_in_box, translate

# library for showing progress bar for runs
from tqdm.auto import tqdm

# library that does parallel processing
from tqdm.contrib.concurrent import process_map

from scipy.spatial import cKDTree

import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

import quadrilateral_mesh
import util

import ctleelab_plothelper.plothelpers as ph


def main(
    u: mda.Universe,
    prop_dict: dict,
    neighborhood_cutoff: float = 10,
    deep_defect_cutoff: float = 1,
    mesh_size: float = 1,
    normal_cutoff: float = 3,
    interval: int = 1000,
    output_path: Path = Path("."),
    prefix: str = "output",
    debug=False,
):

    if debug:
        with mda.Writer("tmp.gro") as fd:
            fd.write(u)

    # +10 or -10 to eliminate create a 10A perimeter of mesh to eliminate false positive defects due to pbc
    x_range_adj = (
        np.floor_divide(
            np.min(u.atoms.positions[:, 0] + neighborhood_cutoff), mesh_size
        )
        * mesh_size,
        np.floor_divide(
            np.max(u.atoms.positions[:, 0] - neighborhood_cutoff), mesh_size
        )
        * mesh_size,
    )
    y_range_adj = (
        np.floor_divide(
            np.min(u.atoms.positions[:, 1]) + neighborhood_cutoff, mesh_size
        )
        * mesh_size,
        np.floor_divide(
            np.max(u.atoms.positions[:, 1]) - neighborhood_cutoff, mesh_size
        )
        * mesh_size,
    )

    # # Create bin edges
    x_edges = np.arange(x_range_adj[0], x_range_adj[1], mesh_size)
    y_edges = np.arange(y_range_adj[0], y_range_adj[1], mesh_size)

    x_centers = 0.5 * (x_edges[:-1] + x_edges[1:])
    y_centers = 0.5 * (y_edges[:-1] + y_edges[1:])
    X_c, Y_c = np.meshgrid(x_centers, y_centers, indexing="ij")

    X, Y = np.mgrid[
        x_range_adj[0] : x_range_adj[1] : mesh_size,
        y_range_adj[0] : y_range_adj[1] : mesh_size,
    ]

    num_frames = u.trajectory.n_frames // interval + 1

    nx, ny = X.shape
    n_tri = 2 * (nx - 1) * (ny - 1)  # two triangles per quad per leaflet

    defect_count_per_frame_upper = np.empty((num_frames, n_tri))
    # defect_count_per_frame_lower = np.empty(
    #     (num_frames, n_tri)
    # )

    # create a numpy array of the triangle centroids
    triangle_centroids_per_frame_upper = np.empty((num_frames, n_tri, 3))
    # triangle_centroids_per_frame_lower = np.empty(
    #     (num_frames, n_tri, 3)
    # )

    triangle_faces = np.empty((n_tri, 3))
    triangle_vertices_per_frame_upper = np.empty((num_frames, len(X.ravel()), 3))
    # create a numpy array of the triangle vertices
    # triangle_vertices_per_frame_lower = np.empty(
    #     (num_frames, n_tri, 3)
    # )

    index_of_frame = 0

    for ts in u.trajectory[::interval]:
        print(index_of_frame * interval)

        print("Start: meshing")
        z_surface, leaflets = quadrilateral_mesh.interpolate_cubic_surface(
            u, "name PO4", X, Y
        )
        print("Finished: meshing")

        if debug:
            with plt.style.context(
                ["ctleelab_plothelper.base", "ctleelab_plothelper.transparent"]
            ):
                fig, ax = plt.subplots(
                    1, 1, figsize=(2, 2), subplot_kw={"projection": "3d"}
                )

                po4 = u.select_atoms("name PO4 NC3 GL1 GL2")
                x = po4.positions[:, 0]
                y = po4.positions[:, 1]
                z = po4.positions[:, 2]

                # mask = np.where(z > np.mean(z))

                ax.scatter(
                    x,
                    y,
                    z,
                    color="blue",
                    s=4,
                    alpha=1,
                    linewidth=0,
                )

                skip = 2
                ax.scatter(
                    X[::skip, ::skip],
                    Y[::skip, ::skip],
                    z_surface["upper"][::skip, ::skip],
                    color="red",
                    s=1,
                    marker=".",
                    alpha=1,
                    linewidth=0,
                )

                # ax.view_init(elev=90, azim=-90, roll=0)

                ax.set_aspect("equal")
                ax.set_xlabel("X (A)")
                ax.set_ylabel("Y (A)")
                ax.set_zlabel("Z (A)")

                ax.xaxis.set_major_locator(MultipleLocator(25))

                fig.savefig("3d-tmp.png")
                # plt.show()
                plt.close(fig)

        print("Start:computing normal vectors of triangle")
        (mesh_faces, mesh_vertices, mesh_centroids, mesh_normals) = (
            util.vectorized_normal_vector_triangles(X, Y, z_surface)
        )

        triangle_centroids_per_frame_upper[index_of_frame] = mesh_centroids["upper"]
        # triangle_centroids_per_frame_lower[index_of_frame] = mesh_centroids["lower"]

        # create a numpy array of the triangle faces
        triangle_vertices_per_frame_upper[index_of_frame] = mesh_vertices["upper"]
        # triangle_vertices_per_frame_lower[index_of_frame] = mesh_vertices["lower"]
        print("Finished:computing normal vectors of triangle")

        print("Start:GL2 grid neighborhood assignment")
        #######################################
        ###create a cKDtree of the centroids###
        #######################################
        tree = {}
        tree["upper"] = cKDTree(mesh_centroids["upper"])
        tree["lower"] = cKDTree(mesh_centroids["lower"])

        # allocate memory for new normal array (N,3)
        smeared_mesh_normals = {}
        smeared_mesh_normals["upper"] = np.empty_like(mesh_normals["upper"])
        smeared_mesh_normals["lower"] = np.empty_like(mesh_normals["lower"])

        for leaflet in ["upper", "lower"]:
            for i in np.arange(0, len(mesh_normals[leaflet])):

                # Aggregate triangles in the neighborhood using KD-tree lookup
                neighborhood_indices = tree[leaflet].query_ball_point(
                    mesh_centroids[leaflet][i],
                    normal_cutoff,
                )

                # Average normals and convert into unit length
                average_normal = np.mean(
                    mesh_normals[leaflet][neighborhood_indices], axis=0
                )
                average_normal /= np.linalg.norm(average_normal)

                smeared_mesh_normals[leaflet][i] = average_normal

        ## To determine which grid is the closest to the GL2 bead, the euclidean distance between the Gl2 bead and the grids centroid was calculated

        defect_counts = {
            "upper": np.zeros(len(mesh_centroids["upper"])),
            "lower": np.zeros(len(mesh_centroids["lower"])),
        }

        def increment_defect_count(leaflet, grid_idx, aliphatic):
            if aliphatic == "n":
                defect_counts[leaflet][grid_idx] += 1
            else:
                defect_counts[leaflet][grid_idx] += 0.001

        # for leaflet in ["upper", "lower"]:
        for leaflet in ["upper"]:
            print(f"Processing {leaflet} leaflet")

            headgroup_beads = leaflets[leaflet].select_atoms("name PO4")
            i, j = np.mgrid[
                0 : len(headgroup_beads),
                0 : len(mesh_centroids[leaflet]),
            ]

            distances = np.linalg.norm(
                mesh_centroids[leaflet][j] - headgroup_beads.positions[i], axis=2
            )
            min_distance_indices = np.argmin(distances, axis=1)

            for i, atom in enumerate(headgroup_beads):
                neighborhood_indices = tree[leaflet].query_ball_point(
                    mesh_centroids[leaflet][min_distance_indices[i]],
                    neighborhood_cutoff,
                )

                r = atom.residue

                center = mesh_centroids[leaflet][min_distance_indices[i]]

                for grid_idx in neighborhood_indices:
                    tri = mesh_faces[grid_idx]
                    vertices = mesh_vertices[leaflet][tri]
                    # normal = mesh_normals[leaflet][grid_idx]
                    normal = smeared_mesh_normals[leaflet][grid_idx]
                    centroid = mesh_centroids[leaflet][grid_idx]

                    for atom in r.atoms:
                        key = f"{r.resname} {atom.name}"
                        if key not in prop_dict:
                            sys.stderr.write(
                                f"WARNING: Atom type {key} not found in property file."
                            )
                            break
                        radius, aliphatic = prop_dict[key]

                        w = atom.position - centroid
                        w_n = np.dot(w, normal)  # component of w normal to the plane

                        if w_n + radius < -deep_defect_cutoff:
                            # print(w_n + radius)
                            # print(atom.position, centroid, normal)
                            continue  # atom is too far from the plane

                        # print(np.linalg.norm(normal), w_n, key)

                        atom_proj = (
                            atom.position - w_n * normal
                        )  # W projected onto plane

                        vertex_to_atom_proj = atom_proj - vertices

                        # increment if atom overlaps with a vertex
                        if (
                            np.min(np.linalg.norm(vertex_to_atom_proj, axis=1))
                            <= radius
                        ):
                            increment_defect_count(leaflet, grid_idx, aliphatic)
                        else:
                            # print('checking edges', key)
                            # Check edges
                            edges = np.vstack(
                                (
                                    vertices[1] - vertices[0],
                                    vertices[2] - vertices[1],
                                    vertices[0] - vertices[2],
                                )
                            )
                            edge_lengths = np.linalg.norm(edges, axis=1)

                            edge_dots = np.einsum(
                                "ij,ij->i", vertex_to_atom_proj, edges
                            )

                            edge_projection_valid = np.bitwise_and(
                                (edge_dots >= 0), (edge_dots <= edge_lengths)
                            )

                            dist_to_edge = (
                                np.linalg.norm(
                                    np.cross(vertex_to_atom_proj, edges), axis=1
                                )
                                / edge_lengths
                            )

                            if np.any((dist_to_edge <= radius) & edge_projection_valid):
                                # print("edge valid")
                                increment_defect_count(leaflet, grid_idx, aliphatic)

        defect_count_per_frame_upper[index_of_frame] = defect_counts["upper"]
        # defect_count_per_frame_lower[index_of_frame] = defect_counts["lower"]

        index_of_frame += 1

    output_path.mkdir(parents=True, exist_ok=True)

    # save defect count and grid centroid coordinates as a .npy
    np.save(
        output_path / f"{prefix}-defect-count-per-frame-upper.npy",
        defect_count_per_frame_upper,
    )
    np.save(
        output_path / f"{prefix}-triangle-centroids-per-frame-upper.npy",
        triangle_centroids_per_frame_upper,
    )

    np.save(
        output_path / f"{prefix}-triangle-faces.npy",
        np.squeeze(mesh_faces),
    )

    np.save(
        output_path / f"{prefix}-triangle-vertices-per-frame-upper.npy",
        triangle_vertices_per_frame_upper,
    )


##this is the interface with the user, so that they can specify inputs
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Arguments for the app")

    parser.add_argument("-i", action="store", dest="gro_filename", help="gro file")
    parser.add_argument(
        "-t",
        action="store",
        dest="traj_filename",
        help="Trajectory file (.xtc)",
    )

    parser.add_argument(
        "-r",
        action="store",
        dest="prop_file",
        help="File for radius (default vdw_radiiFinal2014.txt)",
        default="vdw_radiiFinal2014.txt",
    )
    parser.add_argument(
        "-o",
        action="store",
        dest="output_path",
        default="output",
        help="Path for outputs",
    )
    parser.add_argument(
        "-p",
        action="store",
        dest="prefix",
        default="output",
        help="Prefix for output files",
    )
    parser.add_argument(
        "-d",
        action="store",
        dest="deep_defect_cutoff",
        default=1,
        type=float,
        help="Distance cutoff for deep defects in Angstrom",
    )
    # parser.add_argument(
    #     "-n",
    #     action="store",
    #     dest="indexFile",
    #     help="File for index (Gromacs ndx style). Only Lower/Upper group accepted",
    # )
    parser.add_argument(
        "-v", dest="verbosity", action="store_true", help="Increase the verbosity"
    )
    # parser.add_argument("-f", action="store", dest="frame_num", help="frame number")
    parser.add_argument(
        "-m",
        action="store",
        dest="mesh_size",
        default=1,
        type=float,
        help="size of the mesh in Angstrom",
    )
    parser.add_argument(
        "-nc",
        action="store",
        dest="neighborhood_cutoff",
        default=10,
        help="sets the size of the neighborhood cutoff for grids in angstrom.",
    )
    parser.add_argument(
        "--normal-cutoff",
        action="store",
        dest="normal_cutoff",
        default=3,
        help="sets the size of the normal smearing cutoff for grids in angstrom.",
    )

    args = parser.parse_args()

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    ##error handling if the chosen distance is negative. This distance is key when determining where each bead falls within the mesh grid.
    if args.deep_defect_cutoff < 0.0:
        print("ERROR : The distance for the -d option must be > 0")
        sys.exit()

    prop_dict = util.parse_radius_aliphatic_dict(args.prop_file)

    if args.traj_filename:
        u = mda.Universe(args.gro_filename, args.traj_filename, continous=True)
    else:
        u = mda.Universe(args.gro_filename, continous=True)

    debug = False
    if args.verbosity:
        debug = True

    main(
        u,
        prop_dict,
        neighborhood_cutoff=args.neighborhood_cutoff,
        deep_defect_cutoff=args.deep_defect_cutoff,
        mesh_size=args.mesh_size,
        output_path=Path(args.output_path),
        prefix=args.prefix,
        normal_cutoff=args.normal_cutoff,
        debug=False,
    )

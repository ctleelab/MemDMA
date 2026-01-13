#!/usr/bin/env python3

# general libraries
import argparse
import sys
import numpy as np
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


##this is the interface with the user, so that they can specify inputs
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Arguments for the app")
    ##changed to take .gro files
    parser.add_argument(
        "-i",
        action="store",
        dest="traj_filename",
        help="Trajectory file (e.g., .gro, .xtc)",
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
        dest="outputname",
        default="output",
        help="Name for output file",
    )
    parser.add_argument(
        "-d",
        action="store",
        dest="deep_defect_cutoff",
        default=1,
        type=float,
        help="Distance cutoff for deep defects in Angstrom",
    )
    parser.add_argument(
        "-n",
        action="store",
        dest="indexFile",
        help="File for index (Gromacs ndx style). Only Lower/Upper group accepted",
    )
    parser.add_argument(
        "-v", dest="pdbout", action="store_true", help="Increase the verbosity"
    )
    parser.add_argument("-f", action="store", dest="frame_num", help="frame number")
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

    args = parser.parse_args()

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    ##error handling if the chosen distance is negative. This distance is key when determining where each bead falls within the mesh grid.
    if args.deep_defect_cutoff < 0.0:
        print("ERROR : The distance for the -d option must be > 0")
        sys.exit()

    prop_dict = util.parse_radius_aliphatic_dict(args.prop_file)

    u = mda.Universe(args.traj_filename, continous=True)
    u.atoms.positions += [0, 0, 50]
    u.atoms.wrap()

    with mda.Writer("tmp.gro") as fd:
        fd.write(u)

    print("Start: meshing")
    # X, Y, zupper, zlower, leaflet_residues = quadrilateral_mesh.mesh(
    #     args.traj_filename, "name GL2", args.mesh_size
    # )
    X, Y, z_surface, leaflets, X_c, Y_c = quadrilateral_mesh.interpolate_cubic_surface(
        u, selection="name PO4 NC3", step=args.mesh_size
    )
    # X, Y, z_surface, leaflets, X_c, Y_c = quadrilateral_mesh.get_flat_surface(
    #     u, selection="name PO4 NC3", step=args.mesh_size
    # )
    print("Finished: meshing")

    with plt.style.context(
        ["ctleelab_plothelper.base", "ctleelab_plothelper.transparent"]
    ):
        fig, ax = plt.subplots(1, 1, figsize=(2, 2), subplot_kw={"projection": "3d"})

        po4 = u.select_atoms("name PO4 NC3 GL1 GL2")
        x = po4.positions[:, 0]
        y = po4.positions[:, 1]
        z = po4.positions[:, 2]

        mask = np.where(z > np.mean(z))

        ax.scatter(
            x[mask],
            y[mask],
            z[mask],
            color="blue",
            s=4,
            alpha=1,
            linewidth=0,
        )
        ax.set_aspect("equal")

        skip = 2
        # ax.plot_surface(
        #     X[::skip, ::skip],
        #     Y[::skip, ::skip],
        #     z_surface["upper"][::skip, ::skip],
        #     color="red",
        #     alpha=0.8,
        # )
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

        # ax.plot_surface(
        #     X[::skip, ::skip],
        #     Y[::skip, ::skip],
        #     z_surface["lower"][::skip, ::skip],
        #     color="red",
        #     alpha=0.8,
        # )
        # ax.scatter(
        #     X[::skip, ::skip],
        #     Y[::skip, ::skip],
        #     z_surface["lower"][::skip, ::skip],
        #     color="red",
        #     s=1,
        #     marker=".",
        #     alpha=1,
        #     linewidth=0,
        # )

        # ax.view_init(elev=90, azim=-90, roll=0)

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
    print("Finished:computing normal vectors of triangle")

    print("Start:GL2 grid neighborhood assignment")

    #######################################
    ###create a cKDtree of the centroids###
    #######################################
    tree = {}
    tree["upper"] = cKDTree(mesh_centroids["upper"])
    tree["lower"] = cKDTree(mesh_centroids["lower"])

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
                args.neighborhood_cutoff,
            )

            r = atom.residue

            average_normal = np.mean(
                mesh_normals[leaflet][neighborhood_indices], axis=0
            )
            average_normal /= np.linalg.norm(average_normal)
            # print(average_normal)

            for grid_idx in neighborhood_indices:
                tri = mesh_faces[grid_idx]
                vertices = mesh_vertices[leaflet][tri]
                # normal = mesh_normals[leaflet][grid_idx]
                normal = average_normal
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

                    if w_n + radius < -args.deep_defect_cutoff:
                        # print(w_n + radius)
                        # print(atom.position, centroid, normal)
                        continue  # atom is too far from the plane

                    # print(np.linalg.norm(normal), w_n, key)

                    atom_proj = atom.position - w_n * normal  # W projected onto plane

                    vertex_to_atom_proj = atom_proj - vertices

                    # increment if atom overlaps with a vertex
                    if np.min(np.linalg.norm(vertex_to_atom_proj, axis=1)) <= radius:
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

                        edge_dots = np.einsum("ij,ij->i", vertex_to_atom_proj, edges)

                        edge_projection_valid = np.bitwise_and(
                            (edge_dots >= 0), (edge_dots <= edge_lengths)
                        )

                        dist_to_edge = (
                            np.linalg.norm(np.cross(vertex_to_atom_proj, edges), axis=1)
                            / edge_lengths
                        )

                        if np.any((dist_to_edge <= radius) & edge_projection_valid):
                            # print("edge valid")
                            increment_defect_count(leaflet, grid_idx, aliphatic)

    mask = np.where(defect_counts["upper"] >= 1)
    deep_mask = np.where(defect_counts["upper"] == 0)
    shallow_mask = np.where((defect_counts["upper"] < 1) & (defect_counts["upper"] > 0))

    fig, ax = plt.subplots(figsize=(8, 6))

    ax.set_xlim(0, np.max(mesh_centroids["upper"][:, 0]))
    ax.set_ylim(0, np.max(mesh_centroids["upper"][:, 1]))

    po4 = u.select_atoms("name PO4 NC3 GL1 GL2")
    x = po4.positions[:, 0]
    y = po4.positions[:, 1]
    z = po4.positions[:, 2]
    mask = np.where(z > np.mean(z))
    ax.scatter(x[mask], y[mask], s=60, color="blue")

    ax.set_aspect("equal", adjustable="box")
    fig.savefig("upper-defect-tmp_polaronly.png")
    # sc = ax.scatter(
    #     mesh_centroids["upper"][:, 0][mask],
    #     mesh_centroids["upper"][:, 1][mask],
    #     c=defect_counts["upper"][mask],
    #     cmap="viridis",
    #     s=1,
    # )

    # sc = ax.scatter(
    #     mesh_centroids["upper"][:, 0][mask],
    #     mesh_centroids["upper"][:, 1][mask],
    #     c="black",
    #     s=1,
    # )
    # fig.colorbar(sc, ax=ax, label="Atom count")

    print(defect_counts["upper"][shallow_mask])

    ax.scatter(
        mesh_centroids["upper"][:, 0][shallow_mask],
        mesh_centroids["upper"][:, 1][shallow_mask],
        color="orange",
        s=5,
        alpha=0.5,
        label="shallow",
    )

    ax.scatter(
        mesh_centroids["upper"][:, 0][deep_mask],
        mesh_centroids["upper"][:, 1][deep_mask],
        color="red",
        s=5,
        alpha=0.5,
        label="deep",
    )

    ax.set_aspect("equal", adjustable="box")

    ax.xaxis.set_major_locator(MultipleLocator(25))

    ax.set_xlabel("X (A)")
    ax.set_ylabel("Y (A)")
    ax.set_title(f"Triangle Centroids with Atom Counts")

    plt.tight_layout()
    fig.savefig("upper-defect-tmp.png")

    # np.savetxt(
    #     (
    #         mesh_centroids["upper"][:, 0],
    #         mesh_centroids["upper"][:, 1],
    #         z_surface["upper"],
    #         defect_counts["upper"],
    #     ),
    #     delimiter=",",
    # )

    # # # creates a graph that plots the atom centroid (x,y) color with its count value
    # # neigh_graph.atom_mapping_centroids_2d_zvalue(
    # #     centroid_triangles_u, lipid_type_short, args.frame_num
    # # )
    # # neigh_graph.atom_mapping_centroids_2d(
    # #     centroid_triangles_u,
    # #     triangle_atom_assignment_upper_sorted,
    # #     lipid_type_short,
    # #     args.frame_num,
    # # )

    # # also creates a csv file of (x,y,z, atom count)
    # neigh_graph.combined_atom_mappin_plot(
    #     centroid_triangles_upper,
    #     triangle_atom_assignment_upper_sorted,
    #     lipid_type_short,
    #     args.frame_num,
    # )

    # neigh_graph.atom_mapping_centroids_2d_matplotlib(
    #     centroid_triangles_upper,
    #     triangle_atom_assignment_upper_sorted,
    #     lipid_type_short,
    #     args.frame_num,
    # )

    # print(f"completed frame{args.frame_num}")

#!/usr/bin/env python3

# general libraries
import argparse
import sys
import numpy as np
from scipy import ndimage as scipy_ndimage
from pathlib import Path

import MDAnalysis as mda
from MDAnalysis.transformations.translate import center_in_box, translate

# library for showing progress bar for runs
from tqdm.auto import tqdm

# library that does parallel processing
from tqdm.contrib.concurrent import process_map

from scipy.spatial import cKDTree

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.ticker import MultipleLocator

import quadrilateral_mesh
import util

import ctleelab_plothelper.plothelpers as ph


def _process_residue_mapping(args):
    (
        r,
        vertices,
        normal,
        centroid,
        prop_dict,
        deep_defect_cutoff,
        grid_idx,
        leaflet,
    ) = args
    increment = 0

    def increment_defect_count(aliphatic):
        nonlocal increment
        if aliphatic == "n":
            increment += 1
        else:
            increment += 0.001

    for atom in r.atoms:
        key = f"{r.resname} {atom.name}"
        if key not in prop_dict:
            sys.stderr.write(f"WARNING: Atom type {key} not found in property file.")
            break
        radius, aliphatic = prop_dict[key]

        w = atom.position - centroid
        w_n = np.dot(w, normal)  # component of w normal to the plane

        # w_n is a distance relative to the normal orientation
        if w_n + radius < -deep_defect_cutoff:
            continue  # atom is too far from the plane

        atom_proj = atom.position - w_n * normal  # W projected onto plane

        vertex_to_atom_proj = atom_proj - vertices

        # increment if atom overlaps with a vertex
        if np.min(np.linalg.norm(vertex_to_atom_proj, axis=1)) <= radius:
            increment_defect_count(aliphatic)
        else:
            # print('checking edges', key)
            # Check edges
            edges = np.vstack(
                (
                    vertices[1] - vertices[0],
                    vertices[2] - vertices[1],
                    vertices[3] - vertices[2],
                    vertices[0] - vertices[3],
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
                increment_defect_count(aliphatic)
    return (leaflet, grid_idx, increment)


def main(
    u: mda.Universe,
    prop_dict: dict,
    neighborhood_cutoff: float = 10,
    deep_defect_cutoff: float = 1,
    mesh_size: float = 1,
    normal_cutoff: float = 3,
    interval: int = 1,
    output_path: Path = Path("."),
    prefix: str = "output",
    debug=False,
):
    # if debug:
    #     with mda.Writer("tmp.gro") as fd:
    #         fd.write(u)

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
    n_quad = (nx - 1) * (ny - 1)
    quad_shape = (nx - 1, ny - 1)
    n_tri = 2 * (nx - 1) * (ny - 1)  # two triangles per quad per leaflet

    defects_per_frame_upper = np.empty((num_frames, n_quad))
    defects_per_frame_lower = np.empty((num_frames, n_quad))

    # create a numpy array of the triangle centroids
    centroids_per_frame_upper = np.empty((num_frames, n_quad, 3))
    centroids_per_frame_lower = np.empty((num_frames, n_quad, 3))

    vertices_per_frame_upper = np.empty((num_frames, len(X.ravel()), 3))
    vertices_per_frame_lower = np.empty((num_frames, len(X.ravel()), 3))

    index_of_frame = 0
    for ts in tqdm(u.trajectory[::interval], position=0, desc="Processing frames"):
        # print("Start: meshing")
        Z, mesh_faces, mesh_vertices, mesh_normals, mesh_centroids, leaflets = (
            quadrilateral_mesh.interpolate_cubic_surface(u, "name PO4", X, Y)
        )
        # print("Finished: meshing")

        if debug:
            with plt.style.context(
                ["ctleelab_plothelper.base", "ctleelab_plothelper.transparent"]
            ):
                fig, ax = plt.subplots(
                    1, 1, figsize=(2, 2), subplot_kw={"projection": "3d"}
                )

                po4 = leaflets["upper"].select_atoms("name PO4 NC3 GL1 GL2")

                ax.scatter(
                    po4.positions[:, 0],
                    po4.positions[:, 1],
                    po4.positions[:, 2],
                    color="blue",
                    s=4,
                    alpha=1,
                    linewidth=0,
                )

                skip = 2
                ax.scatter(
                    X[::skip, ::skip],
                    Y[::skip, ::skip],
                    Z["upper"][::skip, ::skip],
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

        # mesh_faces, mesh_vertices, mesh_centroids, mesh_normals = (
        #     util.vectorized_normal_vector_triangles(X, Y, z_surface)
        # )

        centroids_per_frame_upper[index_of_frame] = mesh_centroids["upper"]
        centroids_per_frame_lower[index_of_frame] = mesh_centroids["lower"]

        # create a numpy array of the triangle faces
        vertices_per_frame_upper[index_of_frame] = mesh_vertices["upper"]
        vertices_per_frame_lower[index_of_frame] = mesh_vertices["lower"]

        # print("Start:GL2 grid neighborhood assignment")
        #######################################
        ###create a cKDtree of the centroids###
        #######################################
        tree = {}
        tree["upper"] = cKDTree(mesh_centroids["upper"])
        tree["lower"] = cKDTree(mesh_centroids["lower"])

        # allocate memory for new normal array (N,3)
        smeared_normals = {}
        smeared_normals["upper"] = np.empty_like(mesh_normals["upper"])
        smeared_normals["lower"] = np.empty_like(mesh_normals["lower"])

        if debug:
            for cutoff in np.arange(0, 10):
                tmp = []
                tmp_n = []
                for i in np.arange(0, len(mesh_normals["upper"])):
                    neighborhood_indices = tree["upper"].query_ball_point(
                        mesh_centroids["upper"][i],
                        cutoff,
                    )

                    # Average normals and convert into unit length
                    average_normal = np.mean(
                        mesh_normals["upper"][neighborhood_indices], axis=0
                    )
                    average_normal /= np.linalg.norm(average_normal)
                    _smeared_normals = average_normal
                    res = np.inner(_smeared_normals, np.array([0, 0, 1]))
                    tmp.append(res)
                    tmp_n.append(len(neighborhood_indices))

                print(
                    f"Normal Cutoff: {cutoff} Mean/min/max alignment: {np.mean(tmp)}, {np.min(tmp)}, {np.max(tmp)}, Avg neighbors: {np.mean(tmp_n)}"
                )

        if normal_cutoff > 0:
            for leaflet in ["upper", "lower"]:
                for i in np.arange(0, len(mesh_normals[leaflet])):
                    # Aggregate faces in the neighborhood using KD-tree lookup
                    neighborhood_indices = tree[leaflet].query_ball_point(
                        mesh_centroids[leaflet][i],
                        normal_cutoff,
                    )

                    # Average normals and convert into unit length
                    average_normal = np.mean(
                        mesh_normals[leaflet][neighborhood_indices], axis=0
                    )
                    average_normal /= np.linalg.norm(average_normal)

                    smeared_normals[leaflet][i] = average_normal
        else:
            smeared_normals["upper"] = mesh_normals["upper"]
            smeared_normals["lower"] = mesh_normals["lower"]

        ## To determine which grid is the closest to the GL2 bead, the euclidean distance between the Gl2 bead and the grids centroid was calculated
        defect_counts = {
            "upper": np.zeros(len(mesh_centroids["upper"])),
            "lower": np.zeros(len(mesh_centroids["lower"])),
        }

        for leaflet in ["upper", "lower"]:
            # print(f"Processing {leaflet} leaflet")

            headgroup_beads = leaflets[leaflet].select_atoms("name PO4")
            i, j = np.mgrid[
                0 : len(headgroup_beads),
                0 : len(mesh_centroids[leaflet]),
            ]

            distances = np.linalg.norm(
                mesh_centroids[leaflet][j] - headgroup_beads.positions[i], axis=2
            )
            min_distance_indices = np.argmin(distances, axis=1)

            jobs = []
            for i, atom in enumerate(headgroup_beads):
                neighborhood_indices = tree[leaflet].query_ball_point(
                    mesh_centroids[leaflet][min_distance_indices[i]],
                    neighborhood_cutoff,
                )
                r = atom.residue
                # center = centroids[leaflet][min_distance_indices[i]]

                for grid_idx in neighborhood_indices:
                    quad = mesh_faces[grid_idx]
                    verts = mesh_vertices[leaflet][quad]
                    # normal = normals[leaflet][grid_idx]
                    normal = smeared_normals[leaflet][grid_idx]
                    # normal = np.array([0, 0, 1])
                    centroid = mesh_centroids[leaflet][grid_idx]

                    jobs.append(
                        (
                            r,
                            verts,
                            normal,
                            centroid,
                            prop_dict,
                            deep_defect_cutoff,
                            grid_idx,
                            leaflet,
                        )
                    )

            r = process_map(
                _process_residue_mapping,
                jobs,
                max_workers=8,
                chunksize=1000,
                position=1,
                leave=False,
                desc="Mapping residues to grids",
            )
            # print(r)

            for leaflet, grid_idx, increment in r:
                if increment > 0:
                    defect_counts[leaflet][grid_idx] += increment

        defects_per_frame_upper[index_of_frame] = defect_counts["upper"]

        print(defect_counts["upper"].shape)

        for leaflet in ["upper", "lower"]:

            fig, ax = ph.fixed_size_subplots(
                1, 1, subwidth=4, subheight=1, rmargin_scale=1.5
            )

            reshaped_defects = defect_counts[leaflet].reshape(quad_shape)
            # masked_polar_defects = np.ma.masked_where(
            #     reshaped_defects < 1, reshaped_defects
            # )

            # sc = ax.pcolor(X, Y, masked_polar_defects, cmap="viridis", shading="auto")

            masked_shallow_defects = np.ma.masked_where(
                (reshaped_defects >= 1) | (reshaped_defects == 0),
                reshaped_defects,
            )
            cmap_single_orange = mcolors.ListedColormap(
                ["orange"], name="constant_orange"
            )
            ax.pcolor(
                X, Y, masked_shallow_defects, cmap=cmap_single_orange, shading="auto"
            )

            masked_deep_defects = np.ma.masked_where(
                reshaped_defects != 0,
                reshaped_defects,
            )
            cmap_single_red = mcolors.ListedColormap(["red"], name="constant_red")
            ax.pcolor(X, Y, masked_deep_defects, cmap=cmap_single_red, shading="auto")

            # polar_positions = (
            #     leaflets[leaflet].select_atoms("name PO4 NC3 GL1 GL2").positions
            # )
            # x = polar_positions[:, 0]
            # y = polar_positions[:, 1]
            # z = polar_positions[:, 2]
            # ax.scatter(x, y, s=4, color="pink")
            output_path.mkdir(parents=True, exist_ok=True)

            ax.set_aspect("equal", adjustable="datalim")
            ax.xaxis.set_major_locator(MultipleLocator(25))
            ax.set_xlabel("X (A)")
            ax.set_ylabel("Y (A)")
            fig.savefig(output_path / f"{index_of_frame}-{leaflet}-defect-tmp.png")

            plt.close(fig)

        # defect_count_per_frame_lower[index_of_frame] = defect_counts["lower"]

        index_of_frame += 1

    # save defect count and grid centroid coordinates as a .npy
    np.save(
        output_path / f"{prefix}-defects-upper.npy",
        defects_per_frame_upper,
    )
    np.save(
        output_path / f"{prefix}-centroids.npy",
        centroids_per_frame_upper,
    )

    np.save(
        output_path / f"{prefix}-faces.npy",
        np.squeeze(mesh_faces),
    )

    np.save(
        output_path / f"{prefix}-vertices-upper.npy",
        vertices_per_frame_upper,
    )
    np.save(
        output_path / f"{prefix}-shape.npy",
        quad_shape,
    )

    np.save(
        output_path / f"{prefix}-Xshape.npy",
        X,
    )
    np.save(
        output_path / f"{prefix}-Yshape.npy",
        Y,
    )
    for leaflet in ["upper", "lower"]:
        fig, ax = ph.fixed_size_subplots(
            1, 1, subwidth=4, subheight=1, rmargin_scale=1.5
        )

        reshaped_defects = defect_counts[leaflet].reshape(quad_shape)
        # masked_polar_defects = np.ma.masked_where(
        #     reshaped_defects < 1, reshaped_defects
        # # )

        # sc = ax.pcolor(X, Y, masked_polar_defects, cmap="viridis", shading="auto")

        masked_shallow_defects = np.ma.masked_where(
            (reshaped_defects >= 1) | (reshaped_defects == 0),
            reshaped_defects,
        )
        cmap_single_orange = mcolors.ListedColormap(["orange"], name="constant_orange")
        ax.pcolor(X, Y, masked_shallow_defects, cmap=cmap_single_orange, shading="auto")

        masked_deep_defects = np.ma.masked_where(
            reshaped_defects != 0,
            reshaped_defects,
        )
        cmap_single_red = mcolors.ListedColormap(["red"], name="constant_red")
        ax.pcolor(X, Y, masked_deep_defects, cmap=cmap_single_red, shading="auto")

        # polar_positions = (
        #     leaflets[leaflet].select_atoms("name PO4 NC3 GL1 GL2").positions
        # )
        # x = polar_positions[:, 0]
        # y = polar_positions[:, 1]
        # z = polar_positions[:, 2]
        # ax.scatter(x, y, s=4, color="pink")

        # ph.add_fixed_colorbar(sc, ax=ax, pad=0.05, aspect=30, label="Polar atom count")

        ax.set_aspect("equal", adjustable="datalim")
        ax.xaxis.set_major_locator(MultipleLocator(25))
        ax.set_xlabel("X (A)")
        ax.set_ylabel("Y (A)")
        fig.savefig(f"{leaflet}-defect-tmp.png")

        plt.close(fig)

        # GENERATE BINARY DEFECT MAP
        binary_defect_map = np.where(reshaped_defects < 1, 1, 0)
        # plot binary defect map
        fig, ax = ph.fixed_size_subplots(
            1, 1, subwidth=4, subheight=1, rmargin_scale=1.5
        )
        sc = ax.pcolor(X, Y, binary_defect_map, cmap="viridis", shading="auto")
        ax.set_aspect("equal", adjustable="datalim")
        ax.xaxis.set_major_locator(MultipleLocator(25))
        ax.set_xlabel("X (A)")
        ax.set_ylabel("Y (A)")
        fig.savefig(f"{leaflet}-binary-defect-tmp.png")
        plt.close(fig)

        # LABEL CONNECTED COMPONENTS IN BINARY DEFECT MAP
        label, nfeat = scipy_ndimage.label(binary_defect_map)

        # plot labeled defects
        fig, ax = ph.fixed_size_subplots(
            1, 1, subwidth=4, subheight=1, rmargin_scale=1.5
        )
        sc = ax.pcolor(X, Y, label, cmap="viridis", shading="auto")
        ax.set_aspect("equal", adjustable="datalim")
        ax.xaxis.set_major_locator(MultipleLocator(25))
        ax.set_xlabel("X (A)")
        ax.set_ylabel("Y (A)")
        fig.savefig(f"{leaflet}-labeled-defect-tmp.png")
        plt.close(fig)

        defect_sizes = []
        for i in range(1, nfeat + 1):
            mask = np.where(label == i)
            defect_sizes.append(np.count_nonzero(mask))
        
        fig, ax = ph.fixed_size_subplots(1,1, subwidth=3, subheight=3)

        ax.hist(defect_sizes, bins=50)

        ax.set_xlabel("Defect size (A^2)")
        ax.set_ylabel("Count")
        fig.savefig(f"{leaflet}-defect-size-dist.png")
        plt.close(fig)



    fig, ax = ph.fixed_size_subplots(1, 1, subwidth=2, subheight=2)

    po4 = u.select_atoms("name PO4").positions[:, 2]
    nc3 = u.select_atoms("name NC3").positions[:, 2]
    gl = u.select_atoms("name GL1 GL2").positions[:, 2]
    other = u.select_atoms("not (name PO4 NC3 GL1 GL2)").positions[:, 2]

    ax.hist(
        [po4, nc3, gl, other],
        bins=50,
        stacked=True,
        label=["PO4", "NC3", "GL1/GL2", "Other"],
        color=["yellow", "blue", "green", "grey"],
        density=True,
        alpha=0.2,
    )

    ax.set_xlabel("Z Position (A)")
    ax.set_ylabel("Counts")
    ax.xaxis.set_major_locator(MultipleLocator(5))
    fig.savefig(output_path / "density.png")
    plt.close(fig)

    # # Determine the average distance between PO4 and GL1 for each lipid
    # distances = []
    # for lipid in u.select_atoms("resname DPPC").residues:
    #     po4 = lipid.atoms.select_atoms("name PO4")
    #     gl1 = lipid.atoms.select_atoms("name GL1")
    #     if len(po4) > 0 and len(gl1) > 0:
    #         distance = np.linalg.norm(po4.center_of_geometry() - gl1.center_of_geometry())
    #         distances.append(distance)

    # if distances:
    #     average_distance = np.mean(distances)
    #     print(f"Average distance between PO4 and GL1: {average_distance:.2f} Å")
    # else:
    #     print("No PO4 or GL1 atoms found for the specified lipid.")
    # print(average_distances)


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
        help="File for radius parameters",
        default="vdw_radii_CGMartini3_openbeta2_updated.txt",
    )
    parser.add_argument(
        "-o",
        action="store",
        dest="output_path",
        default="output",
        help="Output paths",
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
        "--deep-defect-cutoff",
        action="store",
        dest="deep_defect_cutoff",
        default=10,
        type=float,
        help="Distance cutoff for deep defects in Angstrom",
    )

    parser.add_argument(
        "-v",
        "--verbose",
        dest="verbosity",
        action="count",
        help="Increase the verbosity",
    )
    parser.add_argument(
        "-m",
        "--mesh-size",
        action="store",
        dest="mesh_size",
        default=1,
        type=float,
        help="size of the mesh in Angstrom",
    )
    parser.add_argument(
        "--neighborhood-cutoff",
        action="store",
        dest="neighborhood_cutoff",
        default=10,
        type=float,
        help="sets the size of the neighborhood cutoff for grids in angstrom.",
    )
    parser.add_argument(
        "--normal-cutoff",
        action="store",
        dest="normal_cutoff",
        default=0,
        type=float,
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
        u = mda.Universe(args.gro_filename, args.traj_filename)
    else:
        u = mda.Universe(args.gro_filename)

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
        debug=debug,
    )

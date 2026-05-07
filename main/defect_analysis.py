#!/usr/bin/env python3
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, NullFormatter
import MDAnalysis as mda
import numpy as np
from pathlib import Path
import ctleelab_plothelper.plothelpers as ph
import matplotlib.colors as mcolors
import matplotlib.gridspec as gridspec
from scipy import ndimage as scipy_ndimage
from scipy.stats import linregress
from tqdm.auto import tqdm
import argparse
import os
import pandas as pd

defect_type = ["deep", "shallow", "all"]


# pixi inputs
parser = argparse.ArgumentParser(description="Membrane curvature and defect analysis")

parser.add_argument(
    "-s", required=True, action="store", dest="system", help="system name"
)
parser.add_argument(
    "-r",
    type=int,
    default=200,
    dest="restriction",
    help="Restriction value for where to analyze defects (in Angstroms from the center of the membrane)",
)
parser.add_argument(
    "-i", required=True, action="store", dest="input", help="Input folder"
)
parser.add_argument(
    "-o", required=True, action="store", dest="output", help="Output folder"
)
parser.add_argument(
    "-ds",
    type=int,
    default=15,
    action="store",
    dest="defect_size",
    help="Lower limit for defect size to be included in analysis (in A^2)",
)

parser.add_argument(
    "-l", required=True, action="store", dest="leaflet", help="Output folder"
)



args = parser.parse_args()

# access variables
system = args.system
restriction = args.restriction
input_folder = Path(args.input)
analysis_folder = args.output
lower_limit_size_defect = args.defect_size
leaflet = args.leaflet

os.makedirs(analysis_folder, exist_ok=True)


# no cutoff of defect size used, to look at the spatial distribution of all defects across the membrane, and how it correlates with height of the membrane.
def surface_coverage_defects_nocutoff(input_folder, leaflet):
    defect_counts = np.load(input_folder / f"output-defects-{leaflet}.npy")
    quad_centroids_np = np.load(input_folder / f"output-centroids-{leaflet}.npy")
    shape = np.load(input_folder / "output-shape.npy")
    X = np.load(input_folder / "output-Xshape.npy")
    Y = np.load(input_folder / "output-Yshape.npy")

    n_frames = len(quad_centroids_np)
    cutsize = 8

    deep_defect_counts_per_frame_along_x = np.empty((n_frames, shape[0]))
    shallow_defect_counts_per_frame_along_x = np.empty((n_frames, shape[0]))
    all_defect_counts_per_frame_along_x = np.empty((n_frames, shape[0]))

    deep_defect_counts_per_frame_per_grid = np.empty(
        (n_frames, shape[0], shape[1] - cutsize * 2)
    )
    shallow_defect_counts_per_frame_per_grid = np.empty(
        (n_frames, shape[0], shape[1] - cutsize * 2)
    )
    all_defect_counts_per_frame_per_grid = np.empty(
        (n_frames, shape[0], shape[1] - cutsize * 2)
    )

    z_coord_per_frame_along_x = np.empty((n_frames, shape[0]))

    for frame in tqdm(np.arange(0, n_frames)):
        defect_counts_frame = defect_counts[frame]
        centroid_frame = quad_centroids_np[frame]

        z_coord_centroid_frame = centroid_frame[:, 2]
        # Add a value of 1 to the grid where the Boolean is true
        masked_shallow_defects_num = np.where(
            (defect_counts_frame > 0) & (defect_counts_frame < 1),
            1,
            0,
        )

        # Add a value of 1 to the grid where the Boolean is true
        masked_deep_defects_num = np.where(
            (defect_counts_frame == 0),
            1,
            0,
        )

        # all defect
        masked_all_defects_num = np.where(defect_counts_frame < 1, 1, 0)

        def cut_out_y_artifacts(arr, begin, end):
            out = np.hstack(
                (
                    arr[:, :begin],
                    arr[
                        :,
                        end:,
                    ],
                )
            )
            deleted = arr[:, begin:end]
            return out, deleted

        # reshape mask to represent shape of grids
        # gets a value per grid
        reshaped_shallow_defects = masked_shallow_defects_num.reshape(shape)
        reshaped_deep_defects = masked_deep_defects_num.reshape(shape)
        reshaped_all_defects = masked_all_defects_num.reshape(shape)

        # take out any defects with area < 15A
        # here handle the defect based on size

        def filter_defects_by_size(binary_mask, lower_limit_size_defect):
            label_mask, nfeat = scipy_ndimage.label(binary_mask)
            filtered_mask = np.zeros_like(binary_mask, dtype=int)

            for i in range(1, nfeat + 1):
                defect_size = np.count_nonzero(label_mask == i)
                if defect_size >= lower_limit_size_defect:
                    filtered_mask[label_mask == i] = 1

            return filtered_mask, label_mask, nfeat

        # filter each defect type before any averaging
        reshaped_shallow_defects, label_shallow, nfeat_shallow = filter_defects_by_size(
            reshaped_shallow_defects, lower_limit_size_defect
        )

        reshaped_deep_defects, label_deep, nfeat_deep = filter_defects_by_size(
            reshaped_deep_defects, lower_limit_size_defect
        )

        reshaped_all_defects, label_all, nfeat_all = filter_defects_by_size(
            reshaped_all_defects, lower_limit_size_defect
        )

        center = X.shape[1] // 2

        min_index = center - cutsize
        max_index = center + cutsize
        # print(shape, center, min_index, max_index)

        reshaped_shallow_defects, _ = cut_out_y_artifacts(
            reshaped_shallow_defects, min_index, max_index
        )
        reshaped_deep_defects, _ = cut_out_y_artifacts(
            reshaped_deep_defects, min_index, max_index
        )
        reshaped_all_defects, _ = cut_out_y_artifacts(
            reshaped_all_defects, min_index, max_index
        )

        # new shape of array
        new_shape = reshaped_shallow_defects.shape

        # sum down y, to get a count per grid in the x-dimension
        shallow_defect_counts_per_frame_along_x[frame] = np.sum(
            reshaped_shallow_defects, axis=1
        )
        deep_defect_counts_per_frame_along_x[frame] = np.sum(
            reshaped_deep_defects, axis=1
        )
        all_defect_counts_per_frame_along_x[frame] = np.sum(
            reshaped_all_defects, axis=1
        )

        # save the grid value per frame
        shallow_defect_counts_per_frame_per_grid[frame] = reshaped_shallow_defects
        deep_defect_counts_per_frame_per_grid[frame] = reshaped_deep_defects
        all_defect_counts_per_frame_per_grid[frame] = reshaped_all_defects

        # height calculations
        reshaped_z_coord = z_coord_centroid_frame.reshape(shape)
        z_coord_per_frame_along_x[frame] = np.sum(reshaped_z_coord, axis=1)

    # sum of counts at a particular x grid across all frames
    sum_deep_defect_across_frames_along_x = np.sum(
        deep_defect_counts_per_frame_along_x, axis=0
    )
    sum_shallow_defect_across_frames_along_x = np.sum(
        shallow_defect_counts_per_frame_along_x, axis=0
    )
    sum_all_defect_across_frames_along_x = np.sum(
        all_defect_counts_per_frame_along_x, axis=0
    )

    sum_height_across_frames_along_x = np.sum(z_coord_per_frame_along_x, axis=0)

    # #divide by the total number of grids (number of grids in y-dim * number of frames)
    total_grids = new_shape[1] * n_frames

    # #percentage of defects in particular region of the membrane
    systems_dict_deep_spatial_weight_along_x = np.divide(
        sum_deep_defect_across_frames_along_x, total_grids
    )
    systems_dict_shallow_spatial_weight_along_x = np.divide(
        sum_shallow_defect_across_frames_along_x, total_grids
    )
    systems_dict_all_spatial_weight_along_x = np.divide(
        sum_all_defect_across_frames_along_x, total_grids
    )

    # handeling to get average z-coord per x
    avg_z_coord_across_frames = np.divide(sum_height_across_frames_along_x, total_grids)

    # calculate the mean value of defect count add to system dict
    systems_dict_avg_defect_shallow_per_grid = np.mean(
        shallow_defect_counts_per_frame_per_grid, axis=0
    )
    systems_dict_avg_defect_deep_per_grid = np.mean(
        deep_defect_counts_per_frame_per_grid, axis=0
    )
    systems_dict_avg_defect_all_per_grid = np.mean(
        all_defect_counts_per_frame_per_grid, axis=0
    )

    min_avg_z_coord_across_frames = np.min(avg_z_coord_across_frames)
    max_avg_z_coord_across_frames = np.max(avg_z_coord_across_frames)

    avg_z_coord_across_frames_normalized = (
        avg_z_coord_across_frames - min_avg_z_coord_across_frames
    ) / 10  # convert to nm

    systems_dict_height = avg_z_coord_across_frames_normalized

    kept_Y, deleted_Y = cut_out_y_artifacts(Y, min_index, max_index)
    kept_X, deleted_X = cut_out_y_artifacts(X, min_index, max_index)

    # need to adjust the values in the Y and X to account for the deleted rows of data
    adjusted_kept_Y = []
    for i in range(0, len(kept_Y[1])):
        number = int(kept_Y[1][i])
        if number > deleted_Y[0][0]:
            temp = int(number - cutsize * 2)
            adjusted_kept_Y.append(temp)
        else:
            adjusted_kept_Y.append(number)

    # need to create the mesh with the cutout y-values
    x_axis = kept_Y.shape[0]
    X, Y = np.meshgrid(np.arange(x_axis), adjusted_kept_Y, indexing="ij")

    return (
        systems_dict_height,
        systems_dict_shallow_spatial_weight_along_x,
        systems_dict_deep_spatial_weight_along_x,
        systems_dict_all_spatial_weight_along_x,
    )


####################################################################
########frames that represent the average defect coverage  #########
####################################################################

def frame_extraction_avg_surface_defect_coverage(
    defect_restricted_avg_system, n_frames, percent_coverage_defects_per_frame
):
    """
    Compiles a list of frames that is +/- .01 away from the avg defect coverage over the entire membrane

    Args:
    defect_restricted_avg_system: list [avgerage defect size, std]
    n_frames: int number of frames in the trajectory
    percent_coverage_defects_per_frame: list of average defect coverage per frame in the trajectory

    Returns:
    max_frames_avg_defect_coverage_per_system (int): A single frame that represents the avg coverage of defects across the entire membrane for each lipid composition.

    """

    list_frames_close_avg_defect = []
    avg_defect = defect_restricted_avg_system[0]
    upper_bound = avg_defect + .01
    lower_bound = avg_defect - .01
    for frame in np.arange(0, n_frames):
        defect_coverage = percent_coverage_defects_per_frame[frame]
        if (
            lower_bound < defect_coverage < upper_bound
        ):  # storing frames close to the average defect coverage for entire membrane
            list_frames_close_avg_defect.append(frame)

    configurations_types_frames_to_avg_defect = list_frames_close_avg_defect

    max_frames_avg_defect_coverage_per_system  = max(configurations_types_frames_to_avg_defect)
    return max_frames_avg_defect_coverage_per_system


########################################################################################################################################
########combined graph of height of membrane and percentage of defects in particular region of the membrane in the x-dimension #########
########################################################################################################################################
def height_and_defect_probability_combined_graph(
    leaflet, input_folder, analysis_folder
):

    for defect_type_name in defect_type:
        with plt.style.context(
            ["ctleelab_plothelper.base", "ctleelab_plothelper.light"]
        ):
            fig, ax = ph.fixed_size_subplots(
                1, 1, wmargin=1, hmargin=1, subwidth=2, subheight=2, rmargin_scale=4
            )
            # fig, ax = ph.fixed_size_subplots(1, 2, wmargin = 1, hmargin = 1,  subwidth=3, subheight=1, rmargin_scale = 8)
            legend_handles = None
            legend_labels = None

            (
                systems_dict_height,
                systems_dict_shallow_spatial_weight_along_x,
                systems_dict_deep_spatial_weight_along_x,
                systems_dict_all_spatial_weight_along_x,
            ) = surface_coverage_defects_nocutoff(input_folder, leaflet)

            # Load data
            height_membrane = systems_dict_height

            if defect_type_name == "shallow":
                compiled_defect_counts = systems_dict_shallow_spatial_weight_along_x
            elif defect_type_name == "deep":
                compiled_defect_counts = systems_dict_deep_spatial_weight_along_x
            else:
                compiled_defect_counts = systems_dict_all_spatial_weight_along_x

            x_idx = np.arange(len(height_membrane))
            x_idx_centered = (x_idx - (len(height_membrane) / 2)) / 10  # convert to nm

            ax.scatter(
                x_idx_centered, height_membrane, s=12, color="#8100ffff", label="Height"
            )
            ax.set_xlabel("X (nm)")
            ax.set_ylabel("Height (nm)")
            ax.set_xlim(
                -len(height_membrane) / 20, len(height_membrane) / 20
            )  # adjust for nm
            ax.set_ylim(0, 18)  # nm
            ax.set_title(f"{system},  {defect_type_name} defects", fontsize=12)  # nm
            ax.xaxis.set_minor_locator(MultipleLocator(2.5))
            ax.xaxis.set_major_locator(MultipleLocator(5))

            # Secondary y-axis (defect counts)

            ax2 = ax.twinx()
            ax2.bar(
                x_idx_centered,
                compiled_defect_counts,
                width=0.1,
                color="#1f91bc",
                alpha=0.4,
                label="Defect probability",
            )
            ax2.set_ylabel("Defect probability")
            ax2.set_ylim(0, 0.5)

            if legend_handles is None:
                h1, l1 = ax.get_legend_handles_labels()
                h2, l2 = ax2.get_legend_handles_labels()
                legend_handles = h1 + h2
                legend_labels = l1 + l2

        plt.legend(
            legend_handles,
            legend_labels,
            loc="center left",
            bbox_to_anchor=(1.5, 0.5),
            frameon=False,
            fontsize=10,
        )
        plt.subplots_adjust(right=20)

        plt.savefig(
            f"{analysis_folder}/{leaflet}_leaflet_combined_defect_probability_and_height_{defect_type_name}.png"
        )
        # plt.savefig(f"{analysis_folder}/{leaflet}_leaflet_combined_defect_probability_and_height_{defect_type_name}.pdf")
        plt.close()


height_and_defect_probability_combined_graph("upper", input_folder, analysis_folder)


def total_num_defect_restricted_X_cutoff(
    defect_data, shape, centroids, restriction, cutoff
):

    defect_grid = defect_data.reshape(shape)

    # center of grids at 0
    centered_centroids = centroids - np.mean(centroids)

    centroid_X_mask = (
        np.where(
            (centered_centroids >= -restriction) & (centered_centroids <= restriction),
            1,
            0,
        )
    ).squeeze()

    centroid_mask_grid = centroid_X_mask.reshape(shape)

    restricted_region_grids_with_defects = np.where(
        (centroid_mask_grid == 1) & (defect_grid == 1), 1, 0
    )

    # here handle the defect based on size
    all_defect_sizes_cutoff = []
    label_all, nfeat_all = scipy_ndimage.label(restricted_region_grids_with_defects)
    for i in range(1, nfeat_all + 1):
        defect_size = np.count_nonzero(label_all == i)
        if defect_size >= cutoff:
            all_defect_sizes_cutoff.append(defect_size)

    total_grids = len(defect_data)

    defect_grids = np.sum(all_defect_sizes_cutoff)

    percent_defect_coverage = (defect_grids / total_grids) * 100

    # return avg defect size
    defect_avg_size = np.mean(all_defect_sizes_cutoff)

    # return total_num_defects_restricted_X, total_grids_restricted_X
    return defect_grids, total_grids, percent_defect_coverage, defect_avg_size


# restrict the area to be analyzed
# then determine the size of the defects
# then only add the defects that are > cutoff
def surface_coverage_defects_cutoff(
    restriction, cutoff, input_folder, leaflet, analysis_folder
):
    defect_counts = np.load(input_folder / f"output-defects-{leaflet}.npy")
    quad_centroids_np = np.load(input_folder / f"output-centroids-{leaflet}.npy")
    shape = np.load(input_folder / "output-shape.npy")
    X = np.load(input_folder / "output-Xshape.npy")
    Y = np.load(input_folder / "output-Yshape.npy")

    n_frames = len(quad_centroids_np)

    percent_cut_defects = {}
    percent_cut_defects_per_frame = {}
    avg_defect_size_over_frames = {}

    total_defects_restricted_region_per_frame = np.empty(n_frames)
    total_grids_restricted_region_per_frame = np.empty(n_frames)
    percentage_defect_coverage_per_frame = np.empty(n_frames)
    avg_defect_size_per_frame = np.empty(n_frames)

    for frame in np.arange(0, n_frames):

        defect_counts_frame = defect_counts[frame]

        centroids_frame = quad_centroids_np[frame, :, 0:1]

        masked_all_defects_num = np.where(defect_counts_frame < 1, 1, 0)

        total_per_frame, total_grids, percent_defect_coverage, avg_defect_size = (
            total_num_defect_restricted_X_cutoff(
                masked_all_defects_num, shape, centroids_frame, restriction, cutoff
            )
        )

        total_defects_restricted_region_per_frame[frame] = total_per_frame

        total_grids_restricted_region_per_frame[frame] = total_grids

        percentage_defect_coverage_per_frame[frame] = percent_defect_coverage

        # handles frames with no defects
        if np.isnan(avg_defect_size):
            avg_defect_size_per_frame[frame] = 0
        else:
            avg_defect_size_per_frame[frame] = avg_defect_size / 100  # convert to nm^2

    average_defect_coverage = np.mean(percentage_defect_coverage_per_frame)

    # print(f"{composition} avg defect coverage {average_defect_coverage}")
    std_defect_coverage = np.std(percentage_defect_coverage_per_frame)

    average_defect_size = np.mean(avg_defect_size_per_frame)
    average_defect_size_std = np.std(avg_defect_size_per_frame)

    shape = percentage_defect_coverage_per_frame.shape
    num_frames = shape[0]

    percent_cut_defects = [average_defect_coverage, std_defect_coverage, num_frames]
    avg_defect_size_over_frames = [
        average_defect_size,
        average_defect_size_std,
        num_frames,
    ]

    # storing per frame array, to look at convergence
    percent_cut_defects_per_frame = percentage_defect_coverage_per_frame

    # store as a .csv
    rows_percent_cut_defects = []
    avg, std, num_frames = percent_cut_defects
    rows_percent_cut_defects.append(
        {"avg_defect_coverage": avg, "std": std, "num_frames": num_frames}
    )
    df_1 = pd.DataFrame(rows_percent_cut_defects)
    df_1.to_csv(
        f"{analysis_folder}/{leaflet}_leaflet_{system}_defect_avg_coverage{restriction}_cutoff{cutoff}.csv",
        index=False,
    )

    rows_avg_defect_size = []
    avg, std, num_frames = avg_defect_size_over_frames
    rows_avg_defect_size.append(
        {"avg_defect_size (nm²)": avg, "std": std, "num_frames": num_frames}
    )
    df_2 = pd.DataFrame(rows_avg_defect_size)
    df_2.to_csv(
        f"{analysis_folder}/{leaflet}_leaflet_{system}_defect_data_avg_size{restriction}_cutoff{cutoff}.csv",
        index=False,
    )

    return (
        percent_cut_defects,
        avg_defect_size_over_frames,
        percent_cut_defects_per_frame,
        n_frames
    )


(
    defect_coverage_restricted_avg_system_upper,
    defect_size_restricted_avg_system_upper,
    defect_restricted_per_frame_system_upper,
    n_frames
) = surface_coverage_defects_cutoff(
    restriction, lower_limit_size_defect, input_folder, leaflet, analysis_folder
)

# (
#     defect_coverage_restricted_avg_system_lower,
#     defect_size_restricted_avg_system_lower,
#     defect_restricted_per_frame_system_lower,
#     n_frames
# ) = surface_coverage_defects_cutoff(
#     restriction, lower_limit_size_defect, input_folder, "lower", analysis_folder
# )


################################################################
###scatter plot for defect coverage over frames of simulation###
################################################################


def defect_coverage_simulation(
    defect_restricted_per_frame_system, restriction, leaflet
):
    """
    Plots the percent defect coverage for each frame for a given configuration and restriction, to look at convergence of defect coverage over time.
    Can handle two compositions at a time.
    """
    with plt.style.context(["ctleelab_plothelper.base", "ctleelab_plothelper.light"]):

        fig, ax = ph.fixed_size_subplots(1, 1, subwidth=5, subheight=3)

        per_frame_array = defect_restricted_per_frame_system

        frames = np.arange(len(per_frame_array))

        ax.scatter(frames, per_frame_array, s=3)
        ax.set_xlabel("Frame")
        ax.set_ylabel("Percent defect coverage")
        ax.yaxis.set_major_locator(MultipleLocator(2.5))
        ax.set_title(
            f"Defect coverage per frame {leaflet} ({system}, restriction={restriction} Å)"
        )

        ax.legend(frameon=False)

        plt.tight_layout()

        plt.savefig(
            f"{analysis_folder}/{leaflet}_leaflet_frame_defect_coverage_{system}_restriction{restriction}.png",
            dpi=300,
        )


defect_coverage_simulation(
    defect_restricted_per_frame_system_upper, restriction, leaflet
)
# defect_coverage_simulation(
#     defect_restricted_per_frame_system_lower, restriction, "lower"
# )


#################################################################
###plot the probability of finding a defect of a certain size###
##################################################################
def defect_size_probability(input_folder, restriction, lower_limit_size_defect,leaflet):

    major_probs = np.array([1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6])
    defect_counts = np.load(input_folder / f"output-defects-{leaflet}.npy")
    quad_centroids_np = np.load(input_folder / f"output-centroids-{leaflet}.npy")
    n_frames = len(quad_centroids_np)
    shape = np.load(input_folder / "output-shape.npy")

    defect_sizes = []
    bins = np.arange(0, 50000, 1)  # each bin is 1A^2

    for frame in range(n_frames):

        defects_counts_frame = defect_counts[frame].reshape(shape)

        masked_all_defects_num = np.where(defects_counts_frame < 1, 1, 0)

        centroids_frame = quad_centroids_np[frame, :, 0:1]

        centered_centroids = centroids_frame - np.mean(centroids_frame)
        centroid_X_mask = np.where(
            (centered_centroids >= -restriction)
            & (centered_centroids <= restriction),
            1,
            0,
        ).squeeze()

        centroid_mask_grid = centroid_X_mask.reshape(shape)

        restricted_region_grids_with_defects = np.where(
            (centroid_mask_grid == 1) & (masked_all_defects_num == 1), 1, 0
        )

        label_all, nfeat_all = scipy_ndimage.label(
            restricted_region_grids_with_defects
        )

        for i in range(1, nfeat_all + 1):
            defect_size = np.count_nonzero(label_all == i)
            if defect_size >= lower_limit_size_defect:
                defect_sizes.append(defect_size)

        # histogram
        counts, _ = np.histogram(defect_sizes, bins=bins)
        probabilities = counts / counts.sum()
        bin_centers = 0.5 * (bins[1:] + bins[:-1])

        mask = probabilities > 0
        areas_nonzero = bin_centers[mask]
        probabilities_nonzero = probabilities[mask]
        log_probabilities = np.log(probabilities_nonzero)

        # apply PackMem-like fit limits
        masks_limits = (areas_nonzero >= lower_limit_size_defect) & (
            probabilities_nonzero >= 1e-4
        )

        areas_nonzero = areas_nonzero[masks_limits]
        probabilities_nonzero = probabilities_nonzero[masks_limits]
        log_probabilities = np.log(probabilities_nonzero)

        # saves the log probability data per configuration
        strain_defect_sizes_log_transformed_prob = [
            areas_nonzero,
            log_probabilities,
        ]

    with plt.style.context(["ctleelab_plothelper.base", "ctleelab_plothelper.light"]):
        fig, ax = ph.fixed_size_subplots(
            1, 1, subwidth=1.3, subheight=1.4
        )  # small graphs DOPC vs DPPC
        # fig, ax = ph.fixed_size_subplots(1, 1, subwidth=2.2, subheight=1.5) #large graphs +PA, PS, pure
        # colors_presets = [ "#ff6c0A","#8ece33", "#9A36D8"] #large graphs +PA, PS, pure

        areas_nonzero = strain_defect_sizes_log_transformed_prob[0]
        log_probabilities = strain_defect_sizes_log_transformed_prob[1]
        # scatter
        ax.scatter(areas_nonzero, log_probabilities, s=1)

        # axis formatting
        ax.set_yticks(np.log(major_probs))
        ax.set_yticklabels(
            [
                r"$10^{-1}$",
                r"$10^{-2}$",
                r"$10^{-3}$",
                r"$10^{-4}$",
                r"$10^{-5}$",
                r"$10^{-6}$",
            ]
        )

        ax.set_ylabel("Probability")
        ax.set_xlabel("Defect Area (Å²)")
        ax.set_ylim(np.log(5.0e-5), -1)

        ax.set_title(f"{system}")

        plt.legend()

        plt.tight_layout()

        plt.savefig(
            f"{analysis_folder}/{leaflet}_leaflet_{system}_defect_size_probability.png"
        )
        # plt.savefig(
        #     f"{analysis_folder}/reshaped-{system}-defect-size-probability.pdf"
        # )
        plt.close()


defect_size_probability(input_folder, restriction, lower_limit_size_defect,leaflet)




###############################################################
###plot defect map representing the average defect coverage ###
###############################################################

# create defect map plots
def defect_map_avg_defect_coverage(leaflet, input_folder, analysis_folder, restriction, defect_size_cut_off,PackMem_frame_interval,system): 
    with plt.style.context(["ctleelab_plothelper.base", "ctleelab_plothelper.light"]):
        # computes the average defect coverage over the entire simulation
        (
            avg_percent_coverage_defects,
            avg_defect_size_over_frames,
            percent_coverage_defects_per_frame,
            n_frames,

        ) = surface_coverage_defects_cutoff(restriction, defect_size_cut_off, input_folder, leaflet, analysis_folder
        )

        # determines the last frame of the trajectory that has the average defect coverage
        max_frames_avg_defect_coverage_per_system = (
            frame_extraction_avg_surface_defect_coverage(
                avg_percent_coverage_defects,
                n_frames,
                percent_coverage_defects_per_frame,
            )
        )

        fig, ax_array = ph.fixed_size_subplots(
            1,
            1,
            subwidth=1.3,
            subheight=1,
            wmargin=0.5,
            hmargin=0.002,
            rmargin_scale=0.1,
            tmargin_scale=0.2,
        )

        # load data
        defect_counts = np.load(input_folder /f"output-defects-{leaflet}.npy")

        defect_counts_frame = defect_counts[max_frames_avg_defect_coverage_per_system]

        shape = np.load(input_folder/ "output-shape.npy")
        X = (np.load(input_folder / "output-Xshape.npy")) / 10 #convert to nm
        Y = (np.load(input_folder / "output-Yshape.npy")) / 10 #convert to nm

        # center X so midpoint is at 0
        X_center = (X.max() + X.min()) / 2
        X = X - X_center

        Xmin = X.min()
        Xmax = X.max()
        Ymin = Y.min()
        Ymax = Y.max()

        reshaped_defects = defect_counts_frame.reshape(shape)

        # mask for deep and shallow defects
        masked_shallow_defects = (
            np.ma.masked_where(
                (reshaped_defects >= 1) | (reshaped_defects == 0),
                reshaped_defects,
            )
        ) / 10

        masked_deep_defects = (
            np.ma.masked_where(
                reshaped_defects != 0,
                reshaped_defects,
            )
        ) / 10

        # plot
        # cmap_single_orange = mcolors.ListedColormap(["orange"], name="constant_red")
        cmap_single_blue = mcolors.ListedColormap(["#1f91bc"], name="constant_blue")
        ax_array.pcolor(
            X, Y, masked_shallow_defects, cmap=cmap_single_blue, shading="auto"
        )
        ax_array.pcolor(
            X, Y, masked_deep_defects, cmap=cmap_single_blue, shading="auto"
        )

        # formatting
        #ax_array.set_aspect("equal", adjustable="datalim")
        ax_array.xaxis.set_major_locator(MultipleLocator(5))
        ax_array.xaxis.set_minor_locator(MultipleLocator(2.5))
        ax_array.set_xlabel("X (nm)")
        ax_array.set_ylabel("Y (nm)")
        ax_array.set_aspect("equal", adjustable="box")
        ax_array.set_title(system, fontsize=16)  # just use ax, not ax[i,j]

        # save compiled frame
        # plt.savefig(f"{analysis_folder}/{system}-{restriction}-avg-defect-coverage-frame{max_frames_avg_defect_coverage_per_system}.pdf"
        # )
        plt.savefig(
            f"{analysis_folder}/{system}_{leaflet}_restriction{restriction}_avg_defect_coverage_frame{max_frames_avg_defect_coverage_per_system}.png"
        )
        plt.close(fig)


        # MD_frame = max_frames_avg_defect_coverage_per_system*PackMem_frame_interval

        # lipid_configuration = f"{system}-{configuration}"
        # figures_data[lipid_configuration] = [max_frames_avg_defect_coverage_per_system, avg_percent_coverage_defects[0], MD_frame]


    # #export out figures_data as .csv 
    # rows = []
    # for key in figures_data:
    #     frame, average_percent_coverage_defect, MD_image = figures_data[key]
    #     rows.append({
    #         "system": key,
    #         "representative_defect_map_frame": frame,
    #         "avg_defect_coverage": average_percent_coverage_defect, 
    #         "MD_frame": MD_image
    #     })
    # df = pd.DataFrame(rows)
    # df.to_csv(f"{analysis_folder}/reshaped_frames_visualization.csv", index = False)

defect_map_avg_defect_coverage(leaflet, input_folder, analysis_folder , restriction, lower_limit_size_defect,1,system )

# libraries for multiprocessing running this code directly
# from pathlib import Path
# from tqdm.contrib.concurrent import process_map

# Cailyns outline for code
# goal create a mesh from the PO4 beads in each leaflets
# determine the min dimension of the system. Since the system is NVT there should be no change in the box size
# determine the size of the bins
# both in x and y direction
# divide the min dimension by the bin size to get the number of bins
# MembraneSpectralAnalysis

# size of bins 1A or .1nm
# step = .5

# path = Path("/home/casakurai/scratch/packmem/")


def mesh(gro, selection, step):
    import pickle
    import numpy as np
    import MDAnalysis as mda
    from membrane_spectral_analysis.base import MembraneSpectralAnalysis
    import leaflet_finder

    # gro, selection, traj, path = args
    # u = mda.Universe(gro, str(traj) , continous=True)
    u = mda.Universe(gro, continous=True)

    print("units of the file", u.trajectory.units)  # which is nm

    # Average box dimensions over trajectory
    avg_dim_x = np.mean([ts.dimensions[0] for ts in u.trajectory])
    avg_dim_y = np.mean([ts.dimensions[1] for ts in u.trajectory])

    # Determine number of bins
    num_bins_x = int(avg_dim_x / step)
    num_bins_y = int(avg_dim_y / step)

    # Create bin edges
    x_edges = np.linspace(0, avg_dim_x, num_bins_x + 1)
    y_edges = np.linspace(0, avg_dim_y, num_bins_y + 1)

    # Bin centers for plotting
    x_centers = ((0.5 * (x_edges[:-1] + x_edges[1:]))) / 10
    y_centers = ((0.5 * (y_edges[:-1] + y_edges[1:]))) / 10

    X, Y = np.meshgrid(x_centers, y_centers, indexing="ij")

    # Run MembraneSpectralAnalysis to get the quadrilateral mesh
    mc = MembraneSpectralAnalysis(
        u,
        select=selection,
        n_x_bins=num_bins_x,
        n_y_bins=num_bins_y,
        x_range=(0, avg_dim_x),
        y_range=(0, avg_dim_y),
        wrap=True,  # enable wrapping to avoid edge artifacts
        interpolate=True,
    ).run(verbose=True)

    ## pulls out the GL2 beads that belong to the upper and lower leaflet
    z_data_upper = (mc.results["z_surface"]["upper"]) / 10
    z_data_lower = (mc.results["z_surface"]["lower"]) / 10

    ## create a dictionary of
    leaflet, leaflet_residues = leaflet_finder.determine_leaflets(u, selection)

    return (X, Y, z_data_upper, z_data_lower, leaflet_residues)

#libraries for multiprocessing running this code directly
#from pathlib import Path
#from tqdm.contrib.concurrent import process_map

#Cailyns outline for code
#goal create a mesh from the PO4 beads in each leaflets
#determine the min dimension of the system. Since the system is NVT there should be no change in the box size
#determine the size of the bins
#both in x and y direction
#divide the min dimension by the bin size to get the number of bins 
#MembraneSpectralAnalysis 

#size of bins 1A or .1nm
# step = .5

# path = Path("/home/casakurai/scratch/packmem/")

def mesh(gro, selection, step ):
    import pickle
    import numpy as np
    import MDAnalysis as mda
    from membrane_spectral_analysis.base import MembraneSpectralAnalysis
    import leaflet_finder as leaflet_finder 
    
    #gro, selection, traj, path = args
    #u = mda.Universe(gro, str(traj) , continous=True)
    u = mda.Universe(gro, continous=True)

    print("units of the file", u.trajectory.units) #which is nm 
 

    # Average box dimensions over trajectory
    # might not need since I'm in the NVT ensemble
    avg_dim_x = np.mean([ts.dimensions[0] for ts in u.trajectory])
    avg_dim_y = np.mean([ts.dimensions[1] for ts in u.trajectory])
    print(avg_dim_y)


    # Determine number of bins
    num_bins_x = int(avg_dim_x / step)
    num_bins_y = int(avg_dim_y / step)

    # Create bin edges
    x_edges = np.linspace(0, avg_dim_x, num_bins_x + 1)
    y_edges = np.linspace(0, avg_dim_y, num_bins_y + 1)

    # Bin centers for plotting
    # I don't think this is being used in the MembraneSpectralAnalysis function
    x_centers = ((0.5 * (x_edges[:-1] + x_edges[1:])))/10
    y_centers = ((0.5 * (y_edges[:-1] + y_edges[1:])))/10
    X, Y = np.meshgrid(x_centers, y_centers, indexing="ij")

 

    # Checking z-range of lower glycerol bead
    ag = u.atoms.select_atoms("name GL2")
    z_coords = ag.positions[:, 2]            
    #("min z:", np.min(z_coords))
    #print("max z:", np.max(z_coords))
    #print("z spread:", np.ptp(z_coords))

    # Run MembraneSpectralAnalysis to get the quadrilateral mesh 
    # does it turn the units to angstroms internally?
    mc = MembraneSpectralAnalysis(
        u,
        select=selection,
        n_x_bins=num_bins_x,
        n_y_bins=num_bins_y,
        x_range=(0, avg_dim_x),
        y_range=(0, avg_dim_y),
        wrap=True,          # enable wrapping to avoid edge artifacts
        interpolate=True
    ).run(verbose=True)

    

    z_data_upper = (mc.results["z_surface"]["upper"])/10
    z_data_lower = (mc.results["z_surface"]["lower"])/10
    bin_area = mc.bin_areas

    # # Save only the results I care about. No need to save entire mc object
    # with open("membrane_curvature_smooth.pickle", "wb") as handle:
    #     pickle.dump({"X":X, "Y": Y, "Z_upper": z_data_upper ,"Z_lower": z_data_lower}, handle, protocol=pickle.HIGHEST_PROTOCOL)



    leaflet, leaflet_residues = leaflet_finder.determine_leaflets(u, selection)

    # # Save only the results I care about. No need to save entire mc object
    # with open("membrane_curvature_smooth_.5A.pickle", "wb") as handle:
    #     pickle.dump({"X":X, "Y": Y, "Z_upper": z_data_upper ,"Z_lower": z_data_lower, "bin_area": bin_area}, handle, protocol=pickle.HIGHEST_PROTOCOL)
    z_midplane = .5 * (z_data_upper + z_data_lower)

    
    return (X, Y, z_data_upper, z_data_lower,  leaflet, leaflet_residues, z_midplane)
# jobs = []
# systems = ["test"]
# lipid = "test"
# time = "test"


# for sys in systems:
#     #gro file
#     gro = path / "dopcdope_CGMartini.pdb"
#     #traj file
#     #selection
#     selection = f"name GL2"
#     inputs = (gro, selection, path)
#     jobs.append(inputs)

# process_map(mesh, jobs, max_workers=64)

#interpolate over the z_midplane z coordinates 
def midplane_interpolated(z_midplane):
    from scipy.interpolate import RegularGridInterpolator
    x = X[:, 0]
    y = Y[0, :]
    f_mid = RegularGridInterpolator((x,y, z_midplane))
    return f_mid 


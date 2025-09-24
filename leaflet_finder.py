
from MDAnalysis.analysis.leaflet import LeafletFinder
import MDAnalysis as mda

def largest_groups(atoms):
    """
    From a list of sizes, find out the indices of the two largest groups. These should correspond to the two leaflets of the bilayer.

    Keyword arguments:
    atoms -- list of sizes of clusters identified by LeafletFinder
    """
    largest = 0
    second_largest = 0

    for i in atoms:
        if atoms[i] > largest:
            largest_index = i
            largest = atoms[i]

    for i in atoms:
        if atoms[i] > second_largest and i != largest_index:
            second_largest_index = i
            second_largest = atoms[i]

    return (largest_index, second_largest_index)


def determine_leaflets(universe, selection = "all"):
    """
    From a selection of phosphates, determine which belong to the upper and lower leaflets.

    Keyword arguments:
    universe -- an MDAnalysis Universe object
    phosphateSelection -- a string specifying how to select the phosphate atoms (e.g. "name P1")
    """
    leaflets = {}

    # calculate the z value of the phosphates defining the bilayer (assumes bilayer is in x and y..)
    ag = universe.atoms.select_atoms(selection)
    bilayerCentre = ag.center_of_geometry()[2]

    # apply the MDAnalysis LeafletFinder graph-based method to determine the two largest groups which
    #  should correspond to the upper and lower leaflets
    phosphates = LeafletFinder(universe,selection)

    # find the two largest groups - required as the first two returned by LeafletFinder, whilst usually are the largest, this is not always so
    (a, b) = largest_groups(phosphates.sizes())

    # check to see where the first leaflet lies
    if phosphates.group(a).centroid()[2] > bilayerCentre:
        leaflets['upper'] = phosphates.group(a)
        leaflets['lower'] = phosphates.group(b)
    else:
        leaflets['lower'] = phosphates.group(a)
        leaflets['upper'] = phosphates.group(b)
    # --- Map residue IDs to leaflets ---
    #if I have issues with indexing the list might be due to how the resids are being saved as Numpy
    # residue_leaflet_map = {
    #     'upper': list(leaflets['upper'].resids),
    #     'lower': list(leaflets['lower'].resids)
    #}
    ##Cailyn took out because this has the format for each residue number it states upper vs lower, rather than having a list upper vs lower with all residues belonging to the appropriate leaflet
    # residue_leaflet_map ={}
    # for resid in leaflets['upper'].resids:
    #     residue_leaflet_map[resid] = 'upper'
    # for resid in leaflets['lower'].resids:
    #     residue_leaflet_map[resid] = 'lower'
   
    ##Cailyn added
    #make it nicer so two list upper vs lower
    leaflet_residues = {
        'upper': [int(r) for r in leaflets['upper'].resids],
        'lower': [int(r) for r in leaflets['lower'].resids]
    }

    # Return both AtomGroups and the residue map
    return leaflets, leaflet_residues

# #Cailyn: troubleshooting to determine if the residue is being named to the correct upper or lower leaflet 
# from pathlib import Path
# path = Path("/home/casakurai/scratch/packmem/")
# gro = path / "dopcdope_CGMartini.pdb"
# u = mda.Universe(gro, continous=True)
# determine_leaflets(u, "name PO4")

#Cailyn: I made to deal with buckled bilayers, but I didn't need it 
def determine_leaflets_local(universe, selection="all", step=0.5, sigma=2):
    from scipy.ndimage import gaussian_filter
    import numpy as np
    import MDAnalysis as mda

    """
    Determine leaflet membership for a buckled bilayer using local midplane smoothing.
    step: bin size in nm
    sigma: Gaussian smoothing standard deviation in bins
    """
    # Convert to Å for positions (MDAnalysis uses Å)
    step_angstrom = step * 10.0

    ag = universe.atoms.select_atoms(selection)
    positions = ag.positions.copy()

    # Average box size over trajectory
    avg_dim_x = np.mean([ts.dimensions[0] for ts in universe.trajectory])
    avg_dim_y = np.mean([ts.dimensions[1] for ts in universe.trajectory])

    num_bins_x = int(avg_dim_x / step_angstrom)
    num_bins_y = int(avg_dim_y / step_angstrom)

    x_edges = np.linspace(0, avg_dim_x, num_bins_x + 1)
    y_edges = np.linspace(0, avg_dim_y, num_bins_y + 1)

    x_centers = 0.5 * (x_edges[:-1] + x_edges[1:])
    y_centers = 0.5 * (y_edges[:-1] + y_edges[1:])
    X, Y = np.meshgrid(x_centers, y_centers)
    
    


    # Assign atoms to bins
    x_bins = np.digitize(positions[:, 0], x_edges) - 1
    y_bins = np.digitize(positions[:, 1], y_edges) - 1

    bin_atoms = [[[] for _ in range(num_bins_y)] for _ in range(num_bins_x)]

    for idx, (xb, yb) in enumerate(zip(x_bins, y_bins)):
        if 0 <= xb < num_bins_x and 0 <= yb < num_bins_y:
            bin_atoms[xb][yb].append(positions[idx, 2])

    # Calculate mean z in each bin
    local_midplane = np.full((num_bins_x, num_bins_y), np.nan)
    for i in range(num_bins_x):
        for j in range(num_bins_y):
            if bin_atoms[i][j]:
                local_midplane[i, j] = np.mean(bin_atoms[i][j])

    # Fill NaNs by interpolation from nearest non-NaN
    mask = np.isnan(local_midplane)
    if np.any(mask):
        from scipy.interpolate import NearestNDInterpolator
        coords = np.array(np.nonzero(~mask)).T
        values = local_midplane[~mask]
        interp = NearestNDInterpolator(coords, values)
        local_midplane[mask] = interp(np.array(np.nonzero(mask)).T)

    # Smooth the midplane
    local_midplane = gaussian_filter(local_midplane, sigma=sigma)

    # Classify atoms into leaflets
    leaflets = {"upper": [], "lower": []}
    for idx, (xb, yb) in enumerate(zip(x_bins, y_bins)):
        if 0 <= xb < num_bins_x and 0 <= yb < num_bins_y:
            mid_z = local_midplane[xb, yb]
            if positions[idx, 2] > mid_z:
                leaflets["upper"].append(idx)
            else:
                leaflets["lower"].append(idx)

    leaflets["upper"] = universe.atoms[leaflets["upper"]]
    leaflets["lower"] = universe.atoms[leaflets["lower"]]
    return leaflets

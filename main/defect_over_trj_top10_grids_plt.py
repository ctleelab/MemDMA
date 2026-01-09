#read in .txt files
#grid number: coords of triangle centroid, number of atoms mapped
from pathlib import Path
import ast  # safer than eval
import pandas as pd
import numpy as np

system = "DPPC_DPPS_strain.2"
data_path = Path(f"/home/casakurai/scratch/packmem/simulation_data_centered_final/analysis_{system}/extracted_gro_frames")

txt_files = sorted(data_path.glob("*atom_mapping.txt"))
list_defect_percent_area = []

for file in txt_files:
    
    data = pd.read_csv(file, sep="\t")  # header row is read automatically
    print(data.head())

    Grid_Num = ['Grid_Number']
    X = data['X']
    Y = data['Y']
    Z = data['Z']
    atom_mapping_value = data['Atom_Mapping_Value']

    n_top = int(0.10 * Z.size)
    z_max = np.max(Z)

    # rank-based mask
    top_idx = np.argsort(z_max - Z)[:n_top]
    mask = np.zeros_like(Z, dtype=bool)
    mask[top_idx] = True

    #total # of grids with height in the top 10% 
    total_grid_90 = sum(mask)
    print("total_grid_90", total_grid_90)

    #pull out the number of those that have defects <1 
    defect_grid_90 = sum((mask) & (atom_mapping_value < 1))
    print("defect_grid_90", defect_grid_90)


    #percent defect 
    defect_percent_area = (defect_grid_90/total_grid_90)*100

    list_defect_percent_area.append(defect_percent_area)
    print("file_name", file)
    print("defect_perecent_area", defect_percent_area)


# Save results
output_file = f"{data_path}/defect_trj_{system}_top10grids.txt"
with open(output_file, "w") as f:
    for val in list_defect_percent_area:
        f.write(f"{val}\n")

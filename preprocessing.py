#preprocessing 
#create a list that contains {atom number, atom_name, leaflet assignment, coords of atom }
#only create this list once to save on computing time

import re 
import numpy as np
import grid_slicing as gs
import util

#create a dictionary with atom: (atom_name), (leaflet_assignment), (atom_coords)
atoms = {}
def atom_features(lines_gro, list_types_atoms, leaflet_residue_new):
    for i, line_atom in enumerate(lines_gro):            
        if i > 1 and i < (len(lines_gro)-1):
            # Extract residue name and number
            res_name_number = line_atom[0:12].strip()
            match = re.match(r"(\d+)([A-Za-z]+)", res_name_number)
            if match:
                res_number = int(match.group(1))
            
            # Extract atom name and number
            atm_name_number = str(line_atom[12:22].strip())
            atm_name = None
            atm_number = None
            for j in list_types_atoms:
                if atm_name_number.startswith(j):
                    atm_name = j
                    atm_number = atm_name_number[len(j):]
                    break
            
            if atm_name is None or atm_number is None:
                continue  # Skip if atom name or number is not found

            # Extract coordinates of the atom
            x_coord = float(line_atom[20:28].strip())
            y_coord = float(line_atom[28:36].strip())
            z_coord = float(line_atom[36:44].strip())
            atm_coords = np.array([x_coord, y_coord, z_coord])

            # Determine leaflet assignment
            leaflet = gs.find_leaflet(res_number, leaflet_residue_new)

            # Add atom details to the dictionary
            atoms[atm_number] = (atm_name, leaflet, atm_coords)
    return atoms

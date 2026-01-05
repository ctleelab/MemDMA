import re
import numpy as np
import grid_slicing as gs
import util as util


atoms = {}


def atom_features(lines_gro: list, list_types_atoms: list, leaflet_residue_new):
    """
    Create a dictionary that describes all features of an atom.
    Args:
        ##input: grolines, list of the types of atoms in the gro file, which residues are in the upper vs lower leaflet
        lines_gro (list): lines of the gro file
        list_types_atoms (list): list of the types of atoms in the gro file
        leaflet_residue_new (dict): {"upper": [list of residues in upper leaflet],"lower": [list of residues in lower leaflet]}

    Returns:
        atoms (dict): dictionary of atom features  {atm#: 'GL2', 'lower', array([18.222,  4.37 , 14.85 ]), 'DOPC GL2', resnumber }

    """
    for i, line_atom in enumerate(lines_gro):
        if i > 1 and i < (len(lines_gro) - 1):
            # Extract residue name and number
            res_name_number = line_atom[0:12].strip()
            match = re.match(r"(\d+)([A-Za-z]+)", res_name_number)
            if match:
                res_number = int(match.group(1))
                res_name = str(match.group(2))

            # Extract atom name and number
            atm_name_number = str(line_atom[12:22].strip())
            atm_name = None
            atm_number = None
            for j in list_types_atoms:
                if atm_name_number.startswith(j):
                    atm_name = j
                    atm_number = (atm_name_number[len(j) :]).strip()
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

            # lipid and atom name
            lipid_atm_name = str(res_name + " " + atm_name)
            # Add atom details to the dictionary
            atoms[atm_number] = (
                atm_name,
                leaflet,
                atm_coords,
                lipid_atm_name,
                res_number,
            )
    return atoms

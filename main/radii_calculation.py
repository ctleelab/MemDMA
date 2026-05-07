#!/usr/bin/env python3
import sys

# pixi inputs
lipid_name = sys.argv[1]
atom_type = sys.argv[2]
sigma_value = float(sys.argv[3])
print(sigma_value)
polarity = sys.argv[4]


def rmin_from_sigma(lipid_name, atom_type, sigma_value, polarity):
    true_sigma_value = sigma_value/2 #divide by two b/c sigma value is for two beads interacting, but we want the radius of one bead
    radii = round((2 ** (1 / 6) * true_sigma_value),4)
    lipid_name = lipid_name.upper()
    atom_type= atom_type.upper()
    print(lipid_name, atom_type, radii, polarity)

rmin_from_sigma(lipid_name, atom_type, sigma_value, polarity)

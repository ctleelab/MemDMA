from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree

list_types_atoms = [
    "W",
    "NA",
    "CL",
    "NC3",
    "PO4",
    "GL1",
    "GL2",
    "C1A",
    "D2A",
    "C3A",
    "C4A",
    "C1B",
    "D2B",
    "C3B",
    "C4B",
    "C2A",
    "C2B",
]
list_types_atoms_lipid = [
    "NC3",
    "PO4",
    "GL1",
    "GL2",
    "C1A",
    "D2A",
    "C3A",
    "C4A",
    "C1B",
    "D2B",
    "C3B",
    "C4B",
    "C2A",
    "C2B",
]
SIZE = 0.5
lipid = "name NC3 PO4 GL1 GL2 C1A D2A C3A C4A C1B D2B C3B C4B C2A C2B"
membrane_sel = "resname DPPC DOPC"


trajectory_path = Path("/Users/ctlee/PackMem-leelab/tests")
analysis_path = Path("/Users/ctlee/PackMem-leelab/tests")


def parse_radius_aliphatic_dict(filename: str) -> dict:
    property_dict = {}
    for line in open(filename, "r"):
        line = line.split("#")[0].strip()  # Remove comments and whitespace
        if line:
            _tmp = line.split()
            key = f"{_tmp[0]} {_tmp[1]}"
            property_dict[key] = (float(_tmp[2]), _tmp[3])
    return property_dict


def cross(v1, v2):
    return np.array(
        [
            v1[1] * v2[2] - v1[2] * v2[1],
            v1[2] * v2[0] - v1[0] * v2[2],
            v1[0] * v2[1] - v1[1] * v2[0],
        ]
    )


def norm(v):
    return np.sqrt(np.dot(v, v))

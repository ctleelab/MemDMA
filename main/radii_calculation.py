#!/usr/bin/env python3
import sys
import os


def rmin_from_sigma(sigma_value):
    # sigma is for two beads interacting; divide by 2 for the radius of one bead
    true_sigma_value = sigma_value / 2
    return round((2 ** (1 / 6) * true_sigma_value), 4)


def parse_ff_file(ff_file):
    beads = []
    with open(ff_file) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith(";"):
                continue
            parts = line.split(";")[0].split()
            # columns: i type resnr residue atom cgnr charge
            if len(parts) >= 5:
                bead_type = parts[1]
                residue = parts[3]
                atom = parts[4]
                beads.append((bead_type, residue, atom))
    return beads


def parse_self_sigmas(itp_file):
    sigma_dict = {}
    in_nonbond = False
    with open(itp_file) as f:
        for line in f:
            stripped = line.strip()
            if "[ nonbond_params ]" in stripped:
                in_nonbond = True
                continue
            if in_nonbond:
                if stripped.startswith("["):
                    break
                if not stripped or stripped.startswith(";"):
                    continue
                parts = stripped.split()
                if len(parts) >= 5 and parts[0] == parts[1]:
                    sigma_dict[parts[0]] = float(parts[3])
    return sigma_dict


if __name__ == "__main__":
    script_dir = os.path.dirname(os.path.abspath(__file__))
    if len(sys.argv) < 2:
        print(f"Usage: {sys.argv[0]} <ff_file> [itp_file]", file=sys.stderr)
        sys.exit(1)
    ff_file = sys.argv[1]
    itp_file = (
        sys.argv[2]
        if len(sys.argv) > 2
        else os.path.join(script_dir, "..", "toppar", "martini_v3.0_openbeta", "martini_v3.0.0.itp")
    )

    beads = parse_ff_file(ff_file)
    sigma_dict = parse_self_sigmas(itp_file)

    print("# RESNAME ATOMTYPE RADIUS_ANGSTROM ALIPHATIC")
    for bead_type, residue, atom in beads:
        if bead_type not in sigma_dict:
            print(f"WARNING: no self-interaction sigma for bead type '{bead_type}'", file=sys.stderr)
            continue
        sigma = sigma_dict[bead_type]
        radii = rmin_from_sigma(sigma) * 10  # nm to angstrom 
        aliphatic = "a" if bead_type.startswith("C") else "n" #determines if the bead is polar or not
        print(f"{residue}  {atom}   {radii} \t{aliphatic}")

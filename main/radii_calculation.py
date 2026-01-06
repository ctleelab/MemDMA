def rmin_from_sigma(sigma):
    return 2 ** (1 / 6) * sigma


# input the sigma value from the martini_v3.0.0.itp
# did the interactions between the same two beads so that you can easily divide by 2 to represent the radii of one bead
# Q1    Q1  1 4.700000e-01    3.980000e+00
# NC3
print("M3 Q1", rmin_from_sigma(3.98) / 2)

# Q5    Q5  1 4.700000e-01    6.450000e+00
# PO4
print("M3 Q5", rmin_from_sigma(6.45) / 2)

# SN4a  SN4a  1 4.100000e-01    2.600000e+00
# GL1 and GL2
print("M3 SN4a", rmin_from_sigma(2.6) / 2)

# C1    C1  1 4.700000e-01    3.390000e+00
# C2A, C3A, C4A, C1A
print("M3 C1", rmin_from_sigma(3.39) / 2)

# C4h   C4h  1 4.700000e-01    3.520000e+00
# D2A, D2B
print("M3 C4h", rmin_from_sigma(3.52) / 2)

# SC1   SC1  1 4.100000e-01    2.350000e+00
# C1B
print("M3 SC1", rmin_from_sigma(2.35) / 2)

def rmin_from_sigma(sigma):
    return 2 ** (1 / 6) * sigma


# input the sigma value from the martini_v3.0.0.itp
# did the interactions between the same two beads so that you can easily
# divide by 2 to represent the radii of one bead
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
# C2A, C3A, C4A, (C1Aand C1B for DO lipids)
print("M3 C1", rmin_from_sigma(3.39) / 2)

# C4h   C4h  1 4.700000e-01    3.520000e+00
# D2A, D2B
print("M3 C4h", rmin_from_sigma(3.52) / 2)

# SC1   SC1  1 4.100000e-01    2.350000e+00
# C1B and C1A for DP lipids
print("M3 SC1", rmin_from_sigma(2.35) / 2)


#P6    P6  1 4.700000e-01    4.990000e+00
#CNO
print("M3 P6", rmin_from_sigma(4.99) / 2)





#%% 
import numpy as np


#%%
# MARTINI 2
# ; O 	- supra attractive: 	(eps=5.6, s=0.47) 
# ; I 	- attractive: 		(eps=5.0, s=0.47) 
# ; II 	- almost attractive: 	(eps=4.5, s=0.47)
# ; III 	- semi attractive:	(eps=4.0, s=0.47) 
# ; IV 	- intermediate: 	(eps=3.5, s=0.47) 
# ; V 	- almost intermediate: 	(eps=3.1, s=0.47) 
# ; VI 	- semi repulsive: 	(eps=2.7, s=0.47) 
# ; VII 	- almost repulsive: 	(eps=2.3, s=0.47) 
# ; VIII	- repulsive: 		(eps=2.0, s=0.47)
print('M2', 2**(1/6)*4.7/2)
# ; IX 	- super repulsive: 	(eps=2.0, s=0.62)
print('M2 repulse', 2**(1/6)*6.2/2)


def rmin_from_sigma(sigma):
    return 2**(1/6) * sigma



#%%
# MARTINI 3
# ; id 	type 	resnr 	residu 	atom 	cgnr 	charge    (mass)
#    1     Q1   1  DDPC  NC3   1 	1.0
#    2     Q5   1  DDPC  PO4   2 	-1.0
#    3   SN4a   1  DDPC  GL1   3 	0
#    4   SN4a   1  DDPC  GL2   4 	0
#    5    C4h   1  DDPC  D1A   5 	0
#    6    C5h   1  DDPC  D2A   6 	0
#    7    C5h   1  DDPC  D3A   7 	0
#    8    C4h   1  DDPC  D4A   8 	0
#    9    C4h   1  DDPC  D5A   9 	0
#   10    C4h   1  DDPC  D1B  10 	0
#   11    C5h   1  DDPC  D2B  11 	0
#   12    C5h   1  DDPC  D3B  12 	0
#   13    C4h   1  DDPC  D4B  13 	0
#   14    C4h   1  DDPC  D5B  14 	0

#Q1    Q1  1 4.700000e-01    3.980000e+00
print("M3 Q1", rmin_from_sigma(3.98)/2)
#Q5    Q5  1 4.700000e-01    6.450000e+00
print("M3 Q5", rmin_from_sigma(6.45)/2)
#SN4a  SN4a  1 4.100000e-01    2.600000e+00
print("M3 SN4a", rmin_from_sigma(2.6)/2)
#    C1    C1  1 4.700000e-01    3.390000e+00
print("M3 C1", rmin_from_sigma(3.39)/2)
#   C4h   C4h  1 4.700000e-01    3.520000e+00
print("M3 C4h", rmin_from_sigma(3.52)/2)
#   SC1   SC1  1 4.100000e-01    2.350000e+00
print("M3 SC1", rmin_from_sigma(2.35)/2)

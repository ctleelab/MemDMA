## About this repo
- Currently only suitable for Martini systems
- Can analyze buckled and flat membranes

## Repository structure
```bash
│── main/                           
│
│── input/                          #where your add your trajectory files that will be analyzed               
│
│── test/                           # Folder where all test cases and their outputs live 
│   ├──output/                      # output of all test analysis
│   ├── .gro                        # .gro file for test cases buckled and flat membranes 
│   ├── .xtc                        # .xtc file for test cases buckled and flat membranes 
│
│── analysis/                       # outputs from analysis scripts
│
│── misc_scripts/   
│
│── README.md                       # Project documentation
│── .gitignore                      # Files ignored by Git
 ```



## Requirements 
- All packages should be present in the pixi environment

## Set-up pixi 
Install pixi:
```bash
curl -fsSL https://pixi.sh/install.sh | sh
 ```

Enter to pixi environment
```bash
pixi shell
 ```

Exit pixi environment
```bash
exit
 ```

Test cases
```bash
pixi r flat-test
pixi r strain2-test
 ```

## Setup for your lipids of interest
Currently, this can only analyze membranes composed of the following lipids:
- DOPC, DPPC, DOPA, DPPA, DOPS, DPPS

If your lipids are not in the list they must be added. To do so, complete the following for each unique bead type in the lipid. 
1) Calculate the radii of each of the bead by completing the following. 
- In the pixi.toml file change the "lipid_name", "atom_type", "sigma_value", and "polarity"  for your choosen bead. The polarity value is n if the bead is polar and a if the bead is not polar.
- The sigma value is determined between two of the same bead types. And can be found in the martini_v3.0.0.itp file. 
  ```bash
    Q1    Q1  1 4.700000e-01    3.980000e+00
  ```
2) Then run the following: 
  ``` bash  
  pixi r radii-calculation 
  ```
3) Copy the ouputted text which will represent the "RESNAME ATOMTYPE RADIUS_ANGSTROM ALIPHATIC" and copy directly into the "vdw_radii_CGMartini3_openbeta2_updated.txt" file.


## Running your own analysis

# Curvature analysis 
In the pixi.toml file add the following 
```bash 
curvature-analysis-filename = "main/curvature.py input/your_file_name.gro tests/your_file_name.xtc 5 your_file_name 40 output/curvature_analysis/your_file_name"
```
Then execute the following in the terminal 
```bash 
pixi r curvature-analysis-filename
```

# Defect analysis
In the pixi.toml file add the following 
```bash 
defect-analysis-filename = "main/packmem.py -i tests/DOPC_buckled.gro -t tests/DOPC_buckled.xtc -r main/vdw_radii_CGMartini3_openbeta2_updated.txt -d 6.3801 -o tests/output/test_buckled"
```
Then execute the following in the terminal 
```bash 
pixi r defect-analysis-filename
```
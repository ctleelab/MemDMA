# Membrane Defect Mapping and Analysis (MemDMA)
- Geometry-preserving membrane defect mapping for flat and curved membranes
- Currently only suitable for Martini systems

#  Repository structure
```bash
│── main/                           # where main MemDMA code lives                    
│   ├──trajectory_prep/             # preparing the simulation files for mapping and analysis
│  
│── input/                          #where your add your trajectory files (.gro and .xtc) that will be analyzed               
│
│── tests/                          # Folder where all test cases and their outputs live 
│   ├──output/                      # output of all test analysis
│   ├── .gro                        # .gro file for test cases of buckled and flat membranes 
│   ├── .xtc                        # .xtc file for test cases of buckled and flat membranes 
│
│── output/                         # General outputs from executions
│   ├──curvature_analysis/          # Curvature analysis outputs
│   ├──defect_mapping/              # Defect analysis outputs
│   ├──defect_analysis/  
│
│── toppar/                         # Martini3 topology files
│ 
│── README.md                       # Project documentation
│── .gitignore                      # Files ignored by Git
 ```



#  Requirements 
- All packages should be present in the pixi environment

# Set-up 

## Set-up pixi 
Install pixi:
```bash
curl -fsSL https://pixi.sh/install.sh | sh
 ```

Enter to pixi environment
```bash
pixi shell
 ```

Test cases
```bash
pixi r flat-test
pixi r strain2-test
 ```

Exit pixi environment
```bash
exit
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
3) Copy the ouputted text which will represent the "RESNAME ATOMTYPE RADIUS_ANGSTROM ALIPHATIC" and place directly into the "vdw_radii_CGMartini3_openbeta2_updated.txt" file.



## Formating .gro and .xtc files prior to analysis
If the your system contains lipids that are not in the list above you will have to add them. 
In the util.py file add your lipids to the string "lipids".
1) Strip the membrane of ions and waters.
   - In the command,  "strip-lipid-filename" in the pixi.toml file change each instance of your_filename with the name of your file of interest.
   - Then run the following:
   ``` bash 
   pixi r strip-lipid-filename
   ```
2) After stripping the membrane, center the buckled membranes. This makes visualization and analysis of defects easier.
   - In the command,  "center-membrane-filename" in the pixi.toml file change each instance of your_filename with the name of your file of interest.
   - Then run the following:
   ``` bash 
   pixi r center-membrane-filename
   ```
The stripped_centered or strippedx files should be used for all further analysis.



# Running your own analysis

## Defect Mapping
This determines the defects present in the membrane throughout the simulation.

In the pixi.toml replace the "your_file_name" for the file name of your .gro and .xtc

```bash 
defect-mapping-filename = "python main/packmem.py -i input/your_file_name.gro -t input/your_file_name.xtc -o output/defect_mapping/your_file_name"
```
Then execute the following in the terminal 
```bash 
pixi r defect-mapping-filename
```

Output
1) A .png file of the defect map for the last frame of the trajectory.
2) A .npy array of defects, centroids, vertices for the upper and lower leaflet
3) A .npy array of shape, Xshape, and Yshape

## Defect Analysis
Takes the data from the defect mapping and creates visualizations for the lower and upper leaflet independently. 
1) A plot of height vs defect probability for all defect types (i.e shallow, deep, all)
2) A plot of defect size probability.
3) A plot of defect coverage across the simulation.
4) .csv file for the average defect coverage and average defect size. 

```bash
defect-analysis-filename = "python main/aggregate_map_defects_plt.py -input output/defect_mapping/your_file_name -output output/defect_analysis/your_file_name -l your_file_name"
```

Then execute the following in the terminal 
```bash 
pixi r defect-analysis-filename
```

## Curvature analysis 
This determines the curvature of the buckled membrane throughout the simulations by fitting a circle to the peak of the membrane.
  
In the pixi.toml replace the "your_file_name" for the file name of your .gro and .xtc
```bash 
curvature-analysis-filename = "python curvature-analysis-filename = "main/curvature.py -gro input/your_file_name.gro -xtc input/your_file_name.xtc -s your_file_name -out output/curvature_analysis/your_file_name"
```
Then execute the following in the terminal 
```bash 
pixi r curvature-analysis-filename
```

Output 
1) A .png file of the circle of best fit for each frame of the trajectory.
2) A .npy file that describes the radius, curvature of the perfect fit circle for each frame.
3) A .csv file that contains a summary of the curvature values.




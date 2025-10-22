Folders
1) Main: contains all relevant scripts
2) tests: contains test .gro files (both a buckled and flat membrane of 100% DOPC)


Running code: 
1) go to main folder
2) run the following:
    ```python base_clean_v2.py -i ../tests/testing_flat_lipid_only.gro  -p param_CGMartini3_openbeta2.txt -r vdw_radii_CGMartini3_openbeta2_inaccuratedata_testing.txt -o output -d 1.0 -t deep```
3) ouput:

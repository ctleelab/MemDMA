#-*- coding: utf-8 -*-
# Funtions for matrix
# R. Gautier A. Bacle 2015

import math
import sys
import PackMem_original.BasicFunctions as bfrg

# matrix size  (square size = 1A)
SIZE=1.0

# calculate limit size (half diagonal) 3D
SIZE_SIDE=0.5*(math.sqrt((math.pow(SIZE, 2))* 3.))

# get radius from [Lipid_name Atom_name]
def get_radius(radius, res_name, atom_name):
    key = res_name+" "+atom_name
    if key not in radius.keys():
        print("ERROR : Association of lipid %s and atom %s not found in parameter file"%(res_name, atom_name))
        sys.exit()
    return radius[key]

# get aliphatic flag [Lipid_name Atom_name]
def get_aliphatic(aliphatic, res_name, atom_name):
    key = res_name+" "+atom_name
    if key not in aliphatic.keys():
        print("ERROR : Association of lipid %s and atom %s not found in parameter file"%(res_name, atom_name))
        sys.exit()
    return aliphatic[key]

# fill the matrix for each atom
##Cailyn: Z-coordinate is not stored. Therefore the acyl chain beads at the top of the membrane are treated the same as the bottom of the leafet or in between
##Cailyn: It's only treated differently if there is no acyl chain bead at all, and the acyl chain bead might be in the other leaflet
##Cailyn: you don't get to understand the true depths of the defects
##Cailyn: I might try this code out, but if I interested in the depths of the defects I think this is not possible with this code
def fill_matrix(matrix,coordtmp,res_number,res_name,atom_name, listX,listY,listZ,
                radius_res, FlagPDtype, aliph_atoms):
    v=5
    ##Cailyn: local neighborhood set up, determining which grids should the atom even effect
    # for each lipid get listZ
    ##Cailyn:cordtmp for the atom
    xLip = coordtmp[0]
    yLip = coordtmp[1]
    zLip = coordtmp[2]

    # works only around X, Y position of the matrix +/- v cases
    ##Cailyn: finds where the atoms falls in the matrix 
    iX,iY = find_X_Y(coordtmp, listX, listY)

    #makes sure that we dont go outside of the grid boundaries
    iX = check_edges(iX, v, len(listX))
    iY = check_edges(iY, v, len(listY))
    
    #only considers a small neighborhood aroun the atom, rather than the entire membrane
    listXM = listX[iX - v : iX + (v + 1)]
    listYM = listY[iY - v : iY + (v + 1)]

    dist_lim = (SIZE+radius_res) ** 2
    dist_meet = (SIZE_SIDE + radius_res) ** 2

    #only looking at grid points that are within size+radius ^2 distance away in the x,y,z dimension
    ##Cailyn: how is it possible to hit the X,Y,Z but then not meet the dist_meet critera
    ##Cailyn: dist_lim will be longer than dis_meet, because dis_meet is the diagonal length
    ##Cailyn: the for loop does the initilization where you will replaced the NA in the matrix with 0 
    ##Cailyn: Then you will go through the if statements to update the matrix
    
    for indZ, sliceZ in enumerate(listZ):
        distZ_to_sliceZ = bfrg.dist_oneAxis(zLip,sliceZ)
        if distZ_to_sliceZ > dist_lim:
            continue
        for indX,sliceX in enumerate(listXM):
            distX_to_sliceX = bfrg.dist_oneAxis(xLip, sliceX)
            if distX_to_sliceX > dist_lim:
                continue
            for indY,sliceY in enumerate(listYM):
                distY_to_sliceY = bfrg.dist_oneAxis(yLip, sliceY)
                if distY_to_sliceY > dist_lim:
                    continue
                
                X = indX + (iX - v)
                Y = indY + (iY - v)
                Z = indZ
                coordCenter = [float(sliceX), float(sliceY), float(sliceZ)]
                distance = bfrg.dist(coordCenter, coordtmp)

                ##Cailyn:initialize the matrix as empty
                if matrix[X][Y]=="NA":
                    matrix[X][Y]=0.
                val_mat = matrix[X][Y]

                ##Cailyn:update the matrix
                if distance <= dist_meet:
                    # shallow PDefects (2) or all PDefects (0)
                    if FlagPDtype == 2 or 0:
                        #determines the contribution depending on atom type 
                        matrix[X][Y] = setDefects(atom_name, aliph_atoms, val_mat)
                    # deep PDefects (1)
                    # 1 is added when polar? 
                    else:
                        matrix[X][Y] +=1.
                        

    return matrix


# fill matrix cell depending on the defect type (Deep, Shallow)
def setDefects(atom_name, list_aliphatic, val_mat):
    if list_aliphatic == "a":
        val_mat += 0.001
    else:
        val_mat += 1
    return val_mat

def initialize_matrix3D(val1, val2, val3, defaut):
    Matrix=list(range(0,val1))
    for i in range(0, val1):
        Matrix[i]=list(range(0,val2))
        for j in range(0, val2):
            Matrix[i][j]=list(range(0,val3))
            for k in range(0, val3):
                Matrix[i][j][k]=defaut
    return Matrix

def initialize_matrix2D(val1, val2, defaut):
    Matrix=list(range(0,val1))
    for i in range(0, val1):
        Matrix[i]=list(range(0,val2))
        for j in range(0, val2):
            Matrix[i][j]=defaut
    return Matrix       

def binarize_matrix(matrix, matrix_ini, val=0.):
    for i in range(0, len(matrix)):
        for j in range(0, len(matrix[0])):
            #print(matrix[i][j])
            if matrix[i][j] == "NA":
                matrix_ini[i][j]=0.
            else:
                if matrix[i][j] > val:
                    matrix_ini[i][j]=1.
                else:
                    matrix_ini[i][j]=0.
    return matrix_ini

def binarize_matrix_whithout0(matrix, matrix_ini, val1=0, val2=0.99):
    for i in range(0, len(matrix)):
        for j in range(0, len(matrix[0])):
            #print(matrix[i][j])
            if matrix[i][j] == "NA":
                matrix_ini[i][j]=0.
            else:
                if (matrix[i][j] > val1 and  matrix[i][j] < val2):
                    matrix_ini[i][j]=0.
                else:
                    matrix_ini[i][j]=1.
    return matrix_ini

def binarize_matrix_whithout0_bis(matrix, matrix_ini, val=0.99):
    for i in range(0, len(matrix)):
        for j in range(0, len(matrix[0])):
            #print(matrix[i][j])
            if matrix[i][j] == "NA":
                matrix_ini[i][j]=0.
            else:
                if (matrix[i][j] > 0. and  matrix[i][j] < val):
                    matrix_ini[i][j]=0.
                else:
                    matrix_ini[i][j]=1.
    return matrix_ini

##Cailyn: strange way to deal with bins that still have NA after labeling
##Cailyn: I believe the clust_edge is just a category of bins that we will ignore. 
def clean_NA_inside(Matrix_labels, clust_edge, Matrix_ini, total_edge):
    ## Cailyn: dictionary matrix_label:count of label 
    ## Cailyn: All NAs mapping to the clust_edge 
    clustPb={}
    for i in range(0, len(Matrix_ini)): ##Caiyn: row length 
        for j in range(0, len(Matrix_ini[0])): ##Cailyn: column length
            if Matrix_ini[i][j] == "NA" and Matrix_labels[i][j] not in clust_edge:
                if Matrix_labels[i][j] in clustPb:
                    clustPb[Matrix_labels[i][j]]+=1
                else:
                    clustPb[Matrix_labels[i][j]]=1
                Matrix_labels[i][j]=clust_edge[0] #Cailyn: stores the Matrix[i][j] that has an NA value to the cluster_edge, and increases the count by one
                total_edge+=1
    return Matrix_labels, total_edge, clustPb
            
def determinelastNA(matrix):
    for i in range(len(matrix)-1, -1, -1):
        print(i, matrix[i][0][0])

# function that determine from x,y coordinates the central value in matrix
def find_X_Y(coord, matX, matY):
    tmpX = int(coord[0])
    tmpY = int(coord[1])
    iX = matX.index(tmpX)
    iY = matY.index(tmpY)
    #print(tmpX, matX, iX)
    return iX,iY

# compute the difference between the last value in matrixZ and the Z for one atom
def find_Z(coordZ, matZ):
    return matZ[-1]-coordZ

def check_edges(val, val_lim1, val_lim2):
    if val < val_lim1:
        return val_lim1
    elif val >= val_lim2-val_lim1:
        return val_lim2-(val_lim1+1)
    else:
        return val

def find_indexfromvalue(value, liste):
    mini=100.
    mini_i=len(liste)
    for i, val in enumerate(liste):
        if math.fabs(val-value) < mini:
            mini=math.fabs(val-value)
            mini_i=i
    return mini_i
    
def modify_matrix(label_matrix, binarized_matrix, label_value):
    for i in range(0, len(label_matrix)):
        for j in range(0, len(label_matrix[0])):
            if label_matrix[i][j] in label_value:
                binarized_matrix[i][j]= 0.0
    return binarized_matrix
    

# -*- coding: utf-8 -*-
# R. Gautier 2014
# Basic functions for packing defects determination

import numpy as np


# load file return list[line(string)]
def read_file(filename):
    try:
        with open(filename) as f:
            data = f.readlines()
    except:
        print("ERROR : Something went wrong with the file %s" % filename)
    return data


def read_radius(filename):
    dico = {}
    data = read_file(filename)
    for i in range(0, len(data)):
        data[i] = data[i].strip()  ##remove any leading/trailing whitespace
        if data[i] != "":  ##skip empty lines
            data[i] = data[i].split()  ##split the line based upon the whitespaces
            cle = data[i][0] + " " + data[i][1]  ##combine lipid and atom name
            dico[cle] = (
                float(data[i][2]) / 10
            )  ##then store the radius value to the lipid atom name {(DOP_D2A: raidus)} #divide by 10 because the radius previously presented 10A was in Angstroms, but now 10 in the context represent 10nm
    return dico


# new 2016-08-18
# open aliphatic table-> dict["RES(3L) ATOMname"]= p (polar) a (aliphatic)
def read_aliphatic(filename):
    dico = {}
    data = read_file(filename)
    for i in range(0, len(data)):
        data[i] = data[i].strip()
        if data[i] != "":
            data[i] = data[i].split()
            cle = data[i][0] + " " + data[i][1]
            dico[cle] = data[i][3]
    return dico


# Euclidian distance
def Edist3D(coord1, coord2):
    magnitude = np.linalg.norm(coord1 - coord2, axis=-1)
    return magnitude

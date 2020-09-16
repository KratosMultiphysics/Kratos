import numpy
import os
import sys
import math
import re, glob, subprocess, time, os, warnings
###############
# Definitions
###############

#cwd = os.getcwd()
#sys.path.append(cwd)
#working_path = cwd
#n_step = 150
#para_path_mod = working_path + '/airfoil_Structured.cntl'
#start_step = 0
#primary_grid_filename = working_path + '/Mesh/NeueGeo_WKA2D_MembraneSeparate_scaliert.grid'
#tau_path = '/work/piquee/Softwares/TAU/taudir_repos.2019.13.08/bin'
#rotate = 0

# Check if the path exists
def CheckIfPathExists(path):
    if not os.path.exists(path):
        raise Exception('Path: "{}" not found'.format(path))

# Write new file from the 'surface.pval*'
def WriteNewFile(filename_original, filename_write, startLine, stopLine):
    f = open(filename_original,'r')
    list_lines_original = f.readlines()
    list_lines = filter(None,list_lines_original)
    with open(filename_write,'a') as fwrite:
        zone_counter = 0
        for elem in list_lines[startLine:stopLine]:
            if "ZONE" in elem:
                zone_counter += 1
                if zone_counter > 1:
                    break
                else:
                    fwrite.write(elem)
            elif elem != "\n":
                fwrite.write(elem)
            
    f.close()
                
# Count lines first zone between 'ZONE' and end Connectivity for the same zone
def Countline(filename_original, word):
    f = open(filename_original,'r')
    list_lines = f.readlines()
    for elem in list_lines:
        if word in elem:
            index = list_lines.index(elem)
    f.close()
    return index

# Gives the length of the file 
def Endline(filename_original):
    f = open(filename_original,'r')
    list_lines = f.readlines()
    endline = len(list_lines)
    return endline

#number1 = Countline(cwd + "/Outputs/airfoilSol.surface.pval.unsteady_i=150_t=7.50000e-01.dat","MEMBRANE_UP")
#print(number1)
#number2 = Countline(cwd + "/Outputs/airfoilSol.surface.pval.unsteady_i=150_t=7.50000e-01.dat","MEMBRANE_DOWN")
#print(number2)

#with open(cwd + "/Outputs/airfoilSol.surface.pval.unsteady_i=150_t=7.50000e-01.dat",'r') as fread:
#    with open(cwd + '/test.dat','w') as fwrite:
#        header1 = fread.readline()
#        fwrite.write(header1)
#        header2 = fread.readline()
#        fwrite.write(header2)

#WriteNewFile(cwd + "/Outputs/airfoilSol.surface.pval.unsteady_i=150_t=7.50000e-01.dat",cwd + '/test.dat',number1, number2)


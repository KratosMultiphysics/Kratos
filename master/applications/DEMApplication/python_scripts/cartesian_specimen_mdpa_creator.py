# To use this script just type in the terminal: python (or python3) cartesian_specimen_mdpa_creator.py name_of_your_case
# It creates a file called name_of_your_caseDEM.mdpa.
# Remember to overwrite the Properties and the boundary conditions if necessary!!

from KratosMultiphysics import *
import KratosMultiphysics.DEMApplication as DEMapp
import sys

if len(sys.argv) < 2:
    raise Exception("You must specify the name of the file to be printed. The string 'DEM.mdpa' will be appended to the provided string. ")

DEMapp.PreUtilities().CreateCartesianSpecimenMdpa(sys.argv[1])

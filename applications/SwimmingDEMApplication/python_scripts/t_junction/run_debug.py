import os
import sys
sys.path.append(os.getcwd())
import json
import fileinput
from KratosMultiphysics import *
import KratosSwimmingDEM as script
import swimming_DEM_algorithm

for i in range(10):
    test = script.Solution()
    test.Run()

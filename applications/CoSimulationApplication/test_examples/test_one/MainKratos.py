from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import json
import KratosMultiphysics
from KratosMultiphysics import *
import KratosMultiphysics.CoSimulationApplication as CoSimulationApplication
from __init__ import *
from co_simulation_analysis import CoSimulationAnalysis

filename = '/home/aditya/work/Kratos/applications/CoSimulationApplication/test_examples/test_one/test_json.json'

#Read JSON data into the datastore variable
if filename:
    with open(filename, 'r') as f:
        json_data = json.load(f)


cs = CoSimulationAnalysis(json_data)
cs.PrintInfo()
cs.Run()
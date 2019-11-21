from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

from KratosMultiphysics.CoSimulationApplication.co_simulation_tools import ImportDataStructure
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools

import numpy as np

from sys import argv

import numpy as np
from copy import deepcopy
import multiprocessing
import os
import subprocess
from os.path import join
import time

def wait_message(dir, message):
    file = join(dir, message + ".msg")
    while not os.path.isfile(file):
        time.sleep(0.01)
    os.remove(file)
    return

# Check number of command line arguments
if len(argv) != 2:
    err_msg = 'Wrong number of input arguments!\n'
    err_msg += 'Use this script in the following way:\n'
    err_msg += '    "python Test.py <parameter-file>.json"\n'
    raise Exception(err_msg)


# Import data structure
parameter_file_name = argv[1]
cs_data_structure = ImportDataStructure(parameter_file_name)

# Import parameters using the data structure
with open(parameter_file_name, 'r') as parameter_file:
    parameters = cs_data_structure.Parameters(parameter_file.read())

settings = parameters['solver_wrappers'][0]['settings']

cmd = 'abaqus cae noGUI=script.py'

dir_csm = settings["working_directory"].GetString()
print(dir_csm)

abaqus_process = subprocess.Popen(cmd, executable='/bin/bash',
                                               shell=True, cwd=dir_csm)

wait_message(dir_csm, "Abaqus_Opened")

abaqus_process.terminate()

print(" \n Detected Abaqus_Opened.msg so reached end ")


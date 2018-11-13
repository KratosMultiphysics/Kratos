from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

from KratosMultiphysics import *
import KratosMultiphysics.FluidDynamicsApplication
#import KratosMultiphysics.IncompressibleFluidApplication
import KratosMultiphysics.DEMApplication
import KratosMultiphysics.SwimmingDEMApplication
import KratosMultiphysics.ExternalSolversApplication
import KratosSwimmingDEM as script
import json
import os

import sys
print(sys.path)

from . import candelier_analysis

def PrintMessage(run_name, radial_error, tolerance):
        run_name += ': '
        error_message = '    relative radial error: '

        messages_to_print = [run_name, error_message]
        max_len = max([len(msg) for msg in {run_name, error_message}])

        if radial_error != None and radial_error < tolerance:
            veredict_msg = 'OK'
        else:
            veredict_msg = 'Fail'

        error_message = error_message.ljust(max_len) + str(radial_error)
        run_name = run_name.ljust(max_len) + veredict_msg
        print(run_name)
        print(error_message)

def PrintOutput(error_names, errors, tolerance):
    width_buffer = 5
    width = max((len(name) + len(str(error)) for name, error in zip(error_names, errors))) + width_buffer
    thick_line = '=' * width
    separator = '-' * width
    print(thick_line)
    print('Candelier tests results (~ 0.03 = effect of the history force)')

    first = True
    for name, error in zip(error_names, errors):
        if first:
            print(thick_line)
            first = False
        else:
            print(separator)
        PrintMessage(name, error, tolerance)
    print(thick_line)

# Setting parameters
tolerance = 1e-4
errors = []
error_names = []
varying_parameters = dict()
varying_parameters['FinalTime'] = 1

def RunCase(varying_parameters, name):
    model = Model()
    parameters = Parameters(json.dumps(varying_parameters))
    with script.Solution(model, candelier_analysis, parameters) as test:
        error_names.append(name)
        errors.append(test.Run())

# No history force benchmark
varying_parameters['basset_force_type'] = 0
RunCase(varying_parameters, 'No history force, Daitche')

# Second-order accurate Daitche benchmark
varying_parameters['basset_force_type'] = 2
RunCase(varying_parameters, 'All forces, Daitche')

# Rotating frame of reference
varying_parameters['frame_of_reference_type'] = 1
varying_parameters['angular_velocity_of_frame_Z'] = 0.5

# No history force benchmark
varying_parameters['basset_force_type'] = 0
RunCase(varying_parameters, 'No history force, Daitche (rotating frame)')

# Second-order accurate Daitche benchmark
varying_parameters['basset_force_type'] = 2
RunCase(varying_parameters, 'All forces, Daitche (rotating frame)')

# Output
PrintOutput(error_names, errors, tolerance)
from KratosMultiphysics import *
import KratosSwimmingDEM as script
import json
import os
import candelier_algorithm

def PrintMessage(run_name, radial_error, tolerance):
        run_name += ': '
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

tolerance = 1e-4
errors = []
error_names = []
varying_parameters = dict()
varying_parameters['FinalTime'] = 1
varying_parameters['time_steps_per_quadrature_step'] = 1
varying_parameters['number_of_exponentials'] = 10
varying_parameters['number_of_quadrature_steps_in_window'] = 10
varying_parameters['basset_force_type'] = 0

# No history force benchmark
parameters = Parameters(json.dumps(varying_parameters))
with script.Solution(candelier_algorithm, parameters) as test:
    error_names.append('No history force, Daitche')
    errors.append(test.alg.Run())

varying_parameters['basset_force_type'] = 2
# Second-order accurate Daitche benchmark
parameters = Parameters(json.dumps(varying_parameters))
with script.Solution(candelier_algorithm, parameters) as test:
    error_names.append('All forces, Daitche')
    errors.append(test.alg.Run())

# Output
width_buffer = 5
width = max([len(name) + len(str(error)) for name, error in zip(error_names, errors)]) + width_buffer
thick_line = '=' * width
separator = '-' * width
print(thick_line)
print('Candelier tests results')
print(thick_line)

for name, error in zip(error_names, errors):
    PrintMessage(name, error, tolerance)
    print(separator)
print(thick_line)

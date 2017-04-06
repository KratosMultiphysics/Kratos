import KratosSwimmingDEM as script
import os

def PrintMessage(run_name, radial_error, tolerance):
        run_name += ': '
        run_name += ': '
        error_message = 'relative radial error: '

        messages_to_print = [run_name, error_message]
        max_len = max([len(msg) for msg in {run_name, error_message}])

        if radial_error < tolerance:
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
varying_parameters['simulation_time'] = 1
varying_parameters['Nq'] = 1
varying_parameters['m'] = 10
varying_parameters['number_of_quadrature_steps_in_window'] = 10
varying_parameters['basset_force_type'] = 0

# No history force benchmark
import candelier_algorithm
test = script.Solution(candelier_algorithm, varying_parameters)
error_names.append('No history force, Daitche')
errors.append(test.Run())

varying_parameters['basset_force_type'] = 2
# Second-order accurate Daitche benchmark
test = script.Solution(candelier_algorithm, varying_parameters)
error_names.append('All forces, Daitche')
errors.append(test.Run())

# Output
print()
print('-----------------------')
print('Candelier tests results')
print('-----------------------')

for i, e in enumerate(errors):
    PrintMessage(error_names[i], e, tolerance)
print('-----------------------')
print()

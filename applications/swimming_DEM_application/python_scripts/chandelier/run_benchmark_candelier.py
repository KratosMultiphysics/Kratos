import KratosSwimmingDEMCandelier as script
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

# No history force benchmark
test = script.Solution(simulation_time = 1.0, basset_force_type = 0)
error_names.append('No history force, Daitche')
errors.append(test.Run())

# Second-order accurate Daitche benchmark
test = script.Solution(simulation_time = 1.0, basset_force_type = 2)
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

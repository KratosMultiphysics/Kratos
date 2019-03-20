import os
import sys
sys.path.append(os.getcwd())
import json
import fileinput
import tracemalloc
import objgraph
objgraph.show_growth(limit=10)
tracemalloc.start(10)
from KratosMultiphysics import *
import KratosSwimmingDEM as script
# import t_junction_analysis
import pre_calculated_fluid_analysis

varying_parameters = dict()
combinations_that_failed = []
errors = []
varying_parameters["fluid_already_calculated"] = True
varying_parameters["do_search_neighbours"] = False
# varying_parameters["full_particle_history_watcher"] = 'ParticlesHistoryWatcher'
varying_parameters["stationary_problem_option"] = True
varying_parameters["store_full_gradient_option"] = True

set_of_material_acceleration_calculation_types = [1]
L = 0.0048
set_of_inlet_radii = [L * 0.01 * i for i in range(1, 6)]

def ReplaceInletMDPAMeanRadiusValue(value):
    newline ='RADIUS ' + str(value) + '\n'
    for line in fileinput.FileInput('trapping_probabilityDEM_Inlet.mdpa', inplace=1):
        if 'RADIUS' in line:
            result = newline
        else:
            result = line
        sys.stdout.write(result)

# outer_globals = ['gc', 'parameters', 'fileinput', 'json', 'sys', 'os', 'material_acceleration_calculation_type', 'radius', 'varying_parameters', 'combinations_that_failed', 'errors', 'set_of_material_acceleration_calculation_types', 'L', 'set_of_inlet_radii', 'ReplaceInletMDPAMeanRadiusValue', 'outer_globals']
print('problems\n', objgraph.show_growth(limit=1000) )
for radius in set_of_inlet_radii:
    # ReplaceInletMDPAMeanRadiusValue(radius)

    for material_acceleration_calculation_type in set_of_material_acceleration_calculation_types:

        varying_parameters["material_acceleration_calculation_type"] = material_acceleration_calculation_type
        parameters = Parameters(json.dumps(varying_parameters))

        with script.Solution(pre_calculated_fluid_analysis, parameters) as test:
            test.Run()

        # try:
        #     test = script.Solution(algorithm, parameters)
        #     test.Run()
        #     del test
        #     del script
        # except:
        #     error = sys.exc_info()
        #     errors.append(error)
        #     combinations_that_failed.append('Combination: fluid_already_calculated = False')

        print('problems\n', objgraph.show_growth(limit=1000) )
        objgraph.show_backrefs(objgraph.by_type('Algorithm'), filename='chain.png')
        objgraph.show_backrefs(objgraph.by_type('Solution'), filename='chain1.png')

        # for var in [str(var) for var in globals() if not var in outer_globals]:
        #     del globals()[var]
        # print(globals())
        # del SimpleMortarMapperProcess2D2NDoubleHistorical

print()
print('****************************************')

if len(combinations_that_failed):
    print('The following combinations produced an error:')
    print()
    for combination, error in zip(combinations_that_failed, errors):
        print(combination)
        print(error)
else:
    print('All combinations run without errors')
print('****************************************')

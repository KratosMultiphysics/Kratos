import os
import sys
sys.path.append(os.getcwd())
from KratosMultiphysics import *
import KratosSwimmingDEM as script
import sys
import ProjectParameters as pp
import DEM_explicit_solver_var as DEM_parameters
import json
varying_parameters = dict()
combinations_that_failed = []
errors = []
varying_parameters["fluid_already_calculated"] = False
varying_parameters["do_search_neighbours"] = False
varying_parameters["full_particle_history_watcher"] = 'ParticlesHistoryWatcher'
varying_parameters["stationary_problem_option"] = True
parameters = Parameters(json.dumps(varying_parameters))

import t_junction_algorithm
with t_junction_algorithm.Algorithm(parameters) as algorithm:
    try:
        test = script.Solution(algorithm, parameters)
        test.alg.Run()
    except:
        error = sys.exc_info()
        errors.append(error)
        combinations_that_failed.append('Combination: fluid_already_calculated = False')

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

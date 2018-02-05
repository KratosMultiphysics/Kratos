from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *
from KratosMultiphysics.FluidDynamicsApplication import *
from KratosMultiphysics.IncompressibleFluidApplication import *
from KratosMultiphysics.SwimmingDEMApplication import *
import rotating_ale_algorithm

import json

import KratosSwimmingDEM as script

varying_parameters = dict()
combinations_that_failed = []
errors = []
varying_parameters["ALE_option"] = True
varying_parameters["fluid_already_calculated"] = False
from math import pi
varying_parameters["angular_velocity_magnitude"] = - 2 * pi
varying_parameters["frame_rotation_axis_initial_point"] = [0., 0., 0.]
varying_parameters["frame_rotation_axis_final_point"] = [0., 0., 1.]
varying_parameters["print_VISCOSITY_option"] = True
varying_parameters["drag_force_type"] = 13
varying_parameters["virtual_mass_force_type"] = 12
varying_parameters["PostNonDimensionalVolumeWear"] = True
varying_parameters["do_search_neighbours"] = False
varying_parameters["stationary_problem_option"] = True
varying_parameters["time_steps_per_stationarity_step"] = 1
varying_parameters["max_pressure_variation_rate_tol"] = 1

parameters = Parameters(json.dumps(varying_parameters))

with script.Solution(rotating_ale_algorithm, parameters) as test:
    test.Run()

print('\n****************************************')

if len(combinations_that_failed):
    print('The following combinations produced an error:\n')

    for combination, error in zip(combinations_that_failed, errors):
        print(combination)
        print(error)
else:
    print('All combinations run without errors')
print('****************************************\n')

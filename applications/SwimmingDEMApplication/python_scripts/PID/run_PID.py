from KratosMultiphysics import Parameters
import rotating_ale_algorithm

import json

import KratosSwimmingDEM as script

varying_parameters = dict()
combinations_that_failed = []
errors = []
varying_parameters["custom_fluid"]["ALE_option"] = True
varying_parameters["custom_fluid"]["fluid_already_calculated"] = False
from math import pi
varying_parameters['frame_of_reference']["angular_velocity_magnitude"] = - 2 * pi
varying_parameters['frame_of_reference']["frame_rotation_axis_initial_point"] = [0., 0., 0.]
varying_parameters['frame_of_reference']["frame_rotation_axis_final_point"] = [0., 0., 1.]
varying_parameters["print_VISCOSITY_option"] = False
varying_parameters["PostNonDimensionalVolumeWear"] = True
varying_parameters["custom_dem"]["do_search_dem_neighbours"] = True
varying_parameters["stationary_problem_option"] = False
varying_parameters["time_steps_per_stationarity_step"] = 1
varying_parameters["max_pressure_variation_rate_tol"] = 1
varying_parameters["print_MATERIAL_ACCELERATION_option"] = False
varying_parameters["print_MATERIAL_FLUID_ACCEL_PROJECTED_option"] = True
varying_parameters["print_BASSET_FORCE_option"] = True
varying_parameters["print_VORTICITY_option"] = False
varying_parameters["print_VELOCITY_GRADIENT_option"] = False
varying_parameters["initial_averaging_time"] = 0.005
varying_parameters["stationary_start_time"] = 0.0
varying_parameters["coupling"]["interaction_start_time"] = 0.006
varying_parameters["steps_per_average_step"] = 1
varying_parameters["rotated_stationary_flow_option"] = False
varying_parameters["averaging_has_already_been_done"] = False
varying_parameters["do_write_results_to_hdf5"] = False
varying_parameters["stationarity"]["time_steps_per_analytic_processing_step"] = 1000

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

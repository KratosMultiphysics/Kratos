import json
import traceback
import os

from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *
from KratosMultiphysics.FluidDynamicsApplication import *
from KratosMultiphysics.IncompressibleFluidApplication import *
from KratosMultiphysics.SwimmingDEMApplication import *
import ProjectParameters as pp
import DEM_explicit_solver_var as DEM_parameters
from DEM_procedures import KratosPrint as Say
import case_runner
#import ethier_benchmark_algorithm as algorithm
import product_of_sines_benchmark_algorithm as algorithm

varying_parameters = dict()
combinations_that_failed = []
errors = []
irregular_mesh_sizes = [0.4, 0.2, 0.1]
regular_mesh_n_points = [10, 20, 40]
derivatives_types = [1, 3, 4, 6, 7]
number_of_simulations = len(irregular_mesh_sizes)
number_of_simulations += len(regular_mesh_n_points)
number_of_simulations *= len(derivatives_types)
varying_parameters['mesh_tag'] = 'Altair'
varying_parameters['print_VECTORIAL_ERROR_option'] = True

runner = case_runner.CaseRunner(main_path=os.getcwd(),
                                algorithm=algorithm,
                                total_number_of_simulations=number_of_simulations)
simulation_id = 0

for size in irregular_mesh_sizes + regular_mesh_n_points:
    varying_parameters['size_parameter'] = size
    varying_parameters['regular_mesh_option'] = size in regular_mesh_n_points

    for derivatives_type in derivatives_types:
        simulation_id += 1
        varying_parameters['material_acceleration_calculation_type'] = derivatives_type
        varying_parameters['laplacian_calculation_type'] = derivatives_type
        parameters = Parameters(json.dumps(varying_parameters))
        identification_text = 'mesh size: ' + str(size) + '\nrecovery type: ' + str(derivatives_type) + '\n'
        error = runner.RunCase(parameters, simulation_id, identification_text=identification_text)
        if error:
            errors.append(error)
            combinations_that_failed.append({'size':size, 'type':derivatives_type})
Say()
Say('****************************************')

if combinations_that_failed:
    Say('The following combinations produced an error:')
    Say('combinations_that_failed',combinations_that_failed)
    for combination, error in zip(combinations_that_failed, errors):
        Say(combination)
        Say(error[1:2])
        Say(traceback.print_tb(error[2]))
else:
    Say('All combinations run without errors')
Say('****************************************')

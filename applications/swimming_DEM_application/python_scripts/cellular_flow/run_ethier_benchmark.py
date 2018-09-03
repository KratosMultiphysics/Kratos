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

class MeshType:
    number_of_mesh_types = 0
    def __init__(self, is_regular, tag):
        MeshType.number_of_mesh_types += 1
        self.is_regular = is_regular
        self.tag = tag
        self.sizes = []

    def AddSizes(self, sizes):
        self.sizes += [size for size in sizes]

combinations_that_failed = []
errors = []
irregular_mesh_sizes = [0.4, 0.2, 0.1]
regular_mesh_n_points = [40]
derivatives_types = [1, 3, 5, 6, 7]
field_identifiers = ['ethier']
mesh_tags = ['Altair']

mesh_types = []
if irregular_mesh_sizes:
    mesh_type = MeshType(False, '')
    mesh_type.AddSizes(irregular_mesh_sizes)
    mesh_types.append(mesh_type)

if regular_mesh_n_points:
    for tag in mesh_tags:
        mesh_type = MeshType(True, tag)
        mesh_type.AddSizes(regular_mesh_n_points)
        mesh_types.append(mesh_type)

number_of_simulations = sum(len(mt.sizes) for mt in mesh_types)
number_of_simulations *= len(derivatives_types)
number_of_simulations *= len(field_identifiers)
varying_parameters = dict()
varying_parameters['print_VECTORIAL_ERROR_option'] = True

simulation_id = 0

for field in field_identifiers:
    varying_parameters['field_identifier'] = field
    if field == 'ethier':
        import ethier_benchmark_algorithm as algorithm
    elif field == 'sines':
        import product_of_sines_benchmark_algorithm as algorithm
        varying_parameters['field_period'] = 2.0
    runner = case_runner.CaseRunner(
        main_path=os.getcwd(),
        algorithm=algorithm,
        total_number_of_simulations=number_of_simulations)
    for mesh_type in mesh_types:
        varying_parameters['mesh_tag'] = mesh_type.tag
        varying_parameters['regular_mesh_option'] = mesh_type.is_regular
        for size in mesh_type.sizes:
            varying_parameters['size_parameter'] = size
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

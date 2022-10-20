# import kratos core and applications
import KratosMultiphysics as km
import KratosMultiphysics.StructuralMechanicsApplication as ksm
import KratosMultiphysics.TopologyOptimizationApplication as kto
import KratosMultiphysics.LinearSolversApplication as kls
import os
from importlib import import_module
from KratosMultiphysics.gid_output_process import GiDOutputProcess
from KratosMultiphysics.TopologyOptimizationApplication import (
    topology_optimizer_factory,
)
from KratosMultiphysics import process_factory
from KratosMultiphysics.KratosUnittest import TestCase


parameter_file = open('ProjectParameters.json', 'r')
ProjectParameters = km.Parameters(parameter_file.read())
optimization_file = open('Optimization_Parameters.json', 'r')
OptimizationParameters = km.Parameters(optimization_file.read())

# Optimization results to compare for the unittest
results = km.Parameters(
    """
                    {

                    "compliance"                            : 0.0,
                    "number_of_iterations"                  : 0,
                    "volume_fraction"                       : 0.0
                    
                    } 
                    """
)

echo_level = ProjectParameters['problem_data']['echo_level'].GetInt()
current_model = km.Model()
dimension = 2
model_part = current_model.CreateModelPart('Structure')
solver_module = ProjectParameters['solver_settings']['solver_type'].GetString()
mod = 'KratosMultiphysics.TopologyOptimizationApplication.topology_optimization_simp_static_solver'
solver = import_module(mod).CreateSolver(
    current_model, ProjectParameters['solver_settings']
)
solver.AddVariables()
solver.ImportModelPart()
solver.AddDofs()
model_part.GetProperties()[1].SetValue(km.YOUNG_MODULUS, 1)
model_part.GetProperties()[1].SetValue(km.POISSON_RATIO, 0.3)
model_part.GetProperties()[1].SetValue(km.DENSITY, 1)
model_part.GetProperties()[1].SetValue(kto.YOUNGS_MODULUS_MIN, 1e-09)
model_part.GetProperties()[1].SetValue(kto.YOUNGS_MODULUS_0, 1)


cons_law = ksm.LinearElastic3DLaw()
model_part.GetProperties()[1].SetValue(km.CONSTITUTIVE_LAW, cons_law)

if echo_level > 1:
    print(model_part)
    for properties in model_part.Properties:
        print(properties)

# obtain the list of the processes to be applied (the process order of execution is important)
list_of_processes = process_factory.KratosProcessFactory(
    current_model
).ConstructListOfProcesses(
    ProjectParameters['processes']['constraints_process_list']
)
list_of_processes += process_factory.KratosProcessFactory(
    current_model
).ConstructListOfProcesses(
    ProjectParameters['processes']['loads_process_list']
)
if ProjectParameters.Has('problem_process_list'):
    list_of_processes += process_factory.KratosProcessFactory(
        current_model
    ).ConstructListOfProcesses(
        ProjectParameters['processes']['problem_process_list']
    )
if ProjectParameters.Has('output_process_list'):
    list_of_processes += process_factory.KratosProcessFactory(
        current_model
    ).ConstructListOfProcesses(
        ProjectParameters['processes']['output_process_list']
    )

if echo_level > 1:
    for process in list_of_processes:
        print(process)
for process in list_of_processes:
    process.ExecuteInitialize()

computing_model_part = solver.GetComputingModelPart()
problem_path = os.getcwd()
problem_name = ProjectParameters['problem_data']['problem_name'].GetString()

# initialize GiD  I/O (gid outputs, file_lists)
# from gid_output_process import GiDOutputProcess
output_settings = ProjectParameters['processes']['output_configuration']
gid_output_process = GiDOutputProcess(
    computing_model_part, problem_name, output_settings
)

gid_output_process.ExecuteInitialize()
solver.Initialize()
for process in list_of_processes:
    process.ExecuteBeforeSolutionLoop()
gid_output_process.ExecuteBeforeSolutionLoop()


def solve_structure(opt_itr):
    for process in list_of_processes:
        process.ExecuteInitializeSolutionStep()
    gid_output_process.ExecuteInitializeSolutionStep()

    # solve problem
    solver.InitializeSolutionStep()
    solver.SolveSolutionStep()
    solver.FinalizeSolutionStep()

    for process in list_of_processes:
        process.ExecuteFinalizeSolutionStep()
    gid_output_process.ExecuteFinalizeSolutionStep()
    for process in list_of_processes:
        process.ExecuteFinalizeSolutionStep()
    for process in list_of_processes:
        process.ExecuteBeforeOutputStep()
    if gid_output_process.IsOutputStep():
        gid_output_process.PrintOutput()
    for process in list_of_processes:
        process.ExecuteAfterOutputStep()


def FinalizeKSMProcess():
    for process in list_of_processes:
        process.ExecuteFinalize()
    # ending the problem (time integration finished)
    gid_output_process.ExecuteFinalize()


def Analyzer(controls, response, opt_itr):

    # Save the current number of iterations for comparison reasons of the unittest
    results['number_of_iterations'].SetInt(opt_itr)

    # Create object to analyze structure response functions if required
    response_analyzer = kto.StructureResponseFunctionUtilities(model_part)
    linear_solver = km.python_linear_solver_factory.ConstructSolver(
        ProjectParameters['solver_settings']['linear_solver_settings']
    )
    sensitivity_solver = kto.StructureAdjointSensitivityStrategy(
        model_part,
        linear_solver,
        ProjectParameters['solver_settings']['domain_size'].GetInt(),
    )

    # Compute objective function value Call the Solid Mechanics Application to compute objective function value
    if controls['strain_energy']['calc_func']:
        # Compute structure solution to get displacement field u
        solve_structure(opt_itr)
        # Calculate objective function value based on u and save into container
        response['strain_energy'][
            'func'
        ] = response_analyzer.ComputeStrainEnergy()

        # Save the compliance of the given iteration as results for the unittest
        results['compliance'].SetDouble(response['strain_energy']['func'])

    # Compute constraint function value
    if controls['volume_fraction']['calc_func']:
        target_volume_fraction = OptimizationParameters[
            'optimization_parameters'
        ]['initial_volume_fraction'].GetDouble()
        Volume_fraction = response_analyzer.ComputeVolumeFraction()
        response['volume_fraction']['func'] = (
            Volume_fraction - target_volume_fraction
        )

        # Save the volume fraction of the given iteration as results for the unittest
        results['volume_fraction'].SetDouble(Volume_fraction)

    # Compute sensitivities of objective function
    if controls['strain_energy']['calc_grad']:
        sensitivity_solver.ComputeStrainEnergySensitivities()
    # Compute sensitivities of constraint function
    if controls['volume_fraction']['calc_grad']:
        sensitivity_solver.ComputeVolumeFractionSensitivities()


# optimization
optimizer = kto.topology_optimizer_factory.ConstructOptimizer(
    model_part, OptimizationParameters['optimization_parameters'], Analyzer
)
optimizer.optimize()

# Testing the results of the optimization
TestCase().assertEqual(results['number_of_iterations'].GetInt(), 13)
TestCase().assertAlmostEqual(results['compliance'].GetDouble(), 278.493, 3)
TestCase().assertAlmostEqual(
    results['volume_fraction'].GetDouble(), 0.50010, 5
)

FinalizeKSMProcess()

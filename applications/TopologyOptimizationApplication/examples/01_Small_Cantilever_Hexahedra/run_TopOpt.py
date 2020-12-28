#import kratos core and applications
import KratosMultiphysics as km
import KratosMultiphysics.StructuralMechanicsApplication as ksm
import KratosMultiphysics.TopologyOptimizationApplication as kto
import KratosMultiphysics.LinearSolversApplication as kls
from KratosMultiphysics import process_factory
import os
import OptimizationParameters as opt_parameters
from importlib import import_module

parameter_file = open("ProjectParameters.json",'r')
ProjectParameters = km.Parameters(parameter_file.read())
echo_level = ProjectParameters["problem_data"]["echo_level"].GetInt()
model = km.Model()
main_model_part = model.CreateModelPart("Structure")
main_model_part.ProcessInfo.SetValue(km.DOMAIN_SIZE, 3)
solver_module = ProjectParameters["solver_settings"]["solver_type"].GetString()
mod = 'KratosMultiphysics.TopologyOptimizationApplication.topology_optimization_simp_static_solver'
solver = import_module(mod).CreateSolver(main_model_part, ProjectParameters["solver_settings"])
solver.AddVariables()
solver.ImportModelPart()
solver.AddDofs()
for i in range(ProjectParameters["solver_settings"]["processes_sub_model_part_list"].size()):
    part_name = ProjectParameters["solver_settings"]["processes_sub_model_part_list"][i].GetString()
    if( main_model_part.HasSubModelPart(part_name) ):
        Model.update({part_name: main_model_part.GetSubModelPart(part_name)})

if(echo_level>1):
    print(main_model_part)
    for properties in main_model_part.Properties:
        print(properties)

#obtain the list of the processes to be applied (the process order of execution is important)
list_of_processes  = process_factory.KratosProcessFactory(Model).ConstructListOfProcesses( ProjectParameters["constraints_process_list"] )
list_of_processes += process_factory.KratosProcessFactory(Model).ConstructListOfProcesses( ProjectParameters["loads_process_list"] )
if(ProjectParameters.Has("problem_process_list")):
    list_of_processes += process_factory.KratosProcessFactory(Model).ConstructListOfProcesses( ProjectParameters["problem_process_list"] )
if(ProjectParameters.Has("output_process_list")):
    list_of_processes += process_factory.KratosProcessFactory(Model).ConstructListOfProcesses( ProjectParameters["output_process_list"] )
            
if(echo_level>1):
    for process in list_of_processes:
        print(process)
for process in list_of_processes:
    process.ExecuteInitialize()

computing_model_part = solver.GetComputingModelPart()
problem_path = os.getcwd()
problem_name = ProjectParameters["problem_data"]["problem_name"].GetString()

# initialize GiD  I/O (gid outputs, file_lists)
from gid_output_process import GiDOutputProcess
output_settings = ProjectParameters["output_configuration"]
gid_output = GiDOutputProcess(computing_model_part,
                              problem_name,
                              output_settings)

gid_output.ExecuteInitialize()
solver.Initialize()
for process in list_of_processes:
    process.ExecuteBeforeSolutionLoop()
gid_output.ExecuteBeforeSolutionLoop()

def solve_structure(opt_itr):    
    for process in list_of_processes:
        process.ExecuteInitializeSolutionStep()
    gid_output.ExecuteInitializeSolutionStep()       
    solver.Solve()
    for process in list_of_processes:
        process.ExecuteFinalizeSolutionStep()
    gid_output.ExecuteFinalizeSolutionStep()
    for process in list_of_processes:
        process.ExecuteFinalizeSolutionStep()
    for process in list_of_processes:
        process.ExecuteBeforeOutputStep()
    if(gid_output.IsOutputStep()):
        gid_output.PrintOutput()
    for process in list_of_processes:
        process.ExecuteAfterOutputStep() 

def FinalizeKSMProcess():
    for process in list_of_processes:
        process.ExecuteFinalize()
    # ending the problem (time integration finished)
    gid_output.ExecuteFinalize()
    
def Analyzer(controls, response, opt_itr):
    # Create object to analyze structure response functions if required
    response_analyzer = StructureResponseFunctionUtilities(main_model_part)
    linear_solver = kls.linear_solver_factory.ConstructSolver(ProjectParameters["solver_settings"]["linear_solver_settings"])
    sensitivity_solver = StructureAdjointSensitivityStrategy(main_model_part, linear_solver,ProjectParameters["problem_data"]["domain_size"].GetInt())
    # Compute objective function value Call the Solid Mechanics Application to compute objective function value
    if(controls["strain_energy"]["calc_func"]):
        # Compute structure solution to get displacement field u
        solve_structure(opt_itr)
        # Calculate objective function value based on u and save into container
        response["strain_energy"]["func"] = response_analyzer.ComputeStrainEnergy()
    # Compute constraint function value
    if(controls["volume_fraction"]["calc_func"]):
        target_volume_fraction = opt_parameters.initial_volume_fraction
        response["volume_fraction"]["func"] = response_analyzer.ComputeVolumeFraction() - target_volume_fraction
    # Compute sensitivities of objective function
    if(controls["strain_energy"]["calc_grad"]):        
        sensitivity_solver.ComputeStrainEnergySensitivities()
    # Compute sensitivities of constraint function
    if(controls["volume_fraction"]["calc_grad"]):  
        sensitivity_solver.ComputeVolumeFractionSensitivities()

# optimization
optimizer = kto.topology_optimizer_factory.ConstructOptimizer(main_model_part, opt_parameters, Analyzer)
optimizer.optimize()
FinalizeKSMProcess()

from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Time monitoring
import time as timer
print (timer.ctime())
initial_time = timer.perf_counter()

## Importing modules -----------------------------------------------------------------------------------------

# Import system python
import os

# Import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.SolidMechanicsApplication  as KratosSolid
import KratosMultiphysics.ExternalSolversApplication as KratosSolvers
import KratosMultiphysics.PoromechanicsApplication as KratosPoro

# Parsing the parameters
parameter_file = open("ProjectParameters.json",'r')
ProjectParameters = KratosMultiphysics.Parameters( parameter_file.read())

#Import solver module
solver_module = __import__(ProjectParameters["solver_settings"]["solver_type"].GetString())

# Import process modules
import process_factory
from gid_output_process import GiDOutputProcess

# Import utilities
import poromechanics_fracture_propagation_utility
import poromechanics_cleaning_utility

## Defining variables ----------------------------------------------------------------------------------------

# Number of threads
parallel=KratosMultiphysics.OpenMPUtils()
parallel.SetNumThreads(ProjectParameters["problem_data"]["OMP_threads"].GetInt())

# Problem variables
domain_size = ProjectParameters["problem_data"]["domain_size"].GetInt()
problem_name = ProjectParameters["problem_data"]["problem_name"].GetString()
problem_path = os.getcwd()
echo_level = ProjectParameters["solver_settings"]["echo_level"].GetInt()
buffer_size = ProjectParameters["solver_settings"]["buffer_size"].GetInt()
delta_time = ProjectParameters["problem_data"]["time_step"].GetDouble()
end_time = ProjectParameters["problem_data"]["end_time"].GetDouble()
time = ProjectParameters["problem_data"]["start_time"].GetDouble()
tol = delta_time*1.0e-10
output_settings = ProjectParameters["output_configuration"]
FracturePropagation = ProjectParameters["solver_settings"]["fracture_propagation"].GetBool()

## Model part ------------------------------------------------------------------------------------------------

# Defining the model part
main_model_part = KratosMultiphysics.ModelPart(ProjectParameters["problem_data"]["model_part_name"].GetString())
main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, domain_size)
main_model_part.ProcessInfo.SetValue(KratosMultiphysics.TIME, time)
main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DELTA_TIME, delta_time)
main_model_part.ProcessInfo.SetValue(KratosPoro.TIME_UNIT_CONVERTER, 1.0)
Model = {ProjectParameters["problem_data"]["model_part_name"].GetString() : main_model_part}

# Construct the solver (main setting methods are located in the solver_module)
solver = solver_module.CreateSolver(main_model_part, ProjectParameters["solver_settings"])

# Add problem variables
solver.AddVariables()

# Read model_part (note: the buffer_size is set here)
solver.ImportModelPart()

# Add degrees of freedom
solver.AddDofs()

# Build sub_model_parts (save the list of the submodel part in the object Model)
for i in range(ProjectParameters["solver_settings"]["processes_sub_model_part_list"].size()):
    part_name = ProjectParameters["solver_settings"]["processes_sub_model_part_list"][i].GetString()
    Model.update({part_name : main_model_part.GetSubModelPart(part_name)})

# Print model_part and properties
if(echo_level > 1):
    print(main_model_part)
    for properties in main_model_part.Properties:
        print(properties)


## Initialize ------------------------------------------------------------------------------------------------

# Construct processes to be applied
list_of_processes = process_factory.KratosProcessFactory(Model).ConstructListOfProcesses( ProjectParameters["constraints_process_list"] )
list_of_processes += process_factory.KratosProcessFactory(Model).ConstructListOfProcesses( ProjectParameters["loads_process_list"] )

# Print list of constructed processes
if(echo_level>1):
    for process in list_of_processes:
        print(process)

# Initialize processes
for process in list_of_processes:
    process.ExecuteInitialize()

# Set TIME and DELTA_TIME and fill the previous steps of the buffer with the initial conditions
time = time - (buffer_size-1)*delta_time
main_model_part.ProcessInfo[KratosMultiphysics.TIME] = time
for step in range(buffer_size-1):
    time = time + delta_time
    main_model_part.CloneTimeStep(time)

# Clean previous post files
poromechanics_cleaning_utility.CleanPreviousFiles(problem_path)

# Initialize GiD I/O
computing_model_part = solver.GetComputingModelPart()
gid_output = GiDOutputProcess(computing_model_part,problem_name,output_settings)
gid_output.ExecuteInitialize()

# Initialize the solver
solver.Initialize()

# ExecuteBeforeSolutionLoop
for process in list_of_processes:
    process.ExecuteBeforeSolutionLoop()
    
## Set results when they are written in a single file
gid_output.ExecuteBeforeSolutionLoop()

# Initialize Fracture Propagation Utility
if FracturePropagation:
    fracture_utility = poromechanics_fracture_propagation_utility.FracturePropagationUtility(domain_size,problem_name)

## Temporal loop ---------------------------------------------------------------------------------------------

while( (time+tol) <= end_time ):
    
    # Update temporal variables
    delta_time = main_model_part.ProcessInfo[KratosMultiphysics.DELTA_TIME]
    time = time + delta_time
    main_model_part.CloneTimeStep(time)
    
    # Update imposed conditions
    for process in list_of_processes:
        process.ExecuteInitializeSolutionStep()

    gid_output.ExecuteInitializeSolutionStep()
    
    # Solve step
    solver.Solve()
    
    gid_output.ExecuteFinalizeSolutionStep()
    
    for process in list_of_processes:
        process.ExecuteFinalizeSolutionStep()
    
    for process in list_of_processes:
        process.ExecuteBeforeOutputStep()
    
    # Write GiD results
    if gid_output.IsOutputStep():
        gid_output.PrintOutput()
    
    for process in list_of_processes:
        process.ExecuteAfterOutputStep()
    
    # Fracture Propagation Utility
    IsRemeshed = False
    if FracturePropagation:
        if fracture_utility.IsPropagationStep():
            main_model_part,solver,list_of_processes,gid_output,IsRemeshed = fracture_utility.CheckPropagation(main_model_part,
                                                                                                               solver,
                                                                                                               list_of_processes,
                                                                                                               gid_output)

## Finalize --------------------------------------------------------------------------------------------------

# Print last step if model was regenerated at the end
if IsRemeshed:
    #solver._CheckConvergence()
    gid_output.ExecuteInitializeSolutionStep()
    gid_output.ExecuteFinalizeSolutionStep()
    gid_output.PrintOutput()

# Finalizing output files
gid_output.ExecuteFinalize()

for process in list_of_processes:
    process.ExecuteFinalize()
    
# Finalizing strategy
solver.Clear()

# Time control
print("Analysis Completed. Elapsed Time = %.3f" % (timer.perf_counter() - initial_time)," seconds.")
print (timer.ctime())

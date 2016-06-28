from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Time monitoring
import time as timer
print (timer.ctime())
initial_time = timer.perf_counter()


## Importing modules -----------------------------------------------------------------------------------------

# Import Kratos
from KratosMultiphysics import *
# Import Applications
from KratosMultiphysics.SolidMechanicsApplication import *
from KratosMultiphysics.ExternalSolversApplication import *
from KratosMultiphysics.PoromechanicsApplication import *

# Parsing the parameters
parameter_file = open("ProjectParameters.json",'r')
ProjectParameters = Parameters( parameter_file.read())

#Import solver module
solver_module = __import__(ProjectParameters["solver_settings"]["solver_type"].GetString())

# Import process modules
import process_factory
from gid_output_process import GiDOutputProcess


## Defining variables ----------------------------------------------------------------------------------------

# Number of threads
parallel=OpenMPUtils()
parallel.SetNumThreads(ProjectParameters["problem_data"]["OMP_threads"].GetInt())

# Problem variables
problem_name = ProjectParameters["problem_data"]["problem_name"].GetString()
#problem_path = os.getcwd()
echo_level = ProjectParameters["problem_data"]["echo_level"].GetInt()
buffer_size = ProjectParameters["solver_settings"]["buffer_size"].GetInt()
delta_time = ProjectParameters["problem_data"]["time_step"].GetDouble()
end_time = ProjectParameters["problem_data"]["end_time"].GetDouble()
step = 0
time = ProjectParameters["problem_data"]["start_time"].GetDouble()
tol = delta_time*1.0e-10
output_settings = ProjectParameters["output_configuration"]

## Model part ------------------------------------------------------------------------------------------------

# Defining the model part
main_model_part = ModelPart(ProjectParameters["problem_data"]["model_part_name"].GetString())
main_model_part.ProcessInfo.SetValue(DOMAIN_SIZE, ProjectParameters["problem_data"]["domain_size"].GetInt())
#TODO replace this "model" for real one once available in kratos core
Model = {ProjectParameters["problem_data"]["model_part_name"].GetString() : main_model_part}

# Construct the solver (main setting methods are located in the solver_module)
solver = solver_module.CreateSolver(main_model_part, ProjectParameters["solver_settings"])

# Add problem variables
solver.AddVariables()

# Read model_part (note: the buffer_size is set here)
solver.ImportModelPart()

# Add degrees of freedom
solver.AddDofs()

# Build sub_model_parts
for i in range(ProjectParameters["solver_settings"]["processes_sub_model_part_list"].size()):
    part_name = ProjectParameters["solver_settings"]["processes_sub_model_part_list"][i].GetString()
    Model.update({part_name : main_model_part.GetSubModelPart(part_name)})

# Print control
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
main_model_part.ProcessInfo[TIME] = time
for step in range(buffer_size-1):
    time = time + delta_time
    main_model_part.CloneTimeStep(time)

# Initialize GiD I/O
computing_model_part = solver.GetComputeModelPart()
gid_output = GiDOutputProcess(computing_model_part,problem_name,output_settings)
gid_output.ExecuteInitialize()

# Initialize the solver
solver.Initialize()

# ExecuteBeforeSolutionLoop
for process in list_of_processes:
    process.ExecuteBeforeSolutionLoop()
    
## Set results when they are written in a single file
gid_output.ExecuteBeforeSolutionLoop()

## Temporal loop ---------------------------------------------------------------------------------------------

while( (time+tol) <= end_time ):
    
    # Update temporal variables
    time = time + delta_time
    step = step + 1
    main_model_part.CloneTimeStep(time)
    
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

## Finalize --------------------------------------------------------------------------------------------------

# Finalizing output files
gid_output.ExecuteFinalize()

for process in list_of_processes:
    process.ExecuteFinalize()
    
# Finalizing strategy
solver.Clear()

# Time control
print("Analysis Completed. Elapsed Time = %.3f" % (timer.perf_counter() - initial_time)," seconds.")
print (timer.ctime())

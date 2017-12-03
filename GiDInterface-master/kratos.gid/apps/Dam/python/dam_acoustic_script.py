from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Time monitoring
import time as timer
print(timer.ctime())
initial_time = timer.perf_counter()

## Necessary modules -----------------------------------------------------------------------------------------

# Import system python
import os

# Including kratos path
import KratosMultiphysics
# Including Applications path
import KratosMultiphysics.ExternalSolversApplication as KratosSolvers
import KratosMultiphysics.ConvectionDiffusionApplication as KratosConvDiff
import KratosMultiphysics.SolidMechanicsApplication as KratosSolid
import KratosMultiphysics.PoromechanicsApplication as KratosPoro
import KratosMultiphysics.DamApplication as KratosDam

# Parsing the parameters
parameter_file = open("ProjectParameters.json",'r')
ProjectParameters = KratosMultiphysics.Parameters( parameter_file.read())

# Parallel Configuration
parallel_type = ProjectParameters["problem_data"]["parallel_type"].GetString()
parallel=KratosMultiphysics.OpenMPUtils()
parallel.SetNumThreads(ProjectParameters["problem_data"]["number_of_threads"].GetInt())
if parallel_type == "MPI":
    import KratosMultiphysics.mpi as KratosMPI
    print("MPI parallel configuration. OMP_NUM_THREADS =",parallel.GetNumThreads())
else:
    print("OpenMP parallel configuration. OMP_NUM_THREADS =",parallel.GetNumThreads())


## Defining variables ----------------------------------------------------------------------------------------
problem_path = os.getcwd()
problem_name = ProjectParameters["problem_data"]["problem_name"].GetString()
domain_size = ProjectParameters["problem_data"]["domain_size"].GetInt()
echo_level = ProjectParameters["solver_settings"]["echo_level"].GetInt()
buffer_size = ProjectParameters["solver_settings"]["buffer_size"].GetInt()
delta_time = ProjectParameters["problem_data"]["time_step"].GetDouble()
end_time = ProjectParameters["problem_data"]["end_time"].GetDouble()
time = ProjectParameters["problem_data"]["start_time"].GetDouble()
tol = delta_time*1.0e-10
time_scale = ProjectParameters["problem_data"]["time_scale"].GetString()
# Time Units Converter
if(time_scale=="Months"):               # Factor to pass from months to seconds
    time_unit_converter = 2592000.0
elif(time_scale=="Days"):               # Factor to pass from days to seconds
    time_unit_converter = 86400.0
elif(time_scale=="Hours"):              # Factor to pass from hours to seconds
    time_unit_converter = 3600.0
else:                                       # No changes
    time_unit_converter = 1.0

# Update time variables
delta_time = delta_time * time_unit_converter
end_time = end_time * time_unit_converter
time = time * time_unit_converter
tol = tol * time_unit_converter

## Model part ------------------------------------------------------------------------------------------------

# Defining the model part

main_model_part = KratosMultiphysics.ModelPart(ProjectParameters["problem_data"]["model_part_name"].GetString())
main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, domain_size)
main_model_part.ProcessInfo.SetValue(KratosMultiphysics.TIME, time)
main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DELTA_TIME, delta_time)
main_model_part.ProcessInfo.SetValue(KratosPoro.TIME_UNIT_CONVERTER, time_unit_converter)

# Construct the solver (main setting methods are located in the solver_module)
solver_module = __import__("dam_P_solver")
solver = solver_module.CreateSolver(main_model_part, ProjectParameters["solver_settings"])

# Add variables
solver.AddVariables()

# Read model_part (note: the buffer_size is set here)
solver.ImportModelPart()

# Add degrees of freedom
solver.AddDofs()

# Creation of Kratos model
DamModel = KratosMultiphysics.Model()
DamModel.AddModelPart(main_model_part)

# Build sub_model_parts or submeshes (rearrange parts for the application of custom processes)
## Get the list of the submodel part in the object Model
for i in range(ProjectParameters["solver_settings"]["processes_sub_model_part_list"].size()):
    part_name = ProjectParameters["solver_settings"]["processes_sub_model_part_list"][i].GetString()
    DamModel.AddModelPart(main_model_part.GetSubModelPart(part_name))

# Print control
if(echo_level > 1):
    print(main_model_part)
    for properties in main_model_part.Properties:
        print(properties)


## Initialize ------------------------------------------------------------------------------------------------

# Construct the processes to be applied
import process_factory
list_of_processes =  process_factory.KratosProcessFactory(DamModel).ConstructListOfProcesses( ProjectParameters["constraints_process_list"] )
list_of_processes += process_factory.KratosProcessFactory(DamModel).ConstructListOfProcesses( ProjectParameters["loads_process_list"] )

# Print list of constructed processes
if(echo_level>1):
    for process in list_of_processes:
        print(process)

# Initialize processes
for process in list_of_processes:
    process.ExecuteInitialize()

# Set TIME and DELTA_TIME and fill the previous steps of the buffer with the initial conditions
time = time - (buffer_size-1)*delta_time
main_model_part.ProcessInfo.SetValue(KratosMultiphysics.TIME, time)
for step in range(buffer_size-1):
    time = time + delta_time
    main_model_part.CloneTimeStep(time)

# Initialize GiD I/O
computing_model_part = solver.GetComputingModelPart()
output_settings = ProjectParameters["output_configuration"]
if parallel_type == "OpenMP":
    import poromechanics_cleaning_utility
    poromechanics_cleaning_utility.CleanPreviousFiles(problem_path) # Clean previous post files
    from gid_output_process import GiDOutputProcess
    gid_output = GiDOutputProcess(computing_model_part,
                                  problem_name,
                                  output_settings)
else:
    from gid_output_process_mpi import GiDOutputProcessMPI
    gid_output = GiDOutputProcessMPI(computing_model_part,
                                     problem_name,
                                     output_settings)
gid_output.ExecuteInitialize()

#Initialize the  solver
solver.Initialize()

# ExecuteBeforeSolutionLoop
for process in list_of_processes:
    process.ExecuteBeforeSolutionLoop()

## Set results when they are written in a single file
gid_output.ExecuteBeforeSolutionLoop()

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


## Finalize --------------------------------------------------------------------------------------------------

# Finalizing output files
gid_output.ExecuteFinalize()

for process in list_of_processes:
    process.ExecuteFinalize()

# Finalizing strategy
if parallel_type == "OpenMP":
    solver.Clear()

# Time control
print("Analysis Completed. Elapsed Time = %.3f" % (timer.perf_counter() - initial_time)," seconds.")
print(timer.ctime())

from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
import KratosMultiphysics.ContactStructuralMechanicsApplication as ContactStructuralMechanicsApplication
import KratosMultiphysics.ExternalSolversApplication as ExternalSolversApplication
import os

#### TIME MONITORING START ####

# time control starts
import time as timer
print(timer.ctime())
# measure process time
t0p = timer.clock()
# measure wall time
t0w = timer.time()

#
def StartTimeMeasuring():
    # measure process time
    time_ip = timer.clock()
    return time_ip

def StopTimeMeasuring(time_ip, process, report):
    # measure process time
    time_fp = timer.clock()
    if( report ):
        used_time = time_fp - time_ip
        print("::[KSM Simulation]:: [ %.2f" % round(used_time,2),"s", process," ] ")

#### TIME MONITORING END ####

#### PARSING THE PARAMETERS ####

# Import define_output
parameter_file = open("ProjectParameters.json",'r')
ProjectParameters = KratosMultiphysics.Parameters( parameter_file.read())

# Set echo level
echo_level = ProjectParameters["problem_data"]["echo_level"].GetInt()

#### Model_part settings start ####
# Defining the model_part
main_model_part = KratosMultiphysics.ModelPart(ProjectParameters["problem_data"]["model_part_name"].GetString())
main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, ProjectParameters["problem_data"]["domain_size"].GetInt())

Model = {ProjectParameters["problem_data"]["model_part_name"].GetString(): main_model_part}

# Construct the solver (main setting methods are located in the solver_module)
solver_module = __import__(ProjectParameters["solver_settings"]["solver_type"].GetString())
solver = solver_module.CreateSolver(main_model_part, ProjectParameters["solver_settings"])

# Add variables (always before importing the model part) (it must be integrated in the ImportModelPart)
# If we integrate it in the model part we cannot use combined solvers
solver.AddVariables()

# Read model_part (note: the buffer_size is set here) (restart can be read here)
solver.ImportModelPart()

# Add dofs (always after importing the model part) (it must be integrated in the ImportModelPart)
# If we integrate it in the model part we cannot use combined solvers
solver.AddDofs()

# Build sub_model_parts or submeshes (rearrange parts for the application of custom processes)
# #Get the list of the submodel part in the object Model
for i in range(ProjectParameters["solver_settings"]["processes_sub_model_part_list"].size()):
    part_name = ProjectParameters["solver_settings"]["processes_sub_model_part_list"][i].GetString()
    Model.update({part_name: main_model_part.GetSubModelPart(part_name)})

# Print model_part and properties
if(echo_level>1):
    print("")
    print(main_model_part)
    for properties in main_model_part.Properties:
        print(properties)

# Obtain the list of the processes to be applied
import process_factory
list_of_processes = process_factory.KratosProcessFactory(Model).ConstructListOfProcesses(ProjectParameters["constraints_process_list"])
list_of_processes += process_factory.KratosProcessFactory(Model).ConstructListOfProcesses(ProjectParameters["loads_process_list"])
if (ProjectParameters.Has("list_other_processes") == True): 
    list_of_processes += process_factory.KratosProcessFactory(Model).ConstructListOfProcesses(ProjectParameters["list_other_processes"])
if (ProjectParameters.Has("json_output_process") == True): 
    list_of_processes += process_factory.KratosProcessFactory(Model).ConstructListOfProcesses(ProjectParameters["json_output_process"]) 

for process in list_of_processes:
    process.ExecuteInitialize()

# ### START SOLUTION ####

computing_model_part = solver.GetComputingModelPart()

# ### Output settings start ####
problem_path = os.getcwd()
problem_name = ProjectParameters["problem_data"]["problem_name"].GetString()

# ### Output settings start ####
output_post = ProjectParameters.Has("output_configuration")
if (output_post == True):
    from gid_output_process import GiDOutputProcess
    output_settings = ProjectParameters["output_configuration"]
    gid_output = GiDOutputProcess(computing_model_part,
                                        problem_name,
                                        output_settings)
    gid_output.ExecuteInitialize()
    
# Sets strategies, builders, linear solvers, schemes and solving info, and fills the buffer
solver.Initialize()
solver.SetEchoLevel(echo_level)

if (output_post == True):
    gid_output.ExecuteBeforeSolutionLoop()

for process in list_of_processes:
    process.ExecuteBeforeSolutionLoop()

# #Stepping and time settings (get from process info or solving info)
# Delta time
delta_time = ProjectParameters["problem_data"]["time_step"].GetDouble()
# Start step
main_model_part.ProcessInfo[KratosMultiphysics.TIME_STEPS] = 0
# Start time
time = ProjectParameters["problem_data"]["start_time"].GetDouble()
# End time
end_time = ProjectParameters["problem_data"]["end_time"].GetDouble()

# Solving the problem (time integration)
while(time <= end_time):
    time = time + delta_time
    main_model_part.ProcessInfo[KratosMultiphysics.TIME_STEPS] += 1
    main_model_part.CloneTimeStep(time)

    for process in list_of_processes:
        process.ExecuteInitializeSolutionStep()
        
    if (output_post == True):
        gid_output.ExecuteInitializeSolutionStep()
                
    solver.Clear()
    solver.Solve()

    current_vals = [ev for ev in computing_model_part.ProcessInfo[StructuralMechanicsApplication.EIGENVALUE_VECTOR]]
    print(current_vals)
    
    if (output_post == True):
        gid_output.ExecuteFinalizeSolutionStep()

    for process in list_of_processes:
        process.ExecuteFinalizeSolutionStep()

    for process in list_of_processes:
        process.ExecuteBeforeOutputStep()

    for process in list_of_processes:
        process.ExecuteAfterOutputStep()

    if (output_post == True):
        if gid_output.IsOutputStep():
            gid_output.PrintOutput()

if (output_post == True):
    gid_output.ExecuteFinalize()

for process in list_of_processes:
    process.ExecuteFinalize()

print("::[KSM Simulation]:: Analysis -END- ")
print(" ")

# Check solving information for any problem
#~ solver.InfoCheck() # InfoCheck not implemented yet.

#### END SOLUTION ####

# measure process time
tfp = timer.clock()
# measure wall time
tfw = timer.time()

print("::[KSM Simulation]:: [ Computing Time = (%.2f" % (tfp - t0p)," seconds process time) ( %.2f" % (tfw - t0w)," seconds wall time) ]")

print(timer.ctime())

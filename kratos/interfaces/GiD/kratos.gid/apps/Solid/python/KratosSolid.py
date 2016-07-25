from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

#### TIME MONITORING START ####

# Time control starts
import time as timer
print(timer.ctime())
# Measure process time
t0p = timer.clock()
# Measure wall time
t0w = timer.time()

def StartTimeMeasuring():
    # Measure process time
    time_ip = timer.clock()
    return time_ip

def StopTimeMeasuring(time_ip, process, report):
    # Measure process time
    time_fp = timer.clock()
    if( report ):
        used_time = time_fp - time_ip
        print("::[KSM Simulation]:: [ %.2f" % round(used_time,2),"s", process," ] ")

#### TIME MONITORING END ####

import os

# Import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.SolidMechanicsApplication  as KratosSolid
import KratosMultiphysics.ExternalSolversApplication as KratosSolvers

######################################################################################
######################################################################################
######################################################################################

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
main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DELTA_TIME, ProjectParameters["problem_data"]["time_step"].GetDouble())
main_model_part.ProcessInfo.SetValue(KratosMultiphysics.TIME, ProjectParameters["problem_data"]["start_time"].GetDouble())

print("DeltaTime", main_model_part.ProcessInfo[KratosMultiphysics.DELTA_TIME])

Model = {ProjectParameters["problem_data"]["model_part_name"].GetString() : main_model_part}

#construct the solver (main setting methods are located in the solver_module)
solver_module = __import__(ProjectParameters["solver_settings"]["solver_type"].GetString())
solver = solver_module.CreateSolver(main_model_part, ProjectParameters["solver_settings"])

# Add variables (always before importing the model part) (it must be integrated in the ImportModelPart)
# If we integrate it in the model part we cannot use combined solvers
solver.AddVariables()

# Read model_part (note: the buffer_size is set here) (restart can be read here)
solver.ImportModelPart()

# Add dofs (always after importing the model part) (it must be integrated in the ImportModelPart)
# If we integrate it in the model part we cannot use combined solvers
if((main_model_part.ProcessInfo).Has(KratosMultiphysics.IS_RESTARTED)):
    if(main_model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED] == False):
        solver.AddDofs()
else:
    solver.AddDofs()
    
# Build sub_model_parts or submeshes (rearrange parts for the application of custom processes)
## Get the list of the submodel part in the object Model
for i in range(ProjectParameters["solver_settings"]["processes_sub_model_part_list"].size()):
    part_name = ProjectParameters["solver_settings"]["processes_sub_model_part_list"][i].GetString()
    Model.update({part_name: main_model_part.GetSubModelPart(part_name)})

#print model_part and properties
if(echo_level>1):
    print("")
    print(main_model_part)
    for properties in main_model_part.Properties:
        print(properties)

#### Model_part settings end ####


#### Processes settings start ####

#obtain the list of the processes to be applied

import process_factory
#the process order of execution is important
list_of_processes  = process_factory.KratosProcessFactory(Model).ConstructListOfProcesses( ProjectParameters["constraints_process_list"] )
list_of_processes += process_factory.KratosProcessFactory(Model).ConstructListOfProcesses( ProjectParameters["loads_process_list"] )
if(ProjectParameters.Has("problem_process_list")):
    list_of_processes += process_factory.KratosProcessFactory(Model).ConstructListOfProcesses( ProjectParameters["problem_process_list"] )
if(ProjectParameters.Has("output_process_list")):
    list_of_processes += process_factory.KratosProcessFactory(Model).ConstructListOfProcesses( ProjectParameters["output_process_list"] )
            
#print list of constructed processes
if(echo_level>1):
    for process in list_of_processes:
        print(process)

for process in list_of_processes:
    process.ExecuteInitialize()

#### processes settings end ####

#### START SOLUTION ####

computing_model_part = solver.GetComputeModelPart()

#### Output settings start ####

problem_path = os.getcwd()
problem_name = ProjectParameters["problem_data"]["problem_name"].GetString()

# Initialize GiD  I/O (gid outputs, file_lists)
from gid_output_process import GiDOutputProcess
output_settings = ProjectParameters["output_configuration"]
gid_output = GiDOutputProcess(computing_model_part,
                              problem_name,
                              output_settings)

gid_output.ExecuteInitialize()

#### Output settings end ####

## Sets strategies, builders, linear solvers, schemes and solving info, and fills the buffer
solver.Initialize()
solver.SetEchoLevel(echo_level)

print(" ")
print("::[KSM Simulation]:: Analysis -START- ")

for process in list_of_processes:
    process.ExecuteBeforeSolutionLoop()

# writing a initial state results file or single file (if no restart)
if((main_model_part.ProcessInfo).Has(KratosMultiphysics.IS_RESTARTED)):
    if(main_model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED] == False):
        gid_output.ExecuteBeforeSolutionLoop()


# Set time settings
step       = main_model_part.ProcessInfo[KratosMultiphysics.STEP]
time       = main_model_part.ProcessInfo[KratosMultiphysics.TIME]

end_time   = ProjectParameters["problem_data"]["end_time"].GetDouble()
delta_time = ProjectParameters["problem_data"]["time_step"].GetDouble()


# Solving the problem (time integration)
while(time <= end_time):

    # current time parameters
    # main_model_part.ProcessInfo.GetPreviousStepInfo(1)[KratosMultiphysics.DELTA_TIME] = delta_time
    delta_time = main_model_part.ProcessInfo[KratosMultiphysics.DELTA_TIME]

    time = time + delta_time
    step = step + 1

    main_model_part.ProcessInfo[KratosMultiphysics.STEP] = step
    main_model_part.CloneTimeStep(time) 

    # solve 
    for process in list_of_processes:
        process.ExecuteInitializeSolutionStep()

    gid_output.ExecuteInitializeSolutionStep()
        
    solver.Solve()

    gid_output.ExecuteFinalizeSolutionStep()
        
    for process in list_of_processes:
        process.ExecuteFinalizeSolutionStep()

    for process in list_of_processes:
        process.ExecuteBeforeOutputStep()
    
    # write results: (frequency writing is controlled internally by the gid_io)
    if(gid_output.IsOutputStep()):
        gid_output.PrintOutput()
                      
    for process in list_of_processes:
        process.ExecuteAfterOutputStep()


# Ending the problem (time integration finished)
gid_output.ExecuteFinalize()

for process in list_of_processes:
    process.ExecuteFinalize()

print("::[KSM Simulation]:: Analysis -END- ")
print(" ")

# Check solving information for any problem
#~ solver.InfoCheck() # InfoCheck not implemented yet.

#### END SOLUTION ####

# Measure process time
tfp = timer.clock()
# Measure wall time
tfw = timer.time()

print("::[KSM Simulation]:: [ Computing Time = (%.2f" % (tfp - t0p)," seconds process time) ( %.2f" % (tfw - t0w)," seconds wall time) ]")

print(timer.ctime())





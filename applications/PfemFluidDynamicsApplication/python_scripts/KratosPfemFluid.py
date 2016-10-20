from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#

import define_output

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

#### SET NUMBER OF THREADS ####

def SetParallelSize(num_threads):
    parallel = OpenMPUtils()
    parallel.SetNumThreads(int(num_threads))
    print("::[KPFEM Simulation]:: [OMP USING",num_threads,"THREADS ]")
    #parallel.PrintOMPInfo()
    print(" ")

#
#
import sys
#sys.path.append(ProjectParameters.kratos_path)
import os

from KratosMultiphysics import *
from KratosMultiphysics.ExternalSolversApplication import *
from KratosMultiphysics.SolidMechanicsApplication import *
from KratosMultiphysics.PfemBaseApplication import *
from KratosMultiphysics.PfemFluidDynamicsApplication import *
from KratosMultiphysics.PfemSolidMechanicsApplication import *

# defining variables to be used

variables_dictionary = {"PRESSURE": PRESSURE,
                        "VELOCITY": VELOCITY,
                        "REACTION": REACTION,
                        "DISTANCE": DISTANCE,
                        "FREESURFACE": FREESURFACE,
                        "INTERF": INTERF}


#import define_output
parameter_file = open("ProjectParameters.json",'r')
PProjectParameters =  Parameters(parameter_file.read())

##set echo level
echo_level = PProjectParameters["problem_data"]["echo_level"].GetInt()

print(" ")

# defining the number of threads:
threads = PProjectParameters["problem_data"]["threads"].GetInt()
SetParallelSize(threads)

print(" ")
print("::[KPFEM Simulation]:: [Time Step:", PProjectParameters["problem_data"]["time_step"].GetDouble()," echo:", echo_level,"]")


main_model_part = ModelPart(PProjectParameters["problem_data"]["model_part_name"].GetString())

main_model_part.ProcessInfo.SetValue(DOMAIN_SIZE, PProjectParameters["problem_data"]["domain_size"].GetInt())
main_model_part.ProcessInfo.SetValue(DELTA_TIME, PProjectParameters["problem_data"]["time_step"].GetDouble())
main_model_part.ProcessInfo.SetValue(TIME, PProjectParameters["problem_data"]["start_time"].GetDouble())

Model = {PProjectParameters["problem_data"]["model_part_name"].GetString() : main_model_part}

#construct the solver (main setting methods are located in the solver_module)
solver_module = __import__(PProjectParameters["solver_settings"]["solver_type"].GetString())
solver = solver_module.CreateSolver(main_model_part, PProjectParameters["solver_settings"])

main_model_part.AddNodalSolutionStepVariable(FREESURFACE)
main_model_part.AddNodalSolutionStepVariable(INTERF)
main_model_part.AddNodalSolutionStepVariable(NODAL_H)
main_model_part.AddNodalSolutionStepVariable(SHRINK_FACTOR)
main_model_part.AddNodalSolutionStepVariable(MEAN_ERROR)
main_model_part.AddNodalSolutionStepVariable(CONTACT_FORCE)
main_model_part.AddNodalSolutionStepVariable(RIGID_WALL)

for node in main_model_part.Nodes:
    node.SetNodalStepVariables(FREESURFACE, 0)
    node.SetNodalStepVariables(INTERF, 0)

solver.AddVariables()

# Read model_part (note: the buffer_size is set here) (restart is read here)
solver.ImportModelPart()

solver.AddDofs()


# Build sub_model_parts or submeshes (rearrange parts for the application of custom processes)
## Get the list of the submodel part in the object Model
for i in range(PProjectParameters["solver_settings"]["processes_sub_model_part_list"].size()):
    part_name = PProjectParameters["solver_settings"]["processes_sub_model_part_list"][i].GetString()
    Model.update({part_name: main_model_part.GetSubModelPart(part_name)})

#### Processes settings start ####

#obtain the list of the processes to be applied

import process_factory
#the process order of execution is important
list_of_processes  = process_factory.KratosProcessFactory(Model).ConstructListOfProcesses( PProjectParameters["constraints_process_list"] )
if(PProjectParameters.Has("load_process_list")):
    list_of_processes += process_factory.KratosProcessFactory(Model).ConstructListOfProcesses( PProjectParameters["loads_process_list"] )
if(PProjectParameters.Has("problem_process_list")):
    list_of_processes += process_factory.KratosProcessFactory(Model).ConstructListOfProcesses( PProjectParameters["problem_process_list"] )
if(PProjectParameters.Has("output_process_list")):
    list_of_processes += process_factory.KratosProcessFactory(Model).ConstructListOfProcesses( PProjectParameters["output_process_list"] )

if(echo_level>0):
    for process in list_of_processes:
        print(process)


for process in list_of_processes:
    print("::[fluid PFEM Simulation]:: process.ExecuteInitialize()")
    process.ExecuteInitialize()


#print list of constructed processes
#### START SOLUTION ####

computing_model_part = solver.GetComputingModelPart()


#### Output settings start ####

problem_path = os.getcwd()
problem_name = PProjectParameters["problem_data"]["problem_name"].GetString()

# Initialize GiD  I/O (gid outputs, file_lists)
from gid_output_process import GiDOutputProcess
output_settings = PProjectParameters["output_configuration"]
gid_output = GiDOutputProcess(computing_model_part,
                              problem_name,
                              output_settings)

gid_output.ExecuteInitialize()

## Sets strategies, builders, linear solvers, schemes and solving info, and fills the buffer
solver.Initialize()
solver.SetEchoLevel(echo_level)


print(" ")
print("::[PFEM Fluid Simulation]:: Analysis -START- ")

for process in list_of_processes:
    process.ExecuteBeforeSolutionLoop()

# writing a initial state results file or single file (if no restart)
if((main_model_part.ProcessInfo).Has(IS_RESTARTED)):
    if(main_model_part.ProcessInfo[IS_RESTARTED] == False):
        gid_output.ExecuteBeforeSolutionLoop()

# Set time settings
time       =  PProjectParameters["problem_data"]["start_time"].GetDouble()
final_time   = PProjectParameters["problem_data"]["end_time"].GetDouble()
Dt = PProjectParameters["problem_data"]["time_step"].GetDouble()

# Check tetrahedral mesh for wrong orientation
throw_errors = False
orientation_check = TetrahedralMeshOrientationCheck(main_model_part,throw_errors)
orientation_check.Execute()

print("here the MODEL_PART = ", main_model_part)

out = 0
step = 0

current_id = 0

while(time <= final_time):

    time = time + Dt
    step = step + 1

    main_model_part.ProcessInfo[STEP] = step
    main_model_part.CloneTimeStep(time) 

    # processes to be executed at the begining of the solution step


    print("STEP = ", step)
    print("TIME = ", time)

    if(step >= 2):
        for process in list_of_processes:
            process.ExecuteInitializeSolutionStep()

        gid_output.ExecuteInitializeSolutionStep()
     
        clock_time = StartTimeMeasuring();
        
        if(step == 2):
            process.ComputeInitialAverageMeshParameters()
        
        process.ComputeAverageMeshParameters()
 

        execute_meshing = True
        # remesh domains
        if(execute_meshing):
            print("RemeshDomains()")
            process.RemeshFluidDomains();

        solver.NodalChecksAndAssignations()

        print("Solve()")
        solver.Solve()

        print("write gid results...")
        gid_output.PrintOutput()
        #graph_printer.PrintGraphs(time)

        solver.InitializeStressStrain()
        #fluid_solver.InitializeStressStrain()

        # processes to be executed after witting the output
        for process in list_of_processes:
            process.ExecuteAfterOutputStep()

        current_id=current_id+1
        out = 0
        print(" ...gid results written")


    out = out + Dt
gid_output.ExecuteFinalize()

for process in list_of_processes:
    process.ExecuteFinalize()

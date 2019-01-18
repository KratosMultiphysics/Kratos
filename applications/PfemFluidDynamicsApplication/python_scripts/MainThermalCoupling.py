from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

#Activate it to import in the gdb path:
#import sys
#sys.path.append('/home/cpuigbo/kratos')
#x = input("stopped to allow debug: set breakpoints and press enter to continue");

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
    parallel = KratosMultiphysics.OpenMPUtils()
    parallel.SetNumThreads(int(num_threads))

def GetParallelSize():
    parallel = KratosMultiphysics.OpenMPUtils()
    return parallel.GetNumThreads()
    
#### SET NUMBER OF THREADS ####

# Import system python
import os

# Import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.SolidMechanicsApplication       as KratosSolid
import KratosMultiphysics.ExternalSolversApplication      as KratosSolvers
import KratosMultiphysics.DelaunayMeshingApplication      as KratosDelaunay
import KratosMultiphysics.PfemFluidDynamicsApplication    as KratosPfemFluid
import KratosMultiphysics.ConvectionDiffusionApplication  as KratosConvDiff

######################################################################################
######################################################################################
######################################################################################

#### PARSING THE PARAMETERS ####

# Import input
parameter_file = open("ProjectParameters.json",'r')
ProjectParameters = KratosMultiphysics.Parameters(parameter_file.read())

#set echo level
echo_level = ProjectParameters["problem_data"]["echo_level"].GetInt()

print(" ")

# defining the number of threads:
threads = ProjectParameters["problem_data"]["threads"].GetInt()
SetParallelSize(threads)
num_threads = GetParallelSize()
print("::[KPFEM Simulation]:: [OMP USING",num_threads,"THREADS ]")
#parallel.PrintOMPInfo()


print(" ")
print("::[KPFEM Simulation]:: [Time Step:", ProjectParameters["problem_data"]["time_step"].GetDouble()," echo:", echo_level,"]")

#### Model_part settings start ####

# Defining the model_part
main_model_part = KratosMultiphysics.ModelPart(ProjectParameters["problem_data"]["model_part_name"].GetString())

main_model_part.ProcessInfo.SetValue(KratosMultiphysics.SPACE_DIMENSION, ProjectParameters["problem_data"]["dimension"].GetInt())
main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, ProjectParameters["problem_data"]["dimension"].GetInt())
main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DELTA_TIME, ProjectParameters["problem_data"]["time_step"].GetDouble())
main_model_part.ProcessInfo.SetValue(KratosMultiphysics.TIME, ProjectParameters["problem_data"]["start_time"].GetDouble())
if( ProjectParameters["problem_data"].Has("gravity_vector") ):
    main_model_part.ProcessInfo.SetValue(KratosMultiphysics.GRAVITY_X, ProjectParameters["problem_data"]["gravity_vector"][0].GetDouble())
    main_model_part.ProcessInfo.SetValue(KratosMultiphysics.GRAVITY_Y, ProjectParameters["problem_data"]["gravity_vector"][1].GetDouble())
    main_model_part.ProcessInfo.SetValue(KratosMultiphysics.GRAVITY_Z, ProjectParameters["problem_data"]["gravity_vector"][2].GetDouble())

###TODO replace this "model" for real one once available in kratos core
Model = {ProjectParameters["problem_data"]["model_part_name"].GetString() : main_model_part}

#construct the solver (main setting methods are located in the solver_module)
solver_module = __import__(ProjectParameters["solver_settings"]["solver_type"].GetString())
solver = solver_module.CreateSolver(main_model_part, ProjectParameters["solver_settings"])

# Add variables (always before importing the model part)
solver.AddVariables()

# Add PfemSolidMechanicsApplication Variables
import pfem_variables  
pfem_variables.AddVariables(main_model_part) 

#thermal thing:
import eulerian_convection_diffusion_solver as convection_diffusion_solver_scripts
#we create the thermal settings.
thermal_settings = KratosMultiphysics.ConvectionDiffusionSettings()
thermal_settings.SetDiffusionVariable(KratosMultiphysics.CONDUCTIVITY)
thermal_settings.SetUnknownVariable(KratosMultiphysics.TEMPERATURE)
(main_model_part.ProcessInfo).SetValue(KratosMultiphysics.CONVECTION_DIFFUSION_SETTINGS, thermal_settings)
#adding the variables of the diffusion problem:
convection_diffusion_solver_scripts.AddVariables(main_model_part) # no settings here! (read from the model_part)


# Read model_part (note: the buffer_size is set here) (restart is read here)
solver.ImportModelPart()

# Add dofs (always after importing the model part)
if((main_model_part.ProcessInfo).Has(KratosMultiphysics.IS_RESTARTED)):
    if(main_model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED] == False):
        solver.AddDofs()
else:
    solver.AddDofs()

#adding variables for the thermal problem
convection_diffusion_solver_scripts.AddDofs(main_model_part)

#TEMP!!!!!!!! BOUNDARY CONDITIONS MUST ALSO BE SET BY THE USER!
for node in main_model_part.Nodes:
    node.SetSolutionStepValue(KratosMultiphysics.CONDUCTIVITY,0.0001)
    node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE,2.5)
    if node.Is(KratosMultiphysics.RIGID):
    	if node.X<0.298:
            node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE,0)		
            node.Fix(KratosMultiphysics.TEMPERATURE)
    	if node.X>0.298:
            node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE,5)		
            node.Fix(KratosMultiphysics.TEMPERATURE)


# Build sub_model_parts or submeshes (rearrange parts for the application of custom processes)
## Get the list of the submodel part in the object Model
for i in range(ProjectParameters["solver_settings"]["processes_sub_model_part_list"].size()):
    part_name = ProjectParameters["solver_settings"]["processes_sub_model_part_list"][i].GetString()
    if( main_model_part.HasSubModelPart(part_name) ):
        Model.update({part_name: main_model_part.GetSubModelPart(part_name)})

#### Model_part settings end ####


#print model_part and properties
if(echo_level>1):
    print("")
    print(main_model_part)
    for properties in main_model_part.Properties:
        print(properties)

#### Processes settings start ####

#obtain the list of the processes to be applied

import process_handler

process_parameters = KratosMultiphysics.Parameters("{}") 
process_parameters.AddValue("echo_level", ProjectParameters["problem_data"]["echo_level"])
process_parameters.AddValue("constraints_process_list", ProjectParameters["constraints_process_list"])
process_parameters.AddValue("loads_process_list", ProjectParameters["loads_process_list"])
if( ProjectParameters.Has("problem_process_list") ):
    process_parameters.AddValue("problem_process_list", ProjectParameters["problem_process_list"])
if( ProjectParameters.Has("output_process_list") ):
    process_parameters.AddValue("output_process_list", ProjectParameters["output_process_list"])
if( ProjectParameters.Has("processes_sub_model_part_tree_list") ):
    process_parameters.AddValue("processes_sub_model_part_tree_list",ProjectParameters["processes_sub_model_part_tree_list"])

model_processes = process_handler.ProcessHandler(Model, process_parameters)

model_processes.ExecuteInitialize()

#### processes settings end ####


# --PLOT GRAPHS OPTIONS START--###############
#problem_path = os.getcwd() #current path
#plot_active = general_variables.PlotGraphs
#graph_plot = plot_utils.GraphPlotUtility(model_part, problem_path)
# --PLOT GRAPHS OPTIONS END--#################

#### START SOLUTION ####

computing_model_part = solver.GetComputingModelPart()

## Sets strategies, builders, linear solvers, schemes and solving info, and fills the buffer
solver.Initialize()
solver.InitializeStrategy()
solver.SetEchoLevel(echo_level)


#initializing the diffusion problem
diffusion_solver = convection_diffusion_solver_scripts.EulerianConvectionDiffusionSolver(main_model_part, 2) #2d
diffusion_solver.Initialize()




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

# writing a initial state results file
current_id = 0
#if(load_restart == False):
#    if (general_variables.TryToSetTheWeight):
#        if (general_variables.TryToSetConstantWeight):
#            conditions.SetConstantWeight( general_variables.TryToSetWeightVertical, general_variables.TryToSetWeightHorizontal);
#        else:
#            conditions.SetWeight();

# set solver info starting parameters
# solving_info = solving_info_utils.SolvingInfoUtility(model_part, SolverSettings)

print(" ")
print("::[KPFEM Simulation]:: Analysis -START- ")

model_processes.ExecuteBeforeSolutionLoop()

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
while(time < end_time):

    # current time parameters
    # main_model_part.ProcessInfo.GetPreviousSolutionStepInfo()[KratosMultiphysics.DELTA_TIME] = delta_time
    delta_time = main_model_part.ProcessInfo[KratosMultiphysics.DELTA_TIME]
    
    time = time + delta_time
    step = step + 1
    
    main_model_part.ProcessInfo[KratosMultiphysics.STEP] = step
    main_model_part.CloneTimeStep(time) 

    
    print(" [STEP:",step," TIME:",time,"]")

    # processes to be executed at the begining of the solution step
    model_processes.ExecuteInitializeSolutionStep()
      
    gid_output.ExecuteInitializeSolutionStep()

    # solve time step
    clock_time = StartTimeMeasuring();

    solver.InitializeSolutionStep()

    solver.Predict()

    solver.SolveSolutionStep()

    solver.FinalizeSolutionStep()

    StopTimeMeasuring(clock_time,"Solving", False);


    #adding variables for the thermal problem
    #convection_diffusion_solver_scripts.AddDofs(main_model_part)
    for node in main_model_part.Nodes:
    	if node.IsFixed(KratosMultiphysics.TEMPERATURE)==False:
            node.AddDof(thermal_settings.GetUnknownVariable())
 
    #TEMP!!!!!!!! BOUNDARY CONDITIONS MUST ALSO BE SET BY THE USER!
    '''
    for node in main_model_part.Nodes:
    	if node.Is(KratosMultiphysics.RIGID):
            if node.X0<0.298:
                node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE,0)		
                node.Fix(KratosMultiphysics.TEMPERATURE)
            if node.X0>0.298:
                node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE,5)		
                node.Fix(KratosMultiphysics.TEMPERATURE)
    '''

    #initializing the diffusion problem.
    diffusion_solver.Initialize()

    #solving the diffusion problem.
    diffusion_solver.Solve()

    gid_output.ExecuteFinalizeSolutionStep()

    # processes to be executed at the end of the solution step
    model_processes.ExecuteFinalizeSolutionStep()

    # processes to be executed before witting the output      
    model_processes.ExecuteBeforeOutputStep()
     
    # write output results GiD: (frequency writing is controlled internally)
    if(gid_output.IsOutputStep()):
        gid_output.PrintOutput()

    # processes to be executed after witting the output
    model_processes.ExecuteAfterOutputStep()


# Ending the problem (time integration finished)
gid_output.ExecuteFinalize()

model_processes.ExecuteFinalize()

print("::[KPFEM Simulation]:: Analysis -END- ")
print(" ")

# Check solving information for any problem
#~ solver.InfoCheck() # InfoCheck not implemented yet.

#### END SOLUTION ####

# Measure process time
tfp = timer.clock()
# Measure wall time
tfw = timer.time()

print("::[KPFEM Simulation]:: [Elapsed Time = %.2f" % (tfw - t0w),"seconds] (%.2f" % (tfp - t0p),"seconds of cpu/s time)")

print(timer.ctime())

# to create a benchmark: add standard benchmark files and decomment next two lines 
# rename the file to: run_test.py
#from run_test_benchmark_results import *
#WriteBenchmarkResults(model_part)

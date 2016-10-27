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
    print("::[KPFEM Simulation]:: [OMP USING",num_threads,"THREADS ]")
    #parallel.PrintOMPInfo()
    print(" ")

#### SET NUMBER OF THREADS ####

# Import system python
import os

# Import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.SolidMechanicsApplication     as KratosSolid
import KratosMultiphysics.ExternalSolversApplication    as KratosSolvers
import KratosMultiphysics.PfemBaseApplication           as KratosPfemBase
import KratosMultiphysics.ContactMechanicsApplication   as KratosContact
import KratosMultiphysics.PfemSolidMechanicsApplication as KratosPfemSolid
import KratosMultiphysics.PfemFluidDynamicsApplication  as KratosPfemFluid

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

print(" ")
print("::[KPFEM Simulation]:: [Time Step:", ProjectParameters["problem_data"]["time_step"].GetDouble()," echo:", echo_level,"]")

#### Model_part settings start ####

# Defining the model_part
main_model_part = KratosMultiphysics.ModelPart(ProjectParameters["problem_data"]["model_part_name"].GetString())

main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, ProjectParameters["problem_data"]["domain_size"].GetInt())
main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DELTA_TIME, ProjectParameters["problem_data"]["time_step"].GetDouble())
main_model_part.ProcessInfo.SetValue(KratosMultiphysics.TIME, ProjectParameters["problem_data"]["start_time"].GetDouble())

###TODO replace this "model" for real one once available in kratos core
Model = {ProjectParameters["problem_data"]["model_part_name"].GetString() : main_model_part}

#construct the solver (main setting methods are located in the solver_module)
solver_module = __import__(ProjectParameters["solver_settings"]["solver_type"].GetString())
solver = solver_module.CreateSolver(main_model_part, ProjectParameters["solver_settings"])

# Add variables (always before importing the model part)
solver.AddVariables()

# Add PfemSolidMechanicsApplication Variables
main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL);
main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_H);

main_model_part.AddNodalSolutionStepVariable(KratosPfemBase.OFFSET);
main_model_part.AddNodalSolutionStepVariable(KratosPfemBase.SHRINK_FACTOR);
main_model_part.AddNodalSolutionStepVariable(KratosPfemBase.MEAN_ERROR);
main_model_part.AddNodalSolutionStepVariable(KratosPfemBase.RIGID_WALL);

main_model_part.AddNodalSolutionStepVariable(KratosSolid.DETERMINANT_F);

main_model_part.AddNodalSolutionStepVariable(KratosPfemSolid.WALL_TIP_RADIUS);
main_model_part.AddNodalSolutionStepVariable(KratosPfemSolid.WALL_REFERENCE_POINT);

main_model_part.AddNodalSolutionStepVariable(KratosContact.CONTACT_STRESS);


# Read model_part (note: the buffer_size is set here) (restart is read here)
solver.ImportModelPart()

# Add dofs (always after importing the model part)
if((main_model_part.ProcessInfo).Has(KratosMultiphysics.IS_RESTARTED)):
    if(main_model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED] == False):
        solver.AddDofs()
else:
    solver.AddDofs()

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


#initialize graph plot variables for time integration
#if( plot_active == True):
#  mesh_id     = 0 #general_variables.PlotMeshId
#  x_variable  = "DISPLACEMENT"
#  y_variable  = "CONTACT_FORCE"
#  graph_plot.Initialize(x_variable,y_variable,mesh_id)

#graph_write= operation_utils.TimeOperationUtility()
#graph_write_frequency = general_variables.PlotFrequency
#graph_write.InitializeTime(time, end_time, delta_time, graph_write_frequency)

#solving_print = operation_utils.TimeOperationUtility()
#solving_time_frequency = ProjectParameters["output_configuration"]["result_file_configuration"]["output_frequency"].GetDouble()
#solving_print.InitializeTime(time, end_time, delta_time, solving_time_frequency)


# Solving the problem (time integration)
while(time < end_time):

    # current time parameters
    # main_model_part.ProcessInfo.GetPreviousSolutionStepInfo()[KratosMultiphysics.DELTA_TIME] = delta_time
    delta_time = main_model_part.ProcessInfo[KratosMultiphysics.DELTA_TIME]
    
    time = time + delta_time
    step = step + 1
    
    main_model_part.ProcessInfo[KratosMultiphysics.STEP] = step
    main_model_part.CloneTimeStep(time) 

    # print process information:
    #print_info = solving_print.perform_time_operation(time)
    #if(print_info):
    #    solving_info.print_step_info(time,step)
    
    print(" [STEP:",step," TIME:",time,"]")

    # processes to be executed at the begining of the solution step
    for process in list_of_processes:
        process.ExecuteInitializeSolutionStep()
      
    gid_output.ExecuteInitializeSolutionStep()

    # solve time step
    clock_time = StartTimeMeasuring();

    solver.InitializeSolutionStep()

    solver.Predict()

    solver.SolveSolutionStep()

    solver.FinalizeSolutionStep()

    StopTimeMeasuring(clock_time,"Solving", False);

    gid_output.ExecuteFinalizeSolutionStep()

    # plot graphs
    #if(plot_active):
    #    graph_plot.SetStepResult()

    # processes to be executed at the end of the solution step
    for process in list_of_processes:
        process.ExecuteFinalizeSolutionStep()

    # processes to be executed before witting the output      
    for process in list_of_processes:
        process.ExecuteBeforeOutputStep()
     
    # write output results GiD: (frequency writing is controlled internally)
    if(gid_output.IsOutputStep()):
        gid_output.PrintOutput()

    # processes to be executed after witting the output
    for process in list_of_processes:
        process.ExecuteAfterOutputStep()

    # plot graphs
    #if(plot_active):
    #    execute_plot = graph_write.perform_time_operation(time)
    #    if(execute_plot):
    #        current_id = output_print.operation_id()
    #        graph_plot.Plot(current_id)

    # print restart file
    #if(save_restart):
    #    execute_save = restart_print.perform_time_operation(time)
    #    if(execute_save):
    #        clock_time=StartTimeMeasuring();
    #        current_id = output_print.operation_id()
    #        problem_restart.Save(time,step,current_id)
    #        solving_info.set_restart_info(execute_save,current_id)
    #        StopTimeMeasuring(clock_time,"Writing Restart", False)

    #solving_info.update_solving_info()
    #if(print_info):
    #    solving_info.print_solving_info()

# Ending the problem (time integration finished)
gid_output.ExecuteFinalize()

for process in list_of_processes:
    process.ExecuteFinalize()

print("::[KPFEM Simulation]:: Analysis -END- ")
print(" ")

# Check solving information for any problem
#~ solver.InfoCheck() # InfoCheck not implemented yet.

# check solving information for any problem
# solving_info.info_check()

#### END SOLUTION ####

# Measure process time
tfp = timer.clock()
# Measure wall time
tfw = timer.time()

print("::[KPFEM Simulation]:: [Elapsed Time = %.2f" % (tfp - t0p),"seconds] (%.2f" % (tfw - t0w),"seconds of cpu/s time)")

print(timer.ctime())

# to create a benchmark: add standard benchmark files and decomment next two lines 
# rename the file to: run_test.py
#from run_test_benchmark_results import *
#WriteBenchmarkResults(model_part)

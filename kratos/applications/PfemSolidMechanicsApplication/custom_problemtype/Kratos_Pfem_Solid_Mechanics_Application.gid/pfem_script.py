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
import KratosMultiphysics.PfemSolidMechanicsApplication as KratosPfemSolid


#import the python utilities:
import restart_utility                 as restart_utils
import pfem_gid_output_utility         as gid_utils

import pfem_conditions_python_utility  as condition_utils
import list_files_python_utility       as files_utils

import rigid_wall_python_utility       as wall_utils
import graph_plot_python_utility       as plot_utils

import time_operation_utility          as operation_utils
import solving_info_utility            as solving_info_utils

######################################################################################
######################################################################################
######################################################################################

#### PARSING THE PARAMETERS ####

# Import input
parameter_file = open("ProjectParameters.json",'r')
ProjectParameters = KratosMultiphysics.Parameters(parameter_file.read())

#set echo level
echo_level = ProjectParameters["problem_data"]["echo_level"].GetInt()

#### Model_part settings start ####

# Defining the model_part
model_part = KratosMultiphysics.ModelPart(ProjectParameters["problem_data"]["model_part_name"].GetString())

model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, ProjectParameters["problem_data"]["domain_size"].GetInt())
#model_part.ProcessInfo.SetValue(KratosMultiphysics.DELTA_TIME, ProjectParameters["problem_data"]["time_step"].GetDouble())
#model_part.ProcessInfo.SetValue(KratosMultiphysics.TIME, ProjectParameters["problem_data"]["start_time"].GetDouble())

###TODO replace this "model" for real one once available in kratos core
Model = {ProjectParameters["problem_data"]["model_part_name"].GetString() : model_part}

#construct the solver (main setting methods are located in the solver_module)
#solver_module = __import__(ProjectParameters["solver_settings"]["solver_type"].GetString())
#solver = solver_module.CreateSolver(model_part, ProjectParameters["solver_settings"])


#### PARSING CLASSICAL PARAMETERS ####

# Import input
import ProjectParameters as general_variables

# defining the type, the name and the path of the problem:
echo_level   = general_variables.EchoLevel
domain_size  = general_variables.domain_size

problem_type = general_variables.ProblemType
problem_name = general_variables.problem_name

#problem_path = general_variables.problem_path #fixed path
problem_path = os.getcwd() #current path

# defining a model part
#model_part = ModelPart("Solid Domain")

# defining solver settings
SolverSettings = general_variables.SolverSettings

# defining the model size to scale
length_scale = 1.0

# --RIGID WALL OPTIONS START--################
# set rigid wall contact if it is active:
# activated instead of classical contact
# set rigid wall configuration
rigid_wall = wall_utils.RigidWallUtility(model_part, domain_size, general_variables.rigid_wall_config)

# --RIGID WALL OPTIONS END--##################


# --SET MESH MODELER START--##################
#construct meshing domains
#meshing_domains = []
#domains_list = ProjectParameters["meshing_domains"]
#for i in range(0,domains_list.size()):
#    item = domains_list[i]
#    domain_module = __import__(item["python_file_name"].GetString())
#    domain = domain_module.CreateMeshingDomain(model_part,item)
#    meshing_domains.append(domain)
#remesh_domains = general_variables.RemeshDomains
#contact_search = general_variables.FindContacts
#rigid_wall_contact_search = general_variables.FindRigidWallContacts
#modeler = modeler_utils.ModelerUtility(model_part, domain_size, remesh_domains, contact_search, rigid_wall_contact_search)

# print check
#if(echo_level>1):
#    print("::[KPFEM Simulation]:: MESH DOMAINS :", len(general_variables.MeshConditions))
#    for conditions in general_variables.MeshConditions:
#        print("::[KPFEM Simulation]:: --> Domain [", conditions["Subdomain"], "] ", conditions["MeshElement"])
# --SET MESH MODELER END--####################


# --READ AND SET MODEL FILES--###############

# set the restart of the problem
restart_step = general_variables.Restart_Step
problem_restart = restart_utils.RestartUtility(model_part, problem_path, problem_name)

# set the results file list of the problem (managed by the problem_restart and gid_print)
print_lists = general_variables.PrintLists
output_mode = general_variables.GidOutputConfiguration.GiDPostMode
list_files = files_utils.ListFilesUtility(problem_path, problem_name, print_lists, output_mode);
list_files.Initialize(general_variables.file_list);

# --READ AND SET MODEL FILES END--############


# --DEFINE CONDITIONS START--#################
incr_disp = general_variables.Incremental_Displacement
incr_load = general_variables.Incremental_Load
rotation_dofs = SolverSettings.RotationDofs
conditions    = condition_utils.ConditionsUtility(model_part, domain_size, incr_disp, incr_load, rotation_dofs);

# --DEFINE CONDITIONS END--###################


# --GID OUTPUT OPTIONS START--###############
# set gid print options
gid_print = gid_utils.GidOutputUtility(problem_name, general_variables.GidOutputConfiguration)

# --GID OUTPUT OPTIONS END--##################

# --PLOT GRAPHS OPTIONS START--###############

plot_active = general_variables.PlotGraphs
graph_plot = plot_utils.GraphPlotUtility(model_part, problem_path)

# --PLOT GRAPHS OPTIONS END--#################

# start problem initialization:
print(" ")

# defining the number of threads:
num_threads = general_variables.NumberofThreads
SetParallelSize(num_threads)

# --DEFINE MAIN SOLVER START--################

print(" ")
print("::[KPFEM Simulation]:: [Time Step:", general_variables.time_step," echo:", echo_level,"]")

# import solver file
solver_constructor = __import__(SolverSettings.solver_type)

# construct the solver
main_step_solver = solver_constructor.CreateSolver(model_part, SolverSettings)

# --DEFINE MAIN SOLVER END--##################


# --CONFIGURATIONS END--######################
# ----------------------------------------------------------------#


# --START SOLUTION--######################
#
#initialize problem : load restart or initial start
load_restart     = general_variables.LoadRestart
save_restart     = general_variables.SaveRestart

#set buffer size
buffer_size = 3;

#define problem variables:
solver_constructor.AddVariables( model_part, SolverSettings)

# set PfemSolidApplicationVariables
model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL);

model_part.AddNodalSolutionStepVariable(KratosPfemBase.OFFSET);
model_part.AddNodalSolutionStepVariable(KratosPfemBase.SHRINK_FACTOR);

model_part.AddNodalSolutionStepVariable(KratosPfemBase.MEAN_ERROR);
model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_H);

model_part.AddNodalSolutionStepVariable(KratosSolid.DETERMINANT_F);

# if hasattr(SolverSettings, "RigidWalls"):
    # if SolverSettings.RigidWalls == True:
model_part.AddNodalSolutionStepVariable(KratosPfemBase.RIGID_WALL);
model_part.AddNodalSolutionStepVariable(KratosPfemSolid.WALL_TIP_RADIUS);
model_part.AddNodalSolutionStepVariable(KratosPfemSolid.WALL_REFERENCE_POINT);



#--- READ MODEL ------#
if(load_restart == False):
  
  print("::[KPFEM Simulation]:: Reading -START- (MDPA FILE) ")

  #remove results, restart, graph and list previous files
  problem_restart.CleanPreviousFiles()
  list_files.RemoveListFiles()

  #reading the model
  model_part_io = KratosMultiphysics.ModelPartIO(problem_name)
  model_part_io.ReadModelPart(model_part)
 
  #set the buffer size
  model_part.SetBufferSize(buffer_size)
  #Note: the buffer size should be set once the mesh is read for the first time

  print("::[KPFEM Simulation]:: Reading -END- ")

  model_part.ProcessInfo[KratosMultiphysics.LOAD_RESTART] = 0
  model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED] = False

  #set the degrees of freedom
  solver_constructor.AddDofs(model_part, SolverSettings)

  #set the constitutive law
  import constitutive_law_python_utility as constitutive_law_utils

  constitutive_law = constitutive_law_utils.ConstitutiveLawUtility(model_part,domain_size)
  constitutive_law.Initialize()

else:

  print("::[KPFEM Simulation]:: Reading -RESTART- [FILE",restart_step,"]")

  #reading the model from the restart file
  problem_restart.Load(restart_step);

  print("::[KPFEM Simulation]:: Reading -END- ")

  model_part.ProcessInfo[KratosMultiphysics.LOAD_RESTART] = 1
  model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED] = True

  #remove results, restart, graph and list posterior files
  problem_restart.CleanPosteriorFiles(restart_step)
  list_files.ReBuildListFiles()


# --RIGID WALL OPTIONS START--################
# set rigid wall contact if it is active:
# activated instead of classical contact
# set rigid wall configuration
rigid_wall = wall_utils.RigidWallUtility(model_part, domain_size, general_variables.rigid_wall_config)

# --RIGID WALL OPTIONS END--##################


#### Processes settings start ####
#obtain the list of the processes to be applied

import process_factory
#the process order of execution is important
#list_of_processes  = process_factory.KratosProcessFactory(Model).ConstructListOfProcesses( ProjectParameters["constraints_process_list"] )
#list_of_processes += process_factory.KratosProcessFactory(Model).ConstructListOfProcesses( ProjectParameters["loads_process_list"] )
if(ProjectParameters.Has("problem_process_list")):
    list_of_processes = process_factory.KratosProcessFactory(Model).ConstructListOfProcesses( ProjectParameters["problem_process_list"] )
#if(ProjectParameters.Has("output_process_list")):
#    list_of_processes += process_factory.KratosProcessFactory(Model).ConstructListOfProcesses( ProjectParameters["output_process_list"] )
            
#print list of constructed processes
if(echo_level>1):
    for process in list_of_processes:
        print(process)

for process in list_of_processes:
    process.ExecuteInitialize()

#### processes settings end ####


# --BUILD MESH MODELER START--####################
# build mesh modeler
#for domain in meshing_domains:
#    domain.SetImposedWalls(rigid_wall)
#    domain.Initialize()
#modeler.BuildMeshModelers(meshing_domains)

# set rigid walls
#if(rigid_wall_contact_search):
#    modeler.SetRigidWall(rigid_wall)
# --BUILD MESH MODELER END--####################

# --CONTACT SEARCH START--####################
# build mesh modeler
# modeler.BuildContactModeler(general_variables.contact_modeler_config)
# --CONTACT SEARCH END--######################


#--- MODELER INITIALIZATION---#

#set mesh searches and modeler
#modeler.InitializeDomains( load_restart ); ## due to the skin conditions at reloading

# mesh size nodal h search
#if(load_restart == False):
#    modeler.SearchNodalH();

#--- PRINT CONTROL ---#

if(echo_level>=1):
    print("")
    print(model_part)
    print(model_part.Properties[1])

# --INITIALIZE--###########################
#

# set delta time in process info
model_part.ProcessInfo[KratosMultiphysics.DELTA_TIME] = general_variables.time_step

# solver initialize
main_step_solver.Initialize()
main_step_solver.SetRestart(load_restart) #calls strategy initialize if no restart

# initial contact search
# modeler.InitialContactSearch()

#define time steps and loop range of steps
delta_time = model_part.ProcessInfo[KratosMultiphysics.DELTA_TIME]

if(load_restart):  

  buffer_size  = 0

else:

  model_part.ProcessInfo[KratosMultiphysics.TIME]                = 0
  model_part.ProcessInfo[KratosMultiphysics.TIME_STEPS]          = 0
  model_part.ProcessInfo[KratosSolid.PREVIOUS_DELTA_TIME] = delta_time

  conditions.Initialize(delta_time);

#initialize step operations
step     = model_part.ProcessInfo[KratosMultiphysics.STEP]
time     = model_part.ProcessInfo[KratosMultiphysics.TIME]
end_time = general_variables.nsteps * delta_time


output_print = operation_utils.TimeOperationUtility()
gid_time_frequency = general_variables.GiDWriteFrequency
output_print.InitializeTime(time, end_time, delta_time, gid_time_frequency)

restart_print = operation_utils.TimeOperationUtility()
restart_time_frequency = general_variables.RestartFrequency
restart_print.InitializeTime(time, end_time, delta_time, restart_time_frequency)

#mesh_generation = operation_utils.TimeOperationUtility()
#mesh_generation_frequency = modeler.GetRemeshFrequency()
#mesh_generation.InitializeTime(time, end_time, delta_time, mesh_generation_frequency)

contact_search = operation_utils.TimeOperationUtility()
contact_search_frequency = general_variables.contact_modeler_config.contact_search_frequency
contact_search.InitializeTime(time, end_time, delta_time, contact_search_frequency)

rigid_wall_contact_search = operation_utils.TimeOperationUtility()
rigid_wall_contact_search_frequency = 0
rigid_wall_contact_search.InitializeTime(time, end_time, delta_time, rigid_wall_contact_search_frequency)

#initialize graph plot variables for time integration
if( plot_active == True):
  mesh_id     = 0 #general_variables.PlotMeshId
  x_variable  = "DISPLACEMENT"
  y_variable  = "CONTACT_FORCE"
  graph_plot.Initialize(x_variable,y_variable,mesh_id)

graph_write= operation_utils.TimeOperationUtility()
graph_write_frequency = general_variables.PlotFrequency
graph_write.InitializeTime(time, end_time, delta_time, graph_write_frequency)

solving_print = operation_utils.TimeOperationUtility()
solving_time_frequency = gid_time_frequency
solving_print.InitializeTime(time, end_time, delta_time, solving_time_frequency)

# --TIME INTEGRATION--#######################
#
  
#writing a single file
gid_print.initialize_results(model_part)

# filling the buffer
for step in range(0,buffer_size):

  model_part.CloneTimeStep(time)
  model_part.ProcessInfo[KratosMultiphysics.DELTA_TIME] = delta_time
  model_part.ProcessInfo[KratosMultiphysics.STEP] = step-buffer_size


# Set time settings
step       = model_part.ProcessInfo[KratosMultiphysics.STEP]
time       = model_part.ProcessInfo[KratosMultiphysics.TIME]

end_time   = ProjectParameters["problem_data"]["end_time"].GetDouble()
delta_time = ProjectParameters["problem_data"]["time_step"].GetDouble()


# writing a initial state results file
current_id = 0
if(load_restart == False):
    if (general_variables.TryToSetTheWeight):
        if (general_variables.TryToSetConstantWeight):
            conditions.SetConstantWeight( general_variables.TryToSetWeightVertical, general_variables.TryToSetWeightHorizontal);
        else:
            conditions.SetWeight();
    gid_print.write_results(model_part, general_variables.nodal_results, general_variables.gauss_points_results, time, step, current_id)
    list_files.PrintListFiles(current_id);

# set solver info starting parameters
solving_info = solving_info_utils.SolvingInfoUtility(model_part, SolverSettings)


print(" ")
print("::[KPFEM Simulation]:: Analysis -START- ")

# Solving the problem (time integration)
while(time <= end_time):

    # current time parameters
    # model_part.ProcessInfo.GetPreviousStepInfo(1)[KratosMultiphysics.DELTA_TIME] = delta_time
    delta_time = model_part.ProcessInfo[KratosMultiphysics.DELTA_TIME]
    
    time = time + delta_time
    step = step + 1
    
    model_part.ProcessInfo[KratosMultiphysics.STEP] = step
    model_part.CloneTimeStep(time) 


    # print process information:
    print_info = solving_print.perform_time_operation(time)
    if(print_info):
        solving_info.print_step_info(time,step)
        
    # processes to be executed at the begining of the solution step
    execute_rigid_wall_contact_search = rigid_wall_contact_search.perform_time_operation(time)
    if(execute_rigid_wall_contact_search):
        rigid_wall.ExecuteContactSearch()

    for process in list_of_processes:
        process.ExecuteInitializeSolutionStep()
      
    # solving the solid problem 
    clock_time = StartTimeMeasuring();

    # solve time step
    main_step_solver.Solve()

    StopTimeMeasuring(clock_time,"Solving", False);

    rigid_wall.UpdatePosition()

    #plot graphs
    if(plot_active):
        graph_plot.SetStepResult()

    #incremental load
    conditions.SetIncrementalLoad(step, delta_time);
    ##conditions.CorrectBoundaryConditions(step, delta_time); ## function to remove load conditions from the contact...

    # processes to be executed at the end of the solution step
    for process in list_of_processes:
        process.ExecuteFinalizeSolutionStep()

    # processes to be executed before witting the output      
    for process in list_of_processes:
        process.ExecuteBeforeOutputStep()
     
    #print the results at the end of the step
    execute_write = output_print.perform_time_operation(time)
    if(execute_write):
        clock_time=StartTimeMeasuring();
        current_id = output_print.operation_id()
        #print gid output file
        gid_print.write_results(model_part,general_variables.nodal_results,general_variables.gauss_points_results,time,step,current_id)
        #print on list files
        list_files.PrintListFiles(current_id);
        solving_info.set_print_info(execute_write, current_id)
        StopTimeMeasuring(clock_time,"Writing Results", False);
        

    # processes to be executed after witting the output
    for process in list_of_processes:
        process.ExecuteAfterOutputStep()

    # plot graphs
    if(plot_active):
        execute_plot = graph_write.perform_time_operation(time)
        if(execute_plot):
            current_id = output_print.operation_id()
            graph_plot.Plot(current_id)

    # print restart file
    if(save_restart):
        execute_save = restart_print.perform_time_operation(time)
        if(execute_save):
            clock_time=StartTimeMeasuring();
            current_id = output_print.operation_id()
            problem_restart.Save(time,step,current_id)
            solving_info.set_restart_info(execute_save,current_id)
            StopTimeMeasuring(clock_time,"Writing Restart", False)


    solving_info.update_solving_info()
    if(print_info):
        solving_info.print_solving_info()

    conditions.RestartImposedDisp()

# --FINALIZE--############################
#

# writing a single file
gid_print.finalize_results()

for process in list_of_processes:
    process.ExecuteFinalize()

print("::[KPFEM Simulation]:: Analysis -END- ")
print(" ")

# --END--###############################
#

# check solving information for any problem
solving_info.info_check()

#### END SOLUTION ####

# Measure process time
tfp = timer.clock()
# Measure wall time
tfw = timer.time()

print("::[KSM Simulation]:: [ Computing Time = (%.2f" % (tfp - t0p)," seconds process time) ( %.2f" % (tfw - t0w)," seconds wall time) ]")

print(timer.ctime())

# to create a benchmark: add standard benchmark files and decomment next two lines 
# rename the file to: run_test.py
#from run_test_benchmark_results import *
#WriteBenchmarkResults(model_part)

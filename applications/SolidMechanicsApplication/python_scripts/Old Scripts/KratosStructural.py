#Activate it to import in the gdb path:
#import sys
#sys.path.append('/home/jmaria/kratos')
#x = raw_input("stopped to allow debug: set breakpoints and press enter to continue");

##################################################################
### ***************GENERAL MAIN OF THE ANALISYS****************###
##################################################################

##time control starts
from time import *
print ctime()
# measure process time
t0p = clock()
# measure wall time
#t0w = time()

#----------------------------------------------------------------#
######################--CONFIGURATIONS START--####################
#Import the general variables read from the GiD
import ProjectParameters as general_variables

#setting the domain size for the problem to be solved
domain_size = general_variables.domain_size

#including kratos path
from KratosMultiphysics import *

#including Applications paths
from KratosMultiphysics.ExternalSolversApplication import *
from KratosMultiphysics.SolidMechanicsApplication import *
  
#import the python solver:
import solid_mechanics_main_solver as main_solver

#import the python utilities:
import restart_python_utility       as restart_utils
import print_results_python_utility as gid_utils

import conditions_python_utility    as condition_utils
import list_files_python_utility    as files_utils

#import modeler_python_utility       as modeler_utils
#import graph_plot_python_utility    as plot_utils



#------------------------#--FUNCTIONS START--#------------------#
#---------------------------------------------------------------#

######################--PRINTING TIME START --###################
def StartTimeMeasuring():
  # measure process time
  time_ip = clock()
  return time_ip

def StopTimeMeasuring(time_ip,process):
  # measure process time
  time_fp = clock()
  print  process," [ spent time = ",time_fp - time_ip,"] "
######################--PRINTING TIME END --#####################

######################--SET NUMBER OF THREADS --#################
def SetParallelSize(num_threads):
  parallel=OpenMPUtils()
  print "Num Threads = ", num_threads
  parallel.SetNumThreads(int(num_threads)); 
######################--SET NUMBER OF THREADS --#################

#------------------------#--FUNCTIONS END--#--------------------#
#---------------------------------------------------------------#


#defining the number of threads:
num_threads = general_variables.NumberofThreads
SetParallelSize(num_threads)

#defining the type, the name and the path of the problem:
problem_type = general_variables.ProblemType
problem_name = general_variables.problem_name
problem_path = general_variables.problem_path

#defining a model part
model_part = ModelPart("SolidDomain");

#defining the model size to scale
length_scale = 1.0


######################--READ AND SET MODEL FILES--###############

#set the restart of the problem
restart_step     = general_variables.Restart_Step
problem_restart  = restart_utils.RestartUtility(model_part,problem_path,problem_name,restart_step);

#set the results file list of the problem (managed by the problem_restart and gid_print)
print_lists      = general_variables.PrintLists
list_files       = files_utils.ListFilesUtility(problem_path,problem_name,print_lists);
list_files.Initialize(general_variables.file_list);

#initialize problem : load restart or initial start
load_restart     = general_variables.LoadRestart
save_restart     = general_variables.SaveRestart
restart_interval = general_variables.Restart_Interval
rotation_dofs    = general_variables.Rotational_Dofs
problem_restart.Initialize(load_restart,save_restart,restart_interval,main_solver,list_files,rotation_dofs);
  

######################--READ AND SET MODEL FILES END--############


######################--DEFINE MAIN SOLVER START--################

#set time integration solver
echo_level = int(general_variables.echo_level)
main_step_solver = main_solver.SolidMechanicsSolver(model_part,domain_size,echo_level) 

#defining the problem type
solver_type  = general_variables.SolverType
line_search  = general_variables.LineSearch
main_step_solver.SetProblemType(problem_type,solver_type,line_search);

#defining the linear solver
linear_solver_type = general_variables.LinearSolver
solver_tolerance   = general_variables.Linear_Solver_Tolerance
max_iters          = general_variables.Linear_Solver_Max_Iteration
main_step_solver.SetLinearSolver(linear_solver_type,solver_tolerance,max_iters);

#defining the convergence criterion
criterion_type  = general_variables.Convergence_Criteria
convergence_tol = general_variables.Convergence_Tolerance
absolute_tol    = general_variables.Absolute_Tolerance
max_iters       = int(general_variables.Max_Iter)
main_step_solver.SetConvergenceCriterion(criterion_type,convergence_tol,absolute_tol,max_iters,rotation_dofs);

######################--DEFINE MAIN SOLVER END--##################


######################--DEFINE CONDITIONS START--#################

incr_disp  = general_variables.Incremental_Displacement
incr_load  = general_variables.Incremental_Load
conditions = condition_utils.ConditionsUtility(model_part,domain_size,incr_disp,incr_load);

######################--DEFINE CONDITIONS END--###################



######################--GID OUTPUT OPTIONS START--###############

gid_configuration_mode   =  WriteDeformedMeshFlag.WriteDeformed
gid_variables            =  WriteConditionsFlag.WriteElementsOnly
gid_output_mode          =  GiDPostMode.GiD_PostBinary
gid_files_mode           =  MultiFileFlag.MultipleFiles

if(general_variables.GiDWriteMeshFlag == "False"):
  gid_configuration_mode = WriteDeformedMeshFlag.WriteUndeformed
if(general_variables.GiDWriteConditionsFlag == "True"):
  gid_mesh_write_type    = WriteConditionsFlag.WriteConditions
if(general_variables.GiDPostMode == "Ascii"):
  gid_output_mode = GiDPostMode.GiD_PostAscii
if(general_variables.GiDMultiFileFlag == "Single"):
  gid_files_mode = MultiFileFlag.SingleFile

#Force to Multiple files write if it is not a StaticSolver
if(solver_type != "StaticSolver" and general_variables.GiDMultiFileFlag == "Single"):
  gid_files_mode = MultiFileFlag.MultipleFiles

gid_print = gid_utils.PrintResultsUtility(model_part,problem_type,solver_type,problem_name,gid_output_mode,gid_files_mode)

#set gid print options
write_particles   =  general_variables.GiDWriteParticlesFlag
write_deformed    =  general_variables.GiDWriteMeshFlag
write_conditions  =  general_variables.GiDWriteConditionsFlag
write_frequency   =  general_variables.WriteFrequency

gid_print.SetPrintOptions(write_particles,write_deformed,write_conditions,write_frequency)

######################--GID OUTPUT OPTIONS END--##################


######################--PLOT GRAPHS OPTIONS START--###############

#plot_active    = general_variables.PlotGraphs
#plot_frequency = general_variables.PlotFrequency

#graph_plot  = plot_utils.GraphPlotUtility(model_part,problem_path,plot_active,plot_frequency);

#x_var   = "TIME"
#y_var   = "REACTION"
#mesh_id = 1

#plot variables on the domain which is remeshed
#for conditions in general_variables.MeshConditions:
#  if(conditions["Remesh"] == 1):
#    mesh_id =int(conditions["Subdomain"])

#print " Graph Subdomain ", mesh_id
#graph_plot.SetPlotVariables(x_var,y_var,mesh_id);

######################--PLOT GRAPHS OPTIONS END--#################

######################--CONFIGURATIONS END--######################
#----------------------------------------------------------------#



#########################--START SOLUTION--######################
#################################################################


#read model_part / set buffer_size / set dofs / set constitutive_law
buffer_size   =  3;

#--- READ MODEL ------#
#problem_restart.StartModelRead(buffer_size,domain_size,problem_type,rotation_dofs,main_solver,modeler);
problem_restart.StartModelRead(buffer_size,domain_size,problem_type,rotation_dofs,main_solver);

#--- PRINT CONTROL ---#
print model_part
print model_part.Properties


#########################--INITIALIZE--###########################
##################################################################

# solver initialize
load_restart = general_variables.LoadRestart 
main_step_solver.Initialize(load_restart)

# initial contact search
#modeler.InitialContactSearch()


#initialize time integration variables
current_step = 0

#define time and time_step
time_step    =  general_variables.time_step
steps_number =  general_variables.nsteps

#define loop range of steps
istep        = 0
nstep        = int(general_variables.nsteps) + buffer_size 

start_steps  = buffer_size-1
start_time   = start_steps*time_step;

problem_restart.SetTimeVariables(time_step,steps_number)

#problem_restart.InitializeTimeIntegration(gid_print,modeler,graph_plot,conditions);
problem_restart.InitializeTimeIntegration(gid_print,conditions);

#redefine loop range of steps after problem_restart
istep        = problem_restart.step_i
nstep        = problem_restart.step_n
time_step    = problem_restart.TimeStep()
current_step = problem_restart.current_step

#######################--TIME INTEGRATION--#######################
##################################################################

  
#writing the initial mesh
#gid_print.PrintInitialMesh()


for step in range(istep,nstep):
  time = time_step*step
  model_part.CloneTimeStep(time)
  model_part.ProcessInfo[DELTA_TIME] = time_step
  model_part.ProcessInfo[TIME_STEPS] = step

  if( (time-start_time) > 0 ):
    print "STEP = ", step - start_steps - 1
    print "TIME = ", time - start_time
    model_part.ProcessInfo[TIME] = time - start_time
  
  step_printed = False;

  #solving the solid problem
  if(step > start_steps ):
    
    clock_time=StartTimeMeasuring();
    #solve time step non-linear system
    main_step_solver.Solve()
    StopTimeMeasuring(clock_time,"Solving");

    #plot graphs

    #update previous time step
    model_part.ProcessInfo[PREVIOUS_DELTA_TIME] = time_step;
    
    #incremental load
    incr_steps = step-start_steps
    conditions.SetIncrementalLoad(incr_steps,time_step);
      
    #print the results at the end of the step
    if(general_variables.WriteResults == "PreMeshing"):
      clock_time=StartTimeMeasuring();
      step_printed = gid_print.PrintResults(time,current_step,list_files)
      StopTimeMeasuring(clock_time,"Write Results");
       
    #print restart file
    if( step_printed == True ):
      write_id = model_part.ProcessInfo[WRITE_ID];
      #graph_plot.Plot(write_id)

      clock_time=StartTimeMeasuring();
      problem_restart.PrintRestartFile(write_id);
      StopTimeMeasuring(clock_time,"Restart");

    current_step = current_step + 1
    
##########################--FINALIZE--############################
##################################################################

print "Analysis Finalized "

############################--END--###############################
##################################################################

# measure process time
tfp = clock()
# measure wall time
#tfw = time()

print ctime()
#print "Analysis Completed  [Process Time = ", tfp - t0p, "seconds, Wall Time = ", tfw - t0w, " ]"
print "Analysis Completed  [Process Time = ",tfp - t0p,"] "

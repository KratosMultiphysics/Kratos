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

#import the python utilities:
import restart_utility              as restart_utils
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

######################--DEFINE MAIN SOLVER START--################

SolverSettings = general_variables.SolverSettings

#import solver file
solver_constructor = __import__(SolverSettings.solver_type)

#construct the solver
main_step_solver   = solver_constructor.CreateSolver( model_part, SolverSettings) 

######################--DEFINE MAIN SOLVER END--##################


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

problem_restart.Initialize(load_restart,save_restart,restart_interval,list_files);
  
######################--READ AND SET MODEL FILES END--############



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

if(SolverSettings.scheme_type != "StaticSolver" and general_variables.GiDMultiFileFlag == "Single"):
  gid_files_mode = MultiFileFlag.MultipleFiles

gid_print = gid_utils.PrintResultsUtility(model_part,problem_type,SolverSettings.solver_type,problem_name,gid_output_mode,gid_files_mode)

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

#set buffer size
buffer_size = 3;

#define problem variables:
solver_constructor.AddVariables( model_part, SolverSettings)


#--- READ MODEL ------#
if(load_restart == "False"):
  #reading the model 
  model_part_io = ModelPartIO(problem_name)
  model_part_io.ReadModelPart(model_part)
  
  #set the buffer size
  model_part.SetBufferSize(buffer_size)
  #Note: the buffer size should be set once the mesh is read for the first time

  #set the degrees of freedom
  solver_constructor.AddDofs(model_part, SolverSettings)

  #set the constitutive law
  import constitutive_law_python_utility as constitutive_law_utils

  constitutive_law = constitutive_law_utils.ConstitutiveLawUtility(model_part,domain_size);
  constitutive_law.Initialize();


#set mesh searches and modeler
# modeler.InitializeDomains();

#if(load_restart == "False"):
  #find nodal h
  #modeler.SearchNodalH();

#set writing numeration
problem_restart.SetStartSteps(buffer_size);

#--- PRINT CONTROL ---#
print model_part
print model_part.Properties


#########################--INITIALIZE--###########################
##################################################################

# solver initialize
main_step_solver.Initialize()
main_step_solver.SetRestart(load_restart)

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

problem_restart.InitializeTimeIntegration(time_step,steps_number);

#redefine loop range of steps after problem_restart
istep        = problem_restart.step_i
nstep        = problem_restart.step_n
time_step    = problem_restart.TimeStep()
current_step = problem_restart.current_step

#initialize print results variables for time integration
write_id = model_part.ProcessInfo[WRITE_ID];
gid_print.Initialize(start_time,nstep,current_step,write_id)

if(load_restart == "False"):
  conditions.Initialize();

#initialize mesh modeling variables for time integration
#modeler.Initialize(current_step,current_step)

#initialize graph plot variables for time integration
#graph_plot.Initialize(current_step)


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

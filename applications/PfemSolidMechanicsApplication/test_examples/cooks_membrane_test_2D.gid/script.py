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
from KratosMultiphysics.MeshingApplication        import *
from KratosMultiphysics.PfemSolidMechanicsApplication import *

#import the python utilities:
import restart_utility              as restart_utils
import gid_output_utility           as gid_utils

import conditions_python_utility    as condition_utils
import list_files_python_utility    as files_utils

import time_operation_utility       as operation_utils
import modeler_python_utility       as modeler_utils
import rigid_wall_python_utility    as wall_utils
#import graph_plot_python_utility    as plot_utils



#------------------------#--FUNCTIONS START--#------------------#
#---------------------------------------------------------------#

######################--TIME MONITORING START--##################
def StartTimeMeasuring():
  # measure process time
  time_ip = clock()
  return time_ip

def StopTimeMeasuring(time_ip,process):
  # measure process time
  time_fp = clock()
  print  process," [ spent time = ",time_fp - time_ip,"] "
######################--TIME MONITORING END --###################

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


######################--SET MESH MODELER START--##################

remesh_domains = general_variables.RemeshDomains
contact_search = general_variables.FindContacts
modeler = modeler_utils.ModelerUtility(model_part,domain_size,remesh_domains, contact_search);

# Optional : mesh refinement based on tool characteristics 

#(deffault arch=5-10 degrees)
#critical_radius      = 0.00004
#critical_radius      = 0.025
#critical_radius      = general_variables.tip_radius

critical_radius      = general_variables.mesh_modeler_config.critical_radius

if(critical_radius > 5*general_variables.rigid_wall_config.tip_radius):
  critical_radius    = general_variables.rigid_wall_config.tip_radius

# Optional : mesh refinement b#defining the mesh conditions

# print check
print " MESH CONDITIONS :", len(general_variables.MeshConditions)
for conditions in general_variables.MeshConditions:
    print " --> Domain [", conditions["Subdomain"],"] ",  conditions["MeshElement"]

#build mesh modeler
modeler.BuildMeshModeler(general_variables.mesh_modeler_config);
  

######################--CONTACT SEARCH START--####################

#build mesh modeler
modeler.BuildContactModeler(general_variables.contact_modeler_config);

######################--CONTACT SEARCH END--######################


######################--SET MESH MODELER END--####################



######################--RIGID WALL OPTIONS START--################

#set rigid wall contact if it is active: 
#activated instead of classical contact

#set rigid wall configuration
rigid_wall = wall_utils.RigidWallUtility(model_part,domain_size);

rigid_wall.Initialize(general_variables.rigid_wall_config);

######################--RIGID WALL OPTIONS END--##################


######################--READ AND SET MODEL FILES--###############

#set the restart of the problem
restart_step     = general_variables.Restart_Step
problem_restart  = restart_utils.RestartUtility(model_part,problem_path,problem_name);

#set the results file list of the problem (managed by the problem_restart and gid_print)
print_lists      = general_variables.PrintLists
list_files       = files_utils.ListFilesUtility(problem_path,problem_name,print_lists);
list_files.Initialize(general_variables.file_list);
  
######################--READ AND SET MODEL FILES END--############



######################--DEFINE CONDITIONS START--#################

incr_disp  = general_variables.Incremental_Displacement
incr_load  = general_variables.Incremental_Load
conditions = condition_utils.ConditionsUtility(model_part,domain_size,incr_disp,incr_load);

######################--DEFINE CONDITIONS END--###################



######################--GID OUTPUT OPTIONS START--###############

#set gid print options
gid_print = gid_utils.GidOutputUtility(problem_name, general_variables.GidOutputConfiguration)

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

#initialize problem : load restart or initial start
load_restart     = general_variables.LoadRestart
save_restart     = general_variables.SaveRestart

#set buffer size
buffer_size = 3;

#define problem variables:
solver_constructor.AddVariables( model_part, SolverSettings)


#--- READ MODEL ------#
if(load_restart == "False"):

  problem_restart.CleanPreviousFiles(list_files)

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

else:
  
  #reading the model from the restart file
  problem_restart.Load(restart_time);

  problem_restart.CleanPosteriorFiles(time_step,restart_time,list_files)

#set mesh searches and modeler
# modeler.InitializeDomains();

#if(load_restart == "False"):
  #find nodal h
  #modeler.SearchNodalH();


#--- PRINT CONTROL ---#
print model_part
print model_part.Properties[1]


#########################--INITIALIZE--###########################
##################################################################

# solver initialize
main_step_solver.Initialize()
main_step_solver.SetRestart(load_restart)

# initial contact search
#modeler.InitialContactSearch()

#define time steps and loop range of steps
if(load_restart == "True"):  

  istep        = model_part.ProcessInfo[TIME_STEPS]+1
  nstep        = int(general_variables.nsteps) + buffer_size 
  time_step    = model_part.ProcessInfo[DELTA_TIME]
  current_step = istep-nstep

else:

  istep        = 0
  nstep        = int(general_variables.nsteps) + buffer_size 
  time_step    = general_variables.time_step
  current_step = 0

  model_part.ProcessInfo[PREVIOUS_DELTA_TIME] = time_step;

  conditions.Initialize();


#initialize step operations
starting_time  = current_step * time_step
ending_time    = general_variables.nsteps * time_step;
 
output_print = operation_utils.TimeOperationUtility()
gid_time_frequency = general_variables.GiDWriteFrequency
output_print.InitializeTime(starting_time,ending_time,time_step,gid_time_frequency)

restart_print = operation_utils.TimeOperationUtility()
restart_time_frequency = general_variables.RestartFrequency
restart_print.InitializeTime(starting_time,ending_time,time_step,restart_time_frequency)


#initialize mesh modeling variables for time integration
#modeler.Initialize(current_step,current_step)

#initialize graph plot variables for time integration
#graph_plot.Initialize(current_step)


#######################--TIME INTEGRATION--#######################
##################################################################
  
#writing a single file
gid_print.initialize_results(model_part)

#initialize time integration variables
start_steps  = buffer_size-1
start_time   = start_steps*time_step;

for step in range(istep,nstep):

  time = time_step*step
  model_part.CloneTimeStep(time)
  model_part.ProcessInfo[DELTA_TIME] = time_step
  model_part.ProcessInfo[TIME_STEPS] = step

  current_time = time - start_time
  current_step = step - start_steps - 1

  if( current_time > 0 ):
    print "STEP = ", current_step
    print "TIME = ", current_time
    model_part.ProcessInfo[TIME] = current_time
  
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
    incr_steps = current_step + 1
    conditions.SetIncrementalLoad(incr_steps,time_step);
      
    #print the results at the end of the step
    if(general_variables.WriteResults == "PreMeshing"):
      execute_write = output_print.perform_time_operation(current_time)
      if( execute_write == True ):
        clock_time=StartTimeMeasuring();
        #print gid output file
        gid_print.write_results(model_part,general_variables.nodal_results,general_variables.gauss_points_results,current_time,current_step)
        #print on list files
        list_files.PrintListFiles(current_step);
        StopTimeMeasuring(clock_time,"Write Results");
        #plot graphs
        #graph_plot.Plot(current_time)

    #print restart file
    if( save_restart == True ):
      execute_save = restart_print.perform_time_operation(current_time)
      if( execute_save == True ):
        clock_time=StartTimeMeasuring();
        problem_restart.Save(current_time,current_step);
        StopTimeMeasuring(clock_time,"Restart");

    
##########################--FINALIZE--############################
##################################################################

#writing a single file
gid_print.finalize_results()

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

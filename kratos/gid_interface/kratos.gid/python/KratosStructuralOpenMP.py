from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# Activate it to import in the gdb path:
# import sys
# sys.path.append('/home/jmaria/kratos')
# x = raw_input("stopped to allow debug: set breakpoints and press enter to continue");

#
# ***************GENERAL MAIN OF THE ANALISYS****************###
#

# time control starts
from time import *
print(ctime())
# measure process time
t0p = clock()
# measure wall time
t0w = time()

# ----------------------------------------------------------------#
# --CONFIGURATIONS START--####################
# Import the general variables read from the GiD
import ProjectParameters as general_variables

# setting the domain size for the problem to be solved
domain_size = general_variables.domain_size

# including kratos path
from KratosMultiphysics import *

# including Applications paths
from KratosMultiphysics.SolidMechanicsApplication  import *
from KratosMultiphysics.ExternalSolversApplication import *


# import the python utilities:
import restart_utility           as restart_utils
import gid_output_utility        as gid_utils

import conditions_python_utility as condition_utils
import list_files_python_utility as files_utils

import time_operation_utility    as operation_utils
import solving_info_utility      as solving_info_utils

# ------------------------#--FUNCTIONS START--#------------------#
# ---------------------------------------------------------------#
# --TIME MONITORING START--##################
def StartTimeMeasuring():
    # measure process time
    time_ip = clock()
    return time_ip

def StopTimeMeasuring(time_ip, process, report):
    # measure process time
    time_fp = clock()
    if( report ):
        used_time = time_fp - time_ip
        print("::[KSM Simulation]:: [ %.2f" % round(used_time,2),"s", process," ] ")
# --TIME MONITORING END --###################

# --SET NUMBER OF THREADS --#################

def SetParallelSize(num_threads):
    parallel = OpenMPUtils()
    parallel.SetNumThreads(int(num_threads))
    print("::[KSM Simulation]:: OMP USING",num_threads,"THREADS ]")
    #parallel.PrintOMPInfo()
    print(" ")

# --SET NUMBER OF THREADS --#################

# ------------------------#--FUNCTIONS END--#--------------------#
# ---------------------------------------------------------------#

# defining the type, the name and the path of the problem:
echo_level   = general_variables.echo_level
problem_type = general_variables.ProblemType
problem_name = general_variables.problem_name

#problem_path = general_variables.problem_path #fixed path
problem_path  = os.getcwd() #current path

# defining a model part
model_part = ModelPart("Solid Domain")

# defining solver settings
SolverSettings = general_variables.SolverSettings

# defining the model size to scale
length_scale = 1.0

# --READ AND SET MODEL FILES--###############

# set the restart of the problem
restart_step = general_variables.Restart_Step
problem_restart = restart_utils.RestartUtility(model_part, problem_path, problem_name)

# set the results file list of the problem (managed by the problem_restart and gid_print)
print_lists = general_variables.PrintLists
output_mode = general_variables.GidOutputConfiguration.GiDPostMode
list_files = files_utils.ListFilesUtility(problem_path, problem_name, print_lists, output_mode)
list_files.Initialize(general_variables.file_list)

# --READ AND SET MODEL FILES END--############


# --DEFINE CONDITIONS START--#################
incr_disp = general_variables.Incremental_Displacement
incr_load = general_variables.Incremental_Load
rotation_dofs = SolverSettings.RotationDofs
conditions    = condition_utils.ConditionsUtility(model_part, domain_size, incr_disp, incr_load, rotation_dofs)

# --DEFINE CONDITIONS END--###################


# --GID OUTPUT OPTIONS START--###############
# set gid print options
gid_print = gid_utils.GidOutputUtility(problem_name, general_variables.GidOutputConfiguration)

# --GID OUTPUT OPTIONS END--##################

# start problem initialization:
print(" ")

# defining the number of threads:
num_threads = general_variables.NumberofThreads
SetParallelSize(num_threads)

# --DEFINE MAIN SOLVER START--################

# import solver file
solver_constructor = __import__(SolverSettings.solver_type)

# construct the solver
main_step_solver = solver_constructor.CreateSolver(model_part, SolverSettings)

# --DEFINE MAIN SOLVER END--##################


# --CONFIGURATIONS END--######################
# ----------------------------------------------------------------#


# --START SOLUTION--######################
#
# initialize problem : load restart or initial start
load_restart = general_variables.LoadRestart
save_restart = general_variables.SaveRestart

# set buffer size
buffer_size = 3

# define problem variables:
solver_constructor.AddVariables(model_part, SolverSettings)


# --- READ MODEL ------#
if(load_restart == False):

    print("::[KSM Simulation]:: Reading -START- (MDPA FILE) ")

    # remove results, restart, graph and list previous files
    problem_restart.CleanPreviousFiles()
    list_files.RemoveListFiles()

    # reading the model
    model_part_io = ModelPartIO(problem_name)
    model_part_io.ReadModelPart(model_part)

    # set the buffer size
    model_part.SetBufferSize(buffer_size)
    # Note: the buffer size should be set once the mesh is read for the first time

    print("::[KSM Simulation]:: Reading -END- ")

    # set the degrees of freedom
    solver_constructor.AddDofs(model_part, SolverSettings)

    # set the constitutive law
    import constitutive_law_python_utility as constitutive_law_utils

    constitutive_law = constitutive_law_utils.ConstitutiveLawUtility(model_part, domain_size);
    constitutive_law.Initialize();


else:

    print("::[KSM Simulation]:: Reading -RESTART- (RESTART FILE ",restart_step,")- ")

    # reading the model from the restart file
    problem_restart.Load(restart_step);

    print("::[KSM Simulation]:: Reading -END- ")

    # remove results, restart, graph and list posterior files
    problem_restart.CleanPosteriorFiles(restart_step)
    list_files.ReBuildListFiles()


# set mesh searches and modeler
# modeler.InitializeDomains();

# if(load_restart == False):
    # find nodal h
    # modeler.SearchNodalH();


# --- PRINT CONTROL ---#

if(echo_level>1):
    print("")
    print(model_part)
    print(model_part.Properties[1])


# --INITIALIZE--###########################
#

# set delta time in process info
model_part.ProcessInfo[DELTA_TIME] = general_variables.time_step

# solver initialize
main_step_solver.Initialize()
main_step_solver.SetRestart(load_restart) #calls strategy initialize if no restart

# initial contact search
# modeler.InitialContactSearch()

#define time steps and loop range of steps
time_step = model_part.ProcessInfo[DELTA_TIME]

# define time steps and loop range of steps
if(load_restart):

    buffer_size = 0

else:

    model_part.ProcessInfo[TIME]                = 0
    model_part.ProcessInfo[TIME_STEPS]          = 0
    model_part.ProcessInfo[PREVIOUS_DELTA_TIME] = time_step;

    conditions.Initialize(time_step);


# initialize step operations
starting_step  = model_part.ProcessInfo[TIME_STEPS]
starting_time  = model_part.ProcessInfo[TIME]
#ending_step    = general_variables.nsteps
ending_time    = general_variables.end_time


output_print = operation_utils.TimeOperationUtility()
gid_time_frequency = general_variables.GiDWriteFrequency
output_print.InitializeTime(starting_time, ending_time, time_step, gid_time_frequency)

restart_print = operation_utils.TimeOperationUtility()
restart_time_frequency = general_variables.RestartFrequency
restart_print.InitializeTime(starting_time, ending_time, time_step, restart_time_frequency)

solving_print = operation_utils.TimeOperationUtility()
solving_time_frequency = gid_time_frequency
solving_print.InitializeTime(starting_time, ending_time, time_step, solving_time_frequency)


# --TIME INTEGRATION--#######################
#

# writing a single file
gid_print.initialize_results(model_part)

#initialize time integration variables
current_time = starting_time
current_step = starting_step

# filling the buffer
for step in range(0,buffer_size):

  model_part.CloneTimeStep(current_time)
  model_part.ProcessInfo[DELTA_TIME] = time_step
  model_part.ProcessInfo[TIME_STEPS] = step-buffer_size

# writing a initial state results file
current_id = 0
if(load_restart == False): 
    gid_print.write_results(model_part, general_variables.nodal_results, general_variables.gauss_points_results, current_time, current_step, current_id)
    list_files.PrintListFiles(current_id);

# set solver info starting parameters
solving_info = solving_info_utils.SolvingInfoUtility(model_part, SolverSettings)

print(" ")
print("::[KSM Simulation]:: Analysis -START- ")

# solving the problem
while(current_time < ending_time):

  # store previous time step
  model_part.ProcessInfo[PREVIOUS_DELTA_TIME] = time_step
  # set new time step ( it can change when solve is called )
  time_step = model_part.ProcessInfo[DELTA_TIME]

  current_time = current_time + time_step
  current_step = current_step + 1
  model_part.CloneTimeStep(current_time)
  model_part.ProcessInfo[TIME] = current_time

  # print process information:
  print_info = solving_print.perform_time_operation(current_time)
  if(print_info):
      solving_info.print_step_info(current_time,current_step)

  # ssolving the solid problem 
  clock_time = StartTimeMeasuring();     

  # solve time step non-linear system
  main_step_solver.Solve()
  
  StopTimeMeasuring(clock_time, "Solving", False);

  # incremental load
  conditions.SetIncrementalLoad(current_step, time_step);

  # print the results at the end of the step
  execute_write = output_print.perform_time_operation(current_time)
  if(execute_write):
      clock_time = StartTimeMeasuring();
      current_id = output_print.operation_id()
      # print gid output file
      gid_print.write_results(model_part, general_variables.nodal_results, general_variables.gauss_points_results, current_time, current_step, current_id)
      # print on list files
      list_files.PrintListFiles(current_id)
      solving_info.set_print_info(execute_write, current_id)    
      StopTimeMeasuring(clock_time, "Writing Results", False);

  # print restart file
  if(save_restart):
      execute_save = restart_print.perform_time_operation(current_time)
      if(execute_save):
          clock_time = StartTimeMeasuring();
          current_id = restart_print.operation_id()
          problem_restart.Save(current_time, current_step, current_id)
          solving_info.set_restart_info(execute_save,current_id)
          StopTimeMeasuring(clock_time, "Writing Restart", False);


  solving_info.update_solving_info()
  if(print_info):
      solving_info.print_solving_info()
        
  conditions.RestartImposedDisp()

# --FINALIZE--############################
#

# writing a single file
gid_print.finalize_results()

print("::[KSM Simulation]:: Analysis -END- ")
print(" ")

# --END--###############################
#

# check solving information for any problem
solving_info.info_check()

# measure process time
tfp = clock()
# measure wall time
tfw = time()

print("::[KSM Simulation]:: [ Computing Time = (%.2f" % (tfp - t0p)," seconds process time) ( %.2f" % (tfw - t0w)," seconds wall time) ]")

print(ctime())

# to create a benchmark: add standard benchmark files and decomment next two lines 
# rename the file to: run_test.py
#from run_test_benchmark_results import *
#WriteBenchmarkResults(model_part)

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

# ----------------------------------------------------------------#
# --CONFIGURATIONS START--####################
# Import the general variables read from the GiD
import ProjectParameters as general_variables

# setting the domain size for the problem to be solved
domain_size = general_variables.domain_size

# including kratos path
from KratosMultiphysics import *

# including Applications paths
from KratosMultiphysics.ExternalSolversApplication import *
from KratosMultiphysics.SolidMechanicsApplication import *

# import the python utilities:
# import restart_utility as restart_utils
# import gid_output_utility as gid_utils

import conditions_python_utility as condition_utils
# import list_files_python_utility as files_utils

import time_operation_utility as operation_utils
# import modeler_python_utility       as modeler_utils
# import graph_plot_python_utility    as plot_utils


# ------------------------#--FUNCTIONS START--#------------------#
# ---------------------------------------------------------------#
# --TIME MONITORING START--##################
def StartTimeMeasuring():
    # measure process time
    time_ip = clock()
    return time_ip


def StopTimeMeasuring(time_ip, process):
    # measure process time
    time_fp = clock()
    print(" ", process, " [ spent time = ", time_fp - time_ip, "] ")
# --TIME MONITORING END --###################

# --SET NUMBER OF THREADS --#################


def SetParallelSize(num_threads):
    parallel = OpenMPUtils()
    print("Num Threads = ", num_threads)
    parallel.SetNumThreads(int(num_threads))
# --SET NUMBER OF THREADS --#################

# ------------------------#--FUNCTIONS END--#--------------------#
# ---------------------------------------------------------------#


# defining the number of threads:
num_threads = general_variables.NumberofThreads
SetParallelSize(num_threads)

# defining the type, the name and the path of the problem:
problem_type = general_variables.ProblemType
problem_name = general_variables.problem_name
problem_path = general_variables.problem_path

# defining a model part
model_part = ModelPart("SolidDomain")

# defining the model size to scale
length_scale = 1.0

# --DEFINE MAIN SOLVER START--################

SolverSettings = general_variables.SolverSettings

# import solver file
solver_constructor = __import__(SolverSettings.solver_type)

# construct the solver
main_step_solver = solver_constructor.CreateSolver(model_part, SolverSettings)

# --DEFINE MAIN SOLVER END--##################


# --READ AND SET MODEL FILES--###############

# set the restart of the problem
restart_step = general_variables.Restart_Step
# problem_restart = restart_utils.RestartUtility(model_part, problem_path, problem_name)

# set the results file list of the problem (managed by the problem_restart and gid_print)
print_lists = general_variables.PrintLists
output_mode = general_variables.GidOutputConfiguration.GiDPostMode
# list_files = files_utils.ListFilesUtility(problem_path, problem_name, print_lists, output_mode)
# list_files.Initialize(general_variables.file_list)

# --READ AND SET MODEL FILES END--############


# --DEFINE CONDITIONS START--#################
incr_disp = general_variables.Incremental_Displacement
incr_load = general_variables.Incremental_Load
conditions = condition_utils.ConditionsUtility(model_part, domain_size, incr_disp, incr_load)

# --DEFINE CONDITIONS END--###################

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
# reading the model
model_part_io = ModelPartIO(problem_name)
model_part_io.ReadModelPart(model_part)

# set the buffer size
model_part.SetBufferSize(buffer_size)
# Note: the buffer size should be set once the mesh is read for the first time

# set the degrees of freedom
solver_constructor.AddDofs(model_part, SolverSettings)

# set the constitutive law
import constitutive_law_python_utility as constitutive_law_utils

constitutive_law = constitutive_law_utils.ConstitutiveLawUtility(model_part, domain_size);
constitutive_law.Initialize();

# --- PRINT CONTROL ---#
print(model_part)
print(model_part.Properties[1])

# --INITIALIZE--###########################
#

# solver initialize
main_step_solver.Initialize()
main_step_solver.SetRestart(load_restart)

# initial contact search
# modeler.InitialContactSearch()

# define time steps and loop range of steps
istep = 0
nstep = int(general_variables.nsteps) + buffer_size
time_step = general_variables.time_step
current_step = 0
model_part.ProcessInfo[PREVIOUS_DELTA_TIME] = time_step;
conditions.Initialize(time_step);

# initialize step operations
starting_time = current_step * time_step
ending_time = general_variables.nsteps * time_step;

output_print = operation_utils.TimeOperationUtility()
gid_time_frequency = general_variables.GiDWriteFrequency
output_print.InitializeTime(starting_time, ending_time, time_step, gid_time_frequency)

# --TIME INTEGRATION--#######################
#

# initialize time integration variables
start_steps = buffer_size - 1
start_time = start_steps * time_step;

for step in range(istep, nstep):

    time = time_step * step
    model_part.CloneTimeStep(time)
    model_part.ProcessInfo[DELTA_TIME] = time_step
    model_part.ProcessInfo[TIME_STEPS] = step

    current_time = time - start_time
    current_step = step - start_steps - 1

    if(current_time > 0):
        print("STEP = ", current_step)
        print("TIME = ", current_time)
        model_part.ProcessInfo[TIME] = current_time

    step_printed = False;

    # solving the solid problem
    if(step > start_steps):

        clock_time = StartTimeMeasuring();
        # solve time step non-linear system
        main_step_solver.Solve()
        StopTimeMeasuring(clock_time, "Solving");

        # plot graphs

        # update previous time step
        model_part.ProcessInfo[PREVIOUS_DELTA_TIME] = time_step;

        # incremental load
        incr_steps = current_step + 1
        conditions.SetIncrementalLoad(incr_steps, time_step);
        
        conditions.RestartImposedDisp()

# --FINALIZE--############################
#

print("Analysis Finalized ")

# --END--###############################
#

# measure process time
tfp = clock()
print(ctime())
print("Analysis Completed  [Process Time = ", tfp - t0p, "] ")

# benchmarking...
import sys
kratos_benchmarking_path = '../../../../benchmarking'  # kratos_root/benchmarking
sys.path.append(kratos_benchmarking_path)
import benchmarking
if (benchmarking.InBenchmarkingMode()):
    # find max Y rotation
    r_max = 0.0
    for node in model_part.Nodes:
        ir = node.GetSolutionStepValue(ROTATION_Z)
        if(ir > r_max):
            r_max = ir
    # write
    abs_tol = 1e-9
    rel_tol = 1e-5
    benchmarking.Output(r_max, "Z Rotation", abs_tol, rel_tol)

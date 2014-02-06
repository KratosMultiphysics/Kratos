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
# t0w = time()

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
from KratosMultiphysics.MeshingApplication import *
from KratosMultiphysics.PfemSolidMechanicsApplication import *

# import the python utilities:
import restart_utility as restart_utils
import gid_output_utility as gid_utils

import conditions_python_utility as condition_utils
import list_files_python_utility as files_utils

import time_operation_utility as operation_utils
import modeler_python_utility as modeler_utils
import rigid_wall_python_utility as wall_utils
import graph_plot_python_utility as plot_utils


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


# --RIGID WALL OPTIONS START--################
# set rigid wall contact if it is active:
# activated instead of classical contact
# set rigid wall configuration
rigid_wall = wall_utils.RigidWallUtility(model_part, domain_size, general_variables.rigid_wall_config)

# --RIGID WALL OPTIONS END--##################


# --SET MESH MODELER START--##################

remesh_domains = general_variables.RemeshDomains
contact_search = general_variables.FindContacts
modeler = modeler_utils.ModelerUtility(model_part, domain_size, remesh_domains, contact_search)

# print check
print(" MESH CONDITIONS :", len(general_variables.MeshConditions))
for conditions in general_variables.MeshConditions:
    print(" --> Domain [", conditions["Subdomain"], "] ", conditions["MeshElement"])

# build mesh modeler
modeler.BuildMeshModeler(general_variables.mesh_modeler_config)
modeler.SetRigidWall(rigid_wall)

# --CONTACT SEARCH START--####################

# build mesh modeler
modeler.BuildContactModeler(general_variables.contact_modeler_config)

# --CONTACT SEARCH END--######################


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
conditions = condition_utils.ConditionsUtility(model_part, domain_size, incr_disp, incr_load);

# --DEFINE CONDITIONS END--###################


# --GID OUTPUT OPTIONS START--###############
# set gid print options
gid_print = gid_utils.GidOutputUtility(problem_name, general_variables.GidOutputConfiguration)

# --GID OUTPUT OPTIONS END--##################


# --PLOT GRAPHS OPTIONS START--###############

plot_active = general_variables.PlotGraphs
graph_plot = plot_utils.GraphPlotUtility(model_part, problem_path)

# --PLOT GRAPHS OPTIONS END--#################

# --CONFIGURATIONS END--######################
# ----------------------------------------------------------------#


# --START SOLUTION--######################
#
# initialize problem : load restart or initial start
load_restart = general_variables.LoadRestart
save_restart = general_variables.SaveRestart

# set buffer size
buffer_size = 3;

# define problem variables:
solver_constructor.AddVariables(model_part, SolverSettings)

# set PfemSolidApplicationVariables
model_part.AddNodalSolutionStepVariable(NORMAL);

model_part.AddNodalSolutionStepVariable(OFFSET);
model_part.AddNodalSolutionStepVariable(SHRINK_FACTOR);

model_part.AddNodalSolutionStepVariable(MEAN_ERROR);
model_part.AddNodalSolutionStepVariable(NODAL_H);

# if hasattr(SolverSettings, "RigidWalls"):
    # if SolverSettings.RigidWalls == True:
model_part.AddNodalSolutionStepVariable(RIGID_WALL);
model_part.AddNodalSolutionStepVariable(WALL_TIP_RADIUS);
model_part.AddNodalSolutionStepVariable(WALL_REFERENCE_POINT);

model_part.AddNodalSolutionStepVariable(CONTACT_FORCE);
model_part.AddNodalSolutionStepVariable(DETERMINANT_F);

# --- READ MODEL ------#
if(load_restart == False):

    # remove results, restart, graph and list previous files
    problem_restart.CleanPreviousFiles()
    list_files.RemoveListFiles()

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

else:

    # reading the model from the restart file
    problem_restart.Load(restart_step);

    # remove results, restart, graph and list posterior files
    problem_restart.CleanPosteriorFiles(restart_step)
    list_files.ReBuildListFiles()

# set mesh searches and modeler
print("initialize modeler")
modeler.InitializeDomains();

if(load_restart == False):
    # find nodal h
    print("search mesh nodal_h")
    modeler.SearchNodalH();


# --- PRINT CONTROL ---#
print(model_part)
print(model_part.Properties[1])


# --INITIALIZE--###########################
#

# solver initialize
main_step_solver.Initialize()
main_step_solver.SetRestart(load_restart)

# initial contact search
modeler.InitialContactSearch()

# define time steps and loop range of steps
if(load_restart):

    istep = model_part.ProcessInfo[TIME_STEPS] + 1
    nstep = int(general_variables.nsteps)
    time_step = model_part.ProcessInfo[DELTA_TIME]
    current_step = istep
    buffer_size = 0

else:

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

restart_print = operation_utils.TimeOperationUtility()
restart_time_frequency = general_variables.RestartFrequency
restart_print.InitializeTime(starting_time, ending_time, time_step, restart_time_frequency)

contact_search = operation_utils.TimeOperationUtility()
contact_search_frequency = general_variables.contact_modeler_config.contact_search_frequency
contact_search.InitializeTime(starting_time, ending_time, time_step, contact_search_frequency)

rigid_wall_contact_search = operation_utils.TimeOperationUtility()
rigid_wall_contact_search_frequency = 0
rigid_wall_contact_search.InitializeTime(starting_time, ending_time, time_step, rigid_wall_contact_search_frequency)

# initialize mesh modeling variables for time integration
modeler.Initialize(current_step, current_step)

# initialize graph plot variables for time integration
if(plot_active):
    mesh_id = 0  # general_variables.PlotMeshId
    x_variable = "DISPLACEMENT"
    y_variable = "REACTION"
    graph_plot.Initialize(x_variable, y_variable, mesh_id)

graph_write = operation_utils.TimeOperationUtility()
graph_write_frequency = general_variables.PlotFrequency
graph_write.InitializeTime(starting_time, ending_time, time_step, graph_write_frequency)


# --TIME INTEGRATION--#######################
#

# writing a single file
gid_print.initialize_results(model_part)

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

        # processes to be executed at the begining of the solution step
        execute_rigid_wall_contact_search = rigid_wall_contact_search.perform_time_operation(current_time)
        if(execute_rigid_wall_contact_search):
            rigid_wall.ExecuteContactSearch()

        # solve time step non-linear system
        main_step_solver.Solve()
        StopTimeMeasuring(clock_time, "Solving");

        # processes to be executed at the end of the solution step
        rigid_wall.UpdatePosition()

        # plot graphs
        if(plot_active):
            graph_plot.SetStepResult()

        # update previous time step
        model_part.ProcessInfo[PREVIOUS_DELTA_TIME] = time_step;

        # incremental load
        incr_steps = current_step + 1
        conditions.SetIncrementalLoad(incr_steps, time_step);

        # print the results at the end of the step
        if(general_variables.WriteResults == "PreMeshing"):
            execute_write = output_print.perform_time_operation(current_time)
            if(execute_write):
                clock_time = StartTimeMeasuring();
                current_id = output_print.operation_id()
                # print gid output file
                gid_print.write_results(model_part, general_variables.nodal_results, general_variables.gauss_points_results, current_time, current_step, current_id)
                # print on list files
                list_files.PrintListFiles(current_id);
                StopTimeMeasuring(clock_time, "Write Results");

        # remesh domains
        modeler.RemeshDomains(current_step);

        # contact search
        modeler.ContactSearch(current_step);

        # print the results at the end of the step
        if(general_variables.WriteResults == "PostMeshing"):
            execute_write = output_print.perform_time_operation(current_time)
            if(execute_write):
                clock_time = StartTimeMeasuring();
                current_id = output_print.operation_id()
                # print gid output file
                gid_print.write_results(model_part, general_variables.nodal_results, general_variables.gauss_points_results, current_time, current_step, current_id)
                # print on list files
                list_files.PrintListFiles(current_id);
                StopTimeMeasuring(clock_time, "Write Results");

           # plot graphs
        if(plot_active):
            execute_plot = graph_write.perform_time_operation(current_time)
            if(execute_plot):
                current_id = output_print.operation_id()
                graph_plot.Plot(current_id)

        # print restart file
        if(save_restart):
            execute_save = restart_print.perform_time_operation(current_time)
            if(execute_save):
                clock_time = StartTimeMeasuring();
                current_id = output_print.operation_id()
                problem_restart.Save(current_time, current_step, current_id);
                StopTimeMeasuring(clock_time, "Restart");

        conditions.RestartImposedDisp()

# --FINALIZE--############################
#

# writing a single file
gid_print.finalize_results()

print("Analysis Finalized ")

# --END--###############################
#

# measure process time
tfp = clock()
# measure wall time
# tfw = time()

print(ctime())
# print "Analysis Completed  [Process Time = ", tfp - t0p, "seconds, Wall Time = ", tfw - t0w, " ]"
print("Analysis Completed  [Process Time = ", tfp - t0p, "] ")

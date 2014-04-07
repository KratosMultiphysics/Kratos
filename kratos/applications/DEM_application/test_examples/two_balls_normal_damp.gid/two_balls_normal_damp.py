from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import time as timer
import os
import sys
import math


# including kratos path
kratos_path = '../../../..'  # BENCHMARK ##########
kratos_libs_path = '../../../../libs'  # kratos_root/libs
kratos_applications_path = '../../../../applications'  # kratos_root/applications
kratos_benchmarking_path = '../../../../benchmarking'  # BENCHMARK ##########

import sys
sys.path.append(kratos_path)
sys.path.append(kratos_libs_path)
sys.path.append(kratos_applications_path)
sys.path.append(kratos_benchmarking_path)  # BENCHMARK ##########

from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *

from DEM_explicit_solver_var import *
import DEM_explicit_solver_var as Param
import DEM_procedures
proc = DEM_procedures.Procedures(Param)

import benchmarking  # BENCHMARK ##########


def FindNode(node_list, X, Y, Z):
    for node in node_list:
        if((node.X - X) ** 2 + (node.Y - Y) ** 2 + (node.Z -Z)**2 < .000001):
            return node


def BenchmarkCheck(time, node1):  # BENCHMARK ##########
    benchmarking.Output(time, "Time",1e-7)  # BENCHMARK ##########
    benchmarking.Output(node1.GetSolutionStepValue(DISPLACEMENT_Y), "Node Displacement", 0.01)  # BENCHMARK ##########
    benchmarking.Output(node1.GetSolutionStepValue(VELOCITY_Y), "Node Velocity", 0.01)  # BENCHMARK ##########

#---------------------MODEL PART KRATOS AND GID.IO ------------------------------------------------------------------

# Defining a model part for the solid part

my_timer = Timer();
balls_model_part = ModelPart("SolidPart");

# Importing the strategy object

if (Param.ElementType == "SphericParticle3D" or Param.ElementType == "CylinderParticle2D"):
    import sphere_strategy as SolverStrategy
elif (Param.ElementType == "SphericContinuumParticle3D"):
    import continuum_sphere_strategy as SolverStrategy

SolverStrategy.AddVariables(balls_model_part, Param)

# Reading the model_part: binary or ascii, multifile or single --> only binary and single for mpi.

if (Param.OutputFileType == "Binary"):
    gid_mode = GiDPostMode.GiD_PostBinary

else:
    gid_mode = GiDPostMode.GiD_PostAscii

if (Param.Multifile == "multiple_files"):
    multifile = MultiFileFlag.MultipleFiles

else:
    multifile = MultiFileFlag.SingleFile

deformed_mesh_flag = WriteDeformedMeshFlag.WriteDeformed
write_conditions = WriteConditionsFlag.WriteConditions

gid_io = GidIO(problem_name, gid_mode, multifile, deformed_mesh_flag, write_conditions)
model_part_io_solid = ModelPartIO(problem_name)
model_part_io_solid.ReadModelPart(balls_model_part)

# Setting up the buffer size: SHOULD BE DONE AFTER READING!!!

balls_model_part.SetBufferSize(2)

# Adding dofs

SolverStrategy.AddDofs(balls_model_part)

# Creating a solver object

solver = SolverStrategy.ExplicitStrategy(balls_model_part, Param)  # here, solver variables initialize as default


# Creating necessary directories

main_path = os.getcwd()
post_path = str(main_path) + '/' + str(problem_name) + '_Post_Files'
list_path = str(main_path) + '/' + str(problem_name) + '_Post_Lists'
neigh_list_path = str(main_path) + '/' + str(problem_name) + '_Neigh_Lists'
data_and_results = str(main_path) + '/' + str(problem_name) + '_Results_and_Data'
graphs_path = str(main_path) + '/' + str(problem_name) + '_Graphs'
MPI_results = str(main_path) + '/' + str(problem_name) + '_MPI_results'

for directory in [post_path, list_path, neigh_list_path, data_and_results, graphs_path, MPI_results]:

    if not os.path.isdir(directory):
        os.makedirs(str(directory))

os.chdir(list_path)

multifile = open(problem_name + '_all' + '.post.lst', 'w')
multifile_5 = open(problem_name + '_5' + '.post.lst', 'w')
multifile_10 = open(problem_name + '_10' + '.post.lst', 'w')
multifile_50 = open(problem_name + '_50' + '.post.lst', 'w')

multifile.write('Multiple\n')
multifile_5.write('Multiple\n')
multifile_10.write('Multiple\n')
multifile_50.write('Multiple\n')

first_print = True
index_5 = 1
index_10 = 1
index_50 = 1
prev_time = 0.0
control = 0.0


os.chdir(main_path)

export_model_part = balls_model_part

#-------------------------------------------------------------------------------------------------------------------------

#------------------------------------------DEM_PROCEDURES FUNCTIONS & INITIALITZATIONS--------------------------------------------------------

Pressure = 0.0
Pressure = proc.GiDSolverTransfer(balls_model_part, solver)

if (Param.ModelDataInfo == "ON"):
    os.chdir(data_and_results)
    proc.ModelData(balls_model_part, solver)       # calculates the mean number of neighbours the mean radius, etc..
    os.chdir(main_path)

print('Initializing Problem....')

solver.Initialize()

print('Initialitzation Complete' + '\n')

#
node1 = FindNode(balls_model_part.Nodes, 0.0, 1.0, 0.0)  # BENCHMARK ##########
print(node1)  # there is a memory problem with the string        ########## BENCHMARK ##########
#

dt = balls_model_part.ProcessInfo.GetValue(DELTA_TIME)

step = 0
time = 0.0
time_old_print = 0.0
initial_pr_time = timer.clock()
initial_real_time = timer.time()

#-------------------------------------------------------------------------------------------------------------------------------------

#-----------------------SINGLE FILE MESH AND RESULTS INITIALITZATION-------------------------------------------------------------------

os.chdir(post_path)

if (Param.Multifile == "single_file"):

    if (Param.ContactMeshOption == "ON"):
        gid_io.InitializeMesh(0.0)
        gid_io.WriteMesh(contact_model_part.GetMesh());
        gid_io.FinalizeMesh()
        gid_io.InitializeResults(0.0, contact_model_part.GetMesh());

    gid_io.InitializeMesh(0.0)
    gid_io.WriteSphereMesh(balls_model_part.GetMesh())
    gid_io.FinalizeMesh()
    gid_io.InitializeResults(0.0, balls_model_part.GetMesh());

#------------------------------------------------------------------------------------------

#
#
# MAIN LOOP                                            #
#
#
os.chdir(main_path)

print(('Main loop starts at instant: ' + str(initial_pr_time) + '\n'))

total_steps_expected = int(FinalTime / dt)

print(('Total number of TIME STEPs expected in the calculation are: ' + str(total_steps_expected) + ' if time step is kept ' + '\n'))

while (time < FinalTime):

    dt = balls_model_part.ProcessInfo.GetValue(DELTA_TIME)  # Possible modifications of DELTA_TIME
    time = time + dt
    balls_model_part.CloneTimeStep(time)
    balls_model_part.ProcessInfo[TIME_STEPS] = step

    # _SOLVE_#########################################4
    os.chdir(main_path)
    solver.Solve()
    # TIME CONTROL######################################4

    incremental_time = (timer.time() - initial_real_time) - prev_time

    if (incremental_time > ControlTime):
        percentage = 100 * (float(step) / total_steps_expected)

        print('Real time calculation: ' + str(timer.time() - initial_real_time))
        print('Percentage Completed: ' + str(percentage) + ' %')
        print("TIME STEP = " + str(step) + '\n')

        prev_time = (timer.time() - initial_real_time)

    if ((timer.time() - initial_real_time > 60) and first_print):
        first_print = False
        estimated_sim_duration = 60 * (total_steps_expected / step)  # seconds

        print(('The calculation total estimated time is ' + str(estimated_sim_duration) + 'seconds' + '\n'))
        print(('in minutes:' + str(estimated_sim_duration / 60) + 'min.' + '\n'))
        print(('in hours:' + str(estimated_sim_duration / 3600) + 'hrs.' + '\n'))
        print(('in days:' + str(estimated_sim_duration / 86400) + 'days' + '\n'))

        if (estimated_sim_duration / 86400 > 2.0):
            print(('WARNING!!!:       VERY LASTING CALCULATION' + '\n'))

    # CONCRETE_TEST_STUFF#########################################4

    os.chdir(data_and_results)

    total_force = 0.0
    force_node = 0.0

    if (FixVelocitiesOption == 'OFF'):
        TimePercentageFixVelocities = 0.0

    os.chdir(list_path)
    multifile.write(problem_name + '_' + str(time) + '.post.bin\n')
    os.chdir(main_path)

    # ___GiD IO____#########################################4

    time_to_print = time - time_old_print

    if (time_to_print >= Param.OutputTimeStep):
        BenchmarkCheck(time, node1)  # BENCHMARK ##########
        os.chdir(data_and_results)

        if (index_5 == 5):
            multifile_5.write(problem_name + '_' + str(time) + '.post.bin\n')
            index_5 = 0

        if (index_10 == 10):
            multifile_10.write(problem_name + '_' + str(time) + '.post.bin\n')
            index_10 = 0

        if (index_50 == 50):
            multifile_50.write(problem_name + '_' + str(time) + '.post.bin\n')
            index_50 = 0

        index_5 += 1
        index_10 += 1
        index_50 += 1

        if (Multifile == "multiple_files"):
            gid_io.FinalizeResults()

        os.chdir(graphs_path)

        if (PrintNeighbourLists == "ON"):  # Printing neighbours id's
            os.chdir(neigh_list_path)
            neighbours_list = open('neigh_list_' + str(time), 'w')

            for elem in balls_model_part.Elements:
                ID = (elem.Id)
                Neigh_ID = elem.GetValue(NEIGHBOURS_IDS)

            for i in range(len(Neigh_ID)):
                neighbours_list.write(str(ID) + ' ' + str(Neigh_ID[i]) + '\n')

            neighbours_list.close()

        os.chdir(post_path)

        if (Multifile == "multiple_files"):
            gid_io.InitializeMesh(time)
            gid_io.WriteSphereMesh(balls_model_part.GetMesh())
            gid_io.FinalizeMesh()
            gid_io.InitializeResults(time, balls_model_part.GetMesh());

        proc.PrintingVariables(gid_io, export_model_part, time)
        os.chdir(main_path)

        time_old_print = time

    step += 1
#-------------------------------------------------------------------------------------------------------------------------------------


#-----------------------FINALITZATION OPERATIONS--------------------------------------------------------------------------------------

if (Multifile == "single_file"):
    gid_io.FinalizeResults()

multifile.close()
multifile_5.close()
multifile_10.close()
multifile_50.close()
os.chdir(main_path)

elapsed_pr_time = timer.clock() - initial_pr_time
elapsed_real_time = timer.time() - initial_real_time

print('Calculation ends at instant: ' + str(timer.time()))
print('Calculation ends at processing time instant: ' + str(timer.clock()))
print('Elapsed processing time: ' + str(elapsed_pr_time))
print('Elapsed real time: ' + str(elapsed_real_time))

print (my_timer)

print("ANALYSIS COMPLETED")

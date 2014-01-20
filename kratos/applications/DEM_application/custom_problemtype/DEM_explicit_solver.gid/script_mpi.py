from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import time as timer
import os
import sys
import math

#
#
#

from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *
from KratosMultiphysics.MetisApplication import *
from KratosMultiphysics.MPISearchApplication import *
from KratosMultiphysics.mpi import *

import DEM_explicit_solver_var as Param
import DEM_procedures

from DEM_procedures_mpi import *

proc = DEM_procedures.Procedures(Param)

#---------------------MODEL PART KRATOS AND GID.IO ------------------------------------------------------------------

# Defining a model part for the solid part

my_timer = Timer()
solid_model_part = ModelPart("SolidPart")

# Importing the strategy object

if(Param.ElementType == "SphericParticle3D" or Param.ElementType == "CylinderParticle2D"):
    import sphere_strategy as SolverStrategy

elif(Param.ElementType == "SphericContinuumParticle3D"):
    import continuum_sphere_strategy as SolverStrategy

SolverStrategy.AddVariables(solid_model_part, Param)

# Add mpi variables
AddMpiVariables(solid_model_part)

# reading the solid part: binary and single are forced in mpi.

gid_mode = GiDPostMode.GiD_PostBinary
multifile = MultiFileFlag.SingleFile

deformed_mesh_flag = WriteDeformedMeshFlag.WriteDeformed
write_conditions = WriteConditionsFlag.WriteConditions

gid_io = GidIO(problem_name, gid_mode, multifile, deformed_mesh_flag, write_conditions)
model_part_io_solid = ModelPartIO(problem_name)

# Perform the initial partition BEFORE reading
[model_part_io_solid, solid_model_part] = PerformInitialPartition(solid_model_part, model_part_io_solid, problem_name)

MPICommSetup = SetMPICommunicatorProcess(solid_model_part)
MPICommSetup.Execute()

model_part_io_solid.ReadModelPart(solid_model_part)

# Setting up the buffer size: SHOULD BE DONE AFTER READING!!!

solid_model_part.SetBufferSize(2)

# Adding dofs

SolverStrategy.AddDofs(solid_model_part)

# Constructing a creator/destructor object

creator_destructor = ParticleCreatorDestructor()

# Creating a solver object

solver = SolverStrategy.ExplicitStrategy(solid_model_part, creator_destructor, Param)  # here, solver variables initialize as default

first_print = True
index_5 = 1
index_10 = 1
index_50 = 1
prev_time = 0.0
control = 0.0

export_model_part = solid_model_part


#-------------------------------------------------------------------------------------------------------------------------

#------------------------------------------DEM_PROCEDURES FUNCTIONS & INITIALITZATION--------------------------------------------------------

if (Param.ModelDataInfo == "ON"):
    os.chdir(data_and_results)
    proc.ModelData(solid_model_part, solver)       # calculates the mean number of neighbours the mean radius, etc..
    os.chdir(main_path)

if(mpi.rank == 0):
    print('Initialitzating Problem....')

# MPI initialization
mpiutils = MpiUtilities()

solver.search_strategy = MPI_DEMSearch(solid_model_part.GetCommunicator())

mpiutils.Repart(solid_model_part, 0, 1)
mpiutils.CalculateModelNewIds(solid_model_part, 0)  # Only for bars?

solver.Initialize()

# Initialization of physics monitor and of the initial position of the center of mass

physics_calculator = SphericElementGlobalPhysicsCalculator(solid_model_part)

properties_list = []

if(mpi.rank == 0):
    print('Initialitzation Complete' + '\n')

dt = solid_model_part.ProcessInfo.GetValue(DELTA_TIME)

step = 0
time = 0.0
time_old_print = 0.0
initial_pr_time = timer.clock()
initial_real_time = timer.time()

#-------------------------------------------------------------------------------------------------------------------------------------

#-----------------------SINGLE FILE MESH AND RESULTS INITIALITZATION-------------------------------------------------------------------

#

gid_io.ChangeOutputName(problem_name + "_" + str(mpi.rank))

if (ContactMeshOption == "ON"):
    gid_io.InitializeMesh(0.0)
    gid_io.WriteMesh(contact_model_part.GetMesh())
    gid_io.FinalizeMesh()
    gid_io.InitializeResults(0.0, contact_model_part.GetMesh())

gid_io.InitializeMesh(0.0)
gid_io.WriteSphereMesh(solid_model_part.GetMesh())
gid_io.FinalizeMesh()
gid_io.InitializeResults(0.0, solid_model_part.GetMesh())

#------------------------------------------------------------------------------------------

#
#
# MAIN LOOP                                            #
#
#

if(mpi.rank == 0):
    print(('SOLVE starts at instant: ' + str(initial_pr_time) + '\n'))

total_steps_expected = int(Param.FinalTime / dt)

if(mpi.rank == 0):
    print(('Total number of TIME STEPs expected in the calculation are: ' + str(total_steps_expected) + ' if time step is kept ' + '\n'))

while (time < Param.FinalTime):

    dt = solid_model_part.ProcessInfo.GetValue(DELTA_TIME)  # Possible modifications of DELTA_TIME
    time = time + dt
    solid_model_part.CloneTimeStep(time)
    solid_model_part.ProcessInfo[TIME_STEPS] = step

    # _SOLVE_#########################################4
    solver.Solve()
    # TIME CONTROL######################################4

    incremental_time = (timer.time() - initial_real_time) - prev_time

    if (incremental_time > Param.ControlTime):
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

    total_force = 0.0
    force_node = 0.0

    # ___GiD IO____#########################################4

    time_to_print = time - time_old_print

    if (time_to_print >= Param.OutputTimeStep):

        proc.PrintingVariables(gid_io, export_model_part, time)

        time_old_print = time

    step += 1
#-------------------------------------------------------------------------------------------------------------------------------------


#-----------------------FINALITZATION OPERATIONS--------------------------------------------------------------------------------------
# proc.PlotPhysicalProperties(properties_list, graphs_path)

gid_io.FinalizeResults()

if(mpi.rank == 0):

    elapsed_pr_time = timer.clock() - initial_pr_time
    elapsed_real_time = timer.time() - initial_real_time

    print('Calculation ends at instant: ' + str(timer.time()))
    print('Calculation ends at processing time instant: ' + str(timer.clock()))
    print('Elapsed processing time: ' + str(elapsed_pr_time))
    print('Elapsed real time: ' + str(elapsed_real_time))

    print (my_timer)

    print("ANALYSIS COMPLETED")

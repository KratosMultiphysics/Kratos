from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import time as timer
import os
import sys
import math

# BENCHMARK ###
# BENCHMARK ###
# BENCHMARK ###
# BENCHMARK ###
# BENCHMARK ###

# BENCHMARK ###
# BENCHMARK ###
# BENCHMARK ###
# BENCHMARK ###
# BENCHMARK ###

from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *

# BENCHMARK ###
import DEM_explicit_solver_var as DEM_parameters
import DEM_procedures
proc = DEM_procedures.Procedures(DEM_parameters)

#import mesh_motion

# BENCHMARK ###

# BENCHMARK ###
# BENCHMARK ###
# BENCHMARK ###
# BENCHMARK ###
# BENCHMARK ###
# BENCHMARK ###
# BENCHMARK ###
# BENCHMARK ###
# BENCHMARK ###

#---------------------MODEL PART KRATOS AND GID.IO ------------------------------------------------------------------

# Defining a model part for the solid part

my_timer = Timer()
balls_model_part = ModelPart("SolidPart")

RigidFace_model_part = ModelPart("RigidFace_Part")
mixed_model_part = ModelPart("Mixed_Part")

RigidFace_model_part.AddNodalSolutionStepVariable(VELOCITY)
RigidFace_model_part.AddNodalSolutionStepVariable(DISPLACEMENT)
RigidFace_model_part.AddNodalSolutionStepVariable(TOTAL_FORCES)
if(DEM_parameters.PostGroupId):
  RigidFace_model_part.AddNodalSolutionStepVariable(GROUP_ID)

# Importing the strategy object

import sphere_strategy as SolverStrategy

SolverStrategy.AddVariables(balls_model_part, DEM_parameters)

# Reading the model_part: binary or ascii, multifile or single --> only binary and single for mpi.

if (DEM_parameters.OutputFileType == "Binary"):
    gid_mode = GiDPostMode.GiD_PostBinary

else:
    gid_mode = GiDPostMode.GiD_PostAscii

if (DEM_parameters.Multifile == "multiple_files"):
    multifile = MultiFileFlag.MultipleFiles

else:
    multifile = MultiFileFlag.SingleFile

deformed_mesh_flag = WriteDeformedMeshFlag.WriteDeformed
write_conditions = WriteConditionsFlag.WriteConditions

gid_io = GidIO(DEM_parameters.problem_name, gid_mode, multifile, deformed_mesh_flag, write_conditions)
model_part_io_solid = ModelPartIO(DEM_parameters.problem_name + "DEM", True)
model_part_io_solid.ReadModelPart(balls_model_part)

rigidFace_mp_filename = DEM_parameters.problem_name + "DEM_FEM_boundary"
model_part_io_solid = ModelPartIO(rigidFace_mp_filename)
model_part_io_solid.ReadModelPart(RigidFace_model_part)

# Setting up the buffer size: SHOULD BE DONE AFTER READING!!!

balls_model_part.SetBufferSize(1)

# Adding dofs

SolverStrategy.AddDofs(balls_model_part)

# Constructing a creator/destructor object

creator_destructor = ParticleCreatorDestructor()

# Creating a solver object

solver = SolverStrategy.ExplicitStrategy(balls_model_part, RigidFace_model_part, creator_destructor, DEM_parameters)  # here, solver variables initialize as default

# Creating necessary directories

main_path = os.getcwd()
post_path = str(main_path) + '/' + str(DEM_parameters.problem_name) + '_Post_Files'
list_path = str(main_path) + '/' + str(DEM_parameters.problem_name) + '_Post_Lists'
data_and_results = str(main_path) + '/' + str(DEM_parameters.problem_name) + '_Results_and_Data'
graphs_path = str(main_path) + '/' + str(DEM_parameters.problem_name) + '_Graphs'
MPI_results = str(main_path) + '/' + str(DEM_parameters.problem_name) + '_MPI_results'

for directory in [post_path, list_path, data_and_results, graphs_path, MPI_results]:

    if not os.path.isdir(directory):
        os.makedirs(str(directory))

os.chdir(list_path)

multifile = open(DEM_parameters.problem_name + '_all' + '.post.lst', 'w')
multifile_5 = open(DEM_parameters.problem_name + '_5' + '.post.lst', 'w')
multifile_10 = open(DEM_parameters.problem_name + '_10' + '.post.lst', 'w')
multifile_50 = open(DEM_parameters.problem_name + '_50' + '.post.lst', 'w')

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

print('Initializing Problem....')
sys.stdout.flush()

solver.Initialize()


# constructing a model part for the DEM inlet. it contains the DEM elements to be released during the simulation
inlet_option                     = 1
dem_inlet_element_type           = "SphericParticle3D"  # "SphericParticle3D", "SphericSwimmingParticle3D"

if (inlet_option):
    max_node_Id = DEM_procedures.FindMaxNodeIdInModelPart(balls_model_part)
    max_FEM_node_Id = DEM_procedures.FindMaxNodeIdInModelPart(RigidFace_model_part)
    if ( max_FEM_node_Id > max_node_Id):
        max_node_Id = max_FEM_node_Id
    creator_destructor.SetMaxNodeId(max_node_Id)
        
    DEM_inlet_model_part = ModelPart("DEMInletPart")
    DEM_Inlet_filename = DEM_parameters.problem_name + "DEM_Inlet"
    SolverStrategy.AddVariables(DEM_inlet_model_part, DEM_parameters)
    model_part_io_demInlet = ModelPartIO(DEM_Inlet_filename)
    model_part_io_demInlet.ReadModelPart(DEM_inlet_model_part)

    # setting up the buffer size:
    DEM_inlet_model_part.SetBufferSize(1)

    # adding nodal degrees of freedom
    SolverStrategy.AddDofs(DEM_inlet_model_part)
    DEM_inlet_parameters = DEM_inlet_model_part.Properties

    # constructing the inlet and intializing it (must be done AFTER the balls_model_part Initialize)
    DEM_inlet = DEM_Inlet(DEM_inlet_model_part)    
    DEM_inlet.InitializeDEM_Inlet(balls_model_part, creator_destructor)



#-------------------------------------------------------------------------------------------------------------------------

#------------------------------------------DEM_PROCEDURES FUNCTIONS & INITIALIZATIONS--------------------------------------------------------


# Initialization of physics monitor and of the initial position of the center of mass
# physics_calculator = SphericElementGlobalPhysicsCalculator(balls_model_part)
# properties_list = []
print('Initialitzation Complete' + '\n')
sys.stdout.flush()

# BENCHMARK ###
# BENCHMARK ###
# BENCHMARK ###
# BENCHMARK ###


step = 0
time = 0.0
time_old_print = 0.0
initial_pr_time = timer.clock()
initial_real_time = timer.time()

#-------------------------------------------------------------------------------------------------------------------------------------

#-----------------------SINGLE FILE MESH AND RESULTS INITIALITZATION-------------------------------------------------------------------

post_utility = PostUtilities()

os.chdir(post_path)

if (DEM_parameters.Multifile == "single_file"):

    post_utility.AddModelPartToModelPart(mixed_model_part, balls_model_part)
    post_utility.AddModelPartToModelPart(mixed_model_part, RigidFace_model_part)
    gid_io.InitializeMesh(0.0)
    gid_io.WriteMesh(RigidFace_model_part.GetMesh())
    gid_io.WriteSphereMesh(balls_model_part.GetMesh())
    gid_io.FinalizeMesh()
    gid_io.InitializeResults(0.0, mixed_model_part.GetMesh())


# MODEL DATA

if (DEM_parameters.ModelDataInfo == "ON"):
    os.chdir(data_and_results)
    proc.ModelData(balls_model_part, balls_model_part, solver)  # dummy contact model part. (only for continuum)      # calculates the mean number of neighbours the mean radius, etc..
    os.chdir(main_path)


#------------------------------------------------------------------------------------------

#
#
# MAIN LOOP                                            #
#
#
os.chdir(main_path)

dt = balls_model_part.ProcessInfo.GetValue(DELTA_TIME)

total_steps_expected = int(DEM_parameters.FinalTime / dt)

print(('Main loop starts at instant: ' + str(initial_pr_time) + '\n'))

print(('Total number of TIME STEPs expected in the calculation are: ' + str(total_steps_expected) + ' if time step is kept ' + '\n'))
sys.stdout.flush()

mesh_motion = DEMFEMUtilities()

while (time < DEM_parameters.FinalTime):

    dt = balls_model_part.ProcessInfo.GetValue(DELTA_TIME)  # Possible modifications of DELTA_TIME
    time = time + dt
    # balls_model_part.CloneTimeStep(time)
    balls_model_part.ProcessInfo[TIME] = time
    balls_model_part.ProcessInfo[DELTA_TIME] = dt
    balls_model_part.ProcessInfo[TIME_STEPS] = step
    
    #walls movement:
    mesh_motion.MoveAllMeshes(RigidFace_model_part, time)
    
    # _SOLVE_###########################################
    os.chdir(main_path)
    solver.Solve()
    # TIME CONTROL######################################
        
    # adding DEM elements by the inlet:
    if (inlet_option):
        DEM_inlet.CreateElementsFromInletMesh(balls_model_part, DEM_inlet_model_part, creator_destructor, dem_inlet_element_type)  # After solving, to make sure that neighbours are already set.        

    incremental_time = (timer.time() - initial_real_time) - prev_time

    if (incremental_time > DEM_parameters.ControlTime):
        percentage = 100.0 * (float(step) / total_steps_expected)
        
        print("%s %.2f %s" % ("Real time calculation: ", timer.time() - initial_real_time,"s"))      
        print('Simulation time: ' + str(time))
        print("%s %.5f %s" % ("Percentage Completed: ", percentage,"%"))        
        print("Time Step: " + str(step) + '\n')

        sys.stdout.flush()

        prev_time = (timer.time() - initial_real_time)
    
    if ((timer.time() - initial_real_time > 60) and first_print == True and step != 0):
        first_print = False
        estimated_sim_duration = 60.0 * (total_steps_expected / step) # seconds

        print('The calculation total estimated time is ' + str(estimated_sim_duration) + 'seconds' + '\n')
        print('in minutes:'        + str(estimated_sim_duration / 60.0) + 'min.' + '\n')
        print('in hours:'        + str(estimated_sim_duration / 3600.0) + 'hrs.' + '\n')
        print('in days:'        + str(estimated_sim_duration / 86400.0) + 'days' + '\n') 
        sys.stdout.flush()

        if (estimated_sim_duration / 86400 > 2.0):

          print('WARNING!!!:       VERY LASTING CALCULATION' + '\n')

    # CONCRETE_TEST_STUFF#########################################4

    os.chdir(list_path)
    multifile.write(DEM_parameters.problem_name + '_' + str(time) + '.post.bin\n')
    os.chdir(main_path)

    # ___GiD IO____#########################################4

    time_to_print = time - time_old_print

    if (time_to_print >= DEM_parameters.OutputTimeStep):
        
        print("")
        print("*******************  PRINTING RESULTS FOR GID  ***************************")
        sys.stdout.flush()
        
        # BENCHMARK ###
        os.chdir(data_and_results)

        # properties_list = proc.MonitorPhysicalProperties(balls_model_part, physics_calculator, properties_list)

        if (index_5 == 5):
            multifile_5.write(DEM_parameters.problem_name + '_' + str(time) + '.post.bin\n')
            index_5 = 0

        if (index_10 == 10):
            multifile_10.write(DEM_parameters.problem_name + '_' + str(time) + '.post.bin\n')
            index_10 = 0

        if (index_50 == 50):
            multifile_50.write(DEM_parameters.problem_name + '_' + str(time) + '.post.bin\n')
            index_50 = 0

        index_5 += 1
        index_10 += 1
        index_50 += 1

        if (DEM_parameters.Multifile == "multiple_files"):
            gid_io.FinalizeResults()

        os.chdir(post_path)

        if (DEM_parameters.Multifile == "multiple_files"):
            mixed_model_part.Elements.clear()
            mixed_model_part.Nodes.clear()

            post_utility.AddModelPartToModelPart(mixed_model_part, balls_model_part)
            post_utility.AddModelPartToModelPart(mixed_model_part, RigidFace_model_part)

            gid_io.InitializeMesh(time)
            gid_io.WriteSphereMesh(balls_model_part.GetMesh())
            gid_io.WriteMesh(RigidFace_model_part.GetMesh())
            gid_io.FinalizeMesh()

            gid_io.InitializeResults(time, mixed_model_part.GetMesh())

        proc.PrintingGlobalVariables(gid_io, mixed_model_part, time)
        proc.PrintingBallsVariables(gid_io, balls_model_part, time)

        if (DEM_parameters.Multifile == "multiple_files"):
            gid_io.FinalizeResults()

        time_old_print = time

    step += 1
#-------------------------------------------------------------------------------------------------------------------------------------


#-----------------------FINALITZATION OPERATIONS--------------------------------------------------------------------------------------
# proc.PlotPhysicalProperties(properties_list, graphs_path)

if (DEM_parameters.Multifile == "single_file"):
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

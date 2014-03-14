from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import time as timer
import os
import sys
import math

from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *
from KratosMultiphysics.MetisApplication import *
from KratosMultiphysics.MPISearchApplication import *
from KratosMultiphysics.mpi import *

import DEM_explicit_solver_var as DEM_parameters
import DEM_procedures
Procedures = DEM_procedures.Procedures(DEM_parameters)
import DEM_procedures_mpi as DEM_procedures_mpi

import DEM_material_test_script_mpi 
import mesh_motion
import MPIer

#---------------------MODEL PART KRATOS AND GID.IO ------------------------------------------------------------------

# Defining a model part for the solid part

if (mpi.rank == 0):
  MPIClassObject = MPIer.MPIerClass(str(DEM_parameters.problem_name) + "DEM.mdpa")

my_timer = Timer();
balls_model_part = ModelPart("SpheresPart");

RigidFace_model_part   = ModelPart("RigidFace_Part");  
mixed_model_part       = ModelPart("Mixed_Part");

RigidFace_model_part.AddNodalSolutionStepVariable(VELOCITY)
RigidFace_model_part.AddNodalSolutionStepVariable(DISPLACEMENT)
RigidFace_model_part.AddNodalSolutionStepVariable(ELASTIC_FORCES)
RigidFace_model_part.AddNodalSolutionStepVariable(TOTAL_FORCES)
RigidFace_model_part.AddNodalSolutionStepVariable(GROUP_ID)
RigidFace_model_part.AddNodalSolutionStepVariable(EXPORT_GROUP_ID)

# Importing the strategy object

import continuum_sphere_strategy as SolverStrategy

SolverStrategy.AddVariables(balls_model_part, DEM_parameters)

DEM_procedures_mpi.AddMpiVariables(balls_model_part)

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
spheres_mp_filename = DEM_parameters.problem_name + "DEM"

model_part_io_spheres = ModelPartIO(spheres_mp_filename)

# Perform the initial partition BEFORE reading

[model_part_io_spheres, balls_model_part] = DEM_procedures_mpi.PerformInitialPartition(balls_model_part, model_part_io_spheres, spheres_mp_filename)

MPICommSetup = SetMPICommunicatorProcess(balls_model_part)
MPICommSetup.Execute()

def KRATOSprint(message):
    if (mpi.rank == 0):
        print(message)    

model_part_io_spheres.ReadModelPart(balls_model_part)

rigidFace_mp_filename = DEM_parameters.problem_name + "DEM_FEM_boundary"

model_part_io_fem = ModelPartIO(rigidFace_mp_filename)
model_part_io_fem.ReadModelPart(RigidFace_model_part)

# Setting up the buffer size: SHOULD BE DONE AFTER READING!!!

balls_model_part.SetBufferSize(1)

# Adding dofs

SolverStrategy.AddDofs(balls_model_part)

# Constructing a creator/destructor object

creator_destructor = ParticleCreatorDestructor()

# Creating a solver object

solver = SolverStrategy.ExplicitStrategy(balls_model_part, RigidFace_model_part,creator_destructor, DEM_parameters) #here, solver variables initialize as default

# Creating necessary directories

main_path        = os.getcwd()
post_path        = str(main_path) + '/' + str(DEM_parameters.problem_name) + '_Post_Files'
list_path        = str(main_path) + '/' + str(DEM_parameters.problem_name) + '_Post_Lists'
data_and_results = str(main_path) + '/' + str(DEM_parameters.problem_name) + '_Results_and_Data'
graphs_path      = str(main_path) + '/' + str(DEM_parameters.problem_name) + '_Graphs'
MPI_results      = str(main_path) + '/' + str(DEM_parameters.problem_name) + '_MPI_results'

for directory in [post_path, list_path, data_and_results, graphs_path, MPI_results]:

    if not os.path.isdir(directory):
        os.makedirs(str(directory))

os.chdir(list_path)

multifile        = open(DEM_parameters.problem_name + '_all' + '.post.lst', 'w')
multifile_5      = open(DEM_parameters.problem_name + '_5'   + '.post.lst', 'w')
multifile_10     = open(DEM_parameters.problem_name + '_10'  + '.post.lst', 'w')
multifile_50     = open(DEM_parameters.problem_name + '_50'  + '.post.lst', 'w')

multifile.write('Multiple\n')
multifile_5.write('Multiple\n')
multifile_10.write('Multiple\n')
multifile_50.write('Multiple\n')

os.chdir(main_path)

KRATOSprint ("Initializing Problem....")

# MPI initialization
mpiutils = MpiUtilities();

solver.search_strategy = MPI_DEMSearch(balls_model_part.GetCommunicator())

mpiutils.Repart(balls_model_part, 0, 1)
mpiutils.CalculateModelNewIds(balls_model_part, 0)




solver.Initialize()

if ( DEM_parameters.ContactMeshOption =="ON" ) :

  contact_model_part = solver.contact_model_part

#-----------------------SINGLE FILE MESH AND RESULTS INITIALITZATION-------------------------------------------------------------------

post_utility = PostUtilities()

gid_io.ChangeOutputName(DEM_parameters.problem_name + "_" + str(mpi.rank))

os.chdir(post_path)

if (DEM_parameters.Multifile == "single_file"):

  post_utility.AddModelPartToModelPart(mixed_model_part, balls_model_part)
  if (DEM_parameters.ContactMeshOption == "ON"):
      post_utility.AddModelPartToModelPart(mixed_model_part, contact_model_part)
  post_utility.AddModelPartToModelPart(mixed_model_part, RigidFace_model_part)
  gid_io.InitializeMesh(0.0) 
  gid_io.WriteMesh(RigidFace_model_part.GetMesh())
  gid_io.WriteSphereMesh(balls_model_part.GetMesh())
  if (DEM_parameters.ContactMeshOption == "ON"):
      gid_io.WriteMesh(contact_model_part.GetMesh())
  gid_io.FinalizeMesh()
  gid_io.InitializeResults(0.0, mixed_model_part.GetMesh())
  
#------------------------------------------DEM_PROCEDURES FUNCTIONS & INITIALIZATIONS--------------------------------------------------------

#if (DEM_parameters.PredefinedSkinOption == "ON" ):

   #ProceduresSetPredefinedSkin(balls_model_part)

if(DEM_parameters.TestType != "None"):
 
 MaterialTest = DEM_material_test_script_mpi.MaterialTest(DEM_parameters, Procedures, solver, graphs_path, post_path, balls_model_part, RigidFace_model_part)
 
KRATOSprint ("Initialization Complete" + "\n")

step                   = 0
time                   = 0.0
time_old_print         = 0.0
initial_pr_time        = timer.clock()
initial_real_time      = timer.time()

first_print  = True; index_5 = 1; index_10  = 1; index_50  = 1; prev_time  = 0.0; control = 0.0
    
##MODEL DATA 

if (DEM_parameters.ModelDataInfo == "ON"):
    os.chdir(data_and_results)
    if (DEM_parameters.ContactMeshOption == "ON"):
      (coordination_number) = Procedures.ModelData(balls_model_part, contact_model_part, solver)       # calculates the mean number of neighbours the mean radius, etc..
      KRATOSprint ("Coordination Number: " + str(coordination_number) + "\n")
      os.chdir(main_path)
    else:
      KRATOSprint("Activate Contact Mesh for ModelData information")

if(DEM_parameters.Dempack and (DEM_parameters.TestType != "None")):
  
 if(mpi.rank == 0):
    MaterialTest.PrintChart();
 MaterialTest.PrepareDataForGraph()
 
#------------------------------------------------------------------------------------------
 
###########################################################################################
#                                                                                         #
#                                    MAIN LOOP                                            #
#                                                                                         #
###########################################################################################

dt = balls_model_part.ProcessInfo.GetValue(DELTA_TIME) # Possible modifications of DELTA_TIME

total_steps_expected = int(DEM_parameters.FinalTime / dt)

KRATOSprint ("Main loop starts at instant: " + str(initial_pr_time) + "\n")

KRATOSprint ("Total number of TIME STEPs expected in the calculation are: " + str(total_steps_expected) + " if time step is kept " + "\n")

#if(DEM_parameters.PoissonMeasure == "ON"):
  #MaterialTest.PoissonMeasure()
  
while (time < DEM_parameters.FinalTime):

    dt = balls_model_part.ProcessInfo.GetValue(DELTA_TIME) # Possible modifications of DELTA_TIME
    time = time + dt
    
    balls_model_part.ProcessInfo[TIME] = time
    balls_model_part.ProcessInfo[DELTA_TIME] = dt
    balls_model_part.ProcessInfo[TIME_STEPS] = step
    
    RigidFace_model_part.ProcessInfo[TIME] = time
    RigidFace_model_part.ProcessInfo[DELTA_TIME] = dt
    RigidFace_model_part.ProcessInfo[TIME_STEPS] = step
    
    #walls movement:
    mesh_motion.MoveAllMeshes(RigidFace_model_part, time)
    
    #########################_SOLVE_#########################################4

    solver.Solve()
    
    #########################TIME CONTROL######################################4
    
    incremental_time = (timer.time() - initial_real_time) - prev_time

    if (incremental_time > DEM_parameters.ControlTime):
        
        percentage = 100 * (float(step) / total_steps_expected)

        KRATOSprint('Real time calculation: ' + str(timer.time() - initial_real_time))
        KRATOSprint('Simulation time: ' + str(time))
        KRATOSprint("%s %.5f %s" % ("Percentage Completed: ", percentage,"%"))      
        KRATOSprint("Time Step: " + str(step) + '\n')
        sys.stdout.flush()

        prev_time = (timer.time() - initial_real_time)
    
    if ((timer.time() - initial_real_time > 60) and first_print == True and step != 0):
        first_print = False
        estimated_sim_duration = 60.0 * (total_steps_expected / step) # seconds

        KRATOSprint('The calculation total estimated time is ' + str(estimated_sim_duration) + 'seconds' + '\n')
        KRATOSprint('in minutes:'        + str(estimated_sim_duration / 60.0) + 'min.' + '\n')
        KRATOSprint('in hours:'        + str(estimated_sim_duration / 3600.0) + 'hrs.' + '\n')
        KRATOSprint('in days:'        + str(estimated_sim_duration / 86400.0) + 'days' + '\n') 
        sys.stdout.flush()

        if (estimated_sim_duration / 86400 > 2.0):

          KRATOSprint('WARNING!!!:       VERY LASTING CALCULATION' + '\n')

    #########################CONCRETE_TEST_STUFF#########################################
    
    if( DEM_parameters.TestType != "None"):
      
      MaterialTest.MeasureForcesAndPressure()
      
      if(mpi.rank == 0):
        MaterialTest.PrintGraph(step)

    ##########################___GiD IO____#########################################
    
    os.chdir(list_path)    
    multifile.write(DEM_parameters.problem_name + '_' + str(time) + '.post.bin\n')   
    os.chdir(main_path)

    time_to_print = time - time_old_print

    if (time_to_print >= DEM_parameters.OutputTimeStep):

        os.chdir(data_and_results)

        #properties_list = ProceduresMonitorPhysicalProperties(balls_model_part, physics_calculator, properties_list)

        if (index_5 == 5):
            multifile_5.write(DEM_parameters.problem_name + '_' + str(time) + '.post.bin\n')
            index_5 = 0

        if (index_10 == 10):
            multifile_10.write(DEM_parameters.problem_name + '_' + str(time) + '.post.bin\n')
            index_10 = 0

        if (index_50 == 50):
            multifile_50.write(DEM_parameters.problem_name + '_' + str(time) + '.post.bin\n')
            index_50 = 0

        index_5  += 1
        index_10 += 1
        index_50 += 1

        os.chdir(post_path)

        if (DEM_parameters.Multifile == "multiple_files"):
            mixed_model_part.Elements.clear()
            mixed_model_part.Nodes.clear()

            post_utility.AddModelPartToModelPart(mixed_model_part, balls_model_part)
            if (DEM_parameters.ContactMeshOption == "ON"):
                post_utility.AddModelPartToModelPart(mixed_model_part, contact_model_part)
            post_utility.AddModelPartToModelPart(mixed_model_part, RigidFace_model_part)

            gid_io.InitializeMesh(time) 
            gid_io.WriteSphereMesh(balls_model_part.GetMesh())
            if (DEM_parameters.ContactMeshOption == "ON"):
                gid_io.WriteMesh(contact_model_part.GetMesh())
            gid_io.WriteMesh(RigidFace_model_part.GetMesh())
            gid_io.FinalizeMesh()
            
            gid_io.InitializeResults(time, mixed_model_part.GetMesh())
                                
        Procedures.PrintingGlobalVariables(gid_io, mixed_model_part, time)
        Procedures.PrintingBallsVariables(gid_io, balls_model_part, time)
        
        if (DEM_parameters.ContactMeshOption == "ON"):
            Procedures.PrintingContactElementsVariables(gid_io, contact_model_part, time)
        
        os.chdir(main_path)     
              
        if (DEM_parameters.Multifile == "multiple_files"):
            gid_io.FinalizeResults() 

        time_old_print = time
  
    step += 1

    #if((step%500) == 0):
      #if (( DEM_parameters.ContactMeshOption =="ON") and (DEM_parameters.TestType!= "None"))  :
          #MaterialTest.OrientationStudy(contact_model_part, step)
    

#-----------------------FINALIZATION OPERATIONS-------------------------------------------------------------------------------------- 

if (DEM_parameters.Multifile == "single_file"):
    gid_io.FinalizeResults()

if ((DEM_parameters.TestType!= "None") and (mpi.rank == 0) ):
  
  MaterialTest.FinalizeGraphs(DEM_parameters)

multifile.close()
multifile_5.close()
multifile_10.close()
multifile_50.close()
os.chdir(main_path)

elapsed_pr_time     = timer.clock() - initial_pr_time
elapsed_real_time   = timer.time() - initial_real_time

KRATOSprint ('Calculation ends at instant: '                 + str(timer.time()))
KRATOSprint ('Calculation ends at processing time instant: ' + str(timer.clock()))
KRATOSprint ('Elapsed processing time: '                     + str(elapsed_pr_time))
KRATOSprint ('Elapsed real time: '                           + str(elapsed_real_time))

#KRATOSprint (my_timer)

KRATOSprint ("ANALYSIS COMPLETED" )
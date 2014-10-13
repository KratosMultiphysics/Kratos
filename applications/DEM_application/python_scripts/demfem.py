from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Python imports
import time as timer
import os
import sys
import math

# Kratos
from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *
from KratosMultiphysics.SolidMechanicsApplication import * #DEMFEM
from KratosMultiphysics.ExternalSolversApplication import * #DEMFEM

# DEM Application
import DEM_explicit_solver_var as DEM_parameters

#FEM IMPORTS

import ProjectParameters as FEM_general_variables #DEMFEM
import restart_utility as restart_utils
import gid_output_utility as gid_utils
import conditions_python_utility as condition_utils
import list_files_python_utility as files_utils
import time_operation_utility as operation_utils

# TODO: Ungly fix. Change it. I don't like this to be in the main...
# Strategy object
if   (DEM_parameters.ElementType == "SphericPartDEMElement3D"     or DEM_parameters.ElementType == "CylinderPartDEMElement2D"):
    import sphere_strategy as SolverStrategy
elif (DEM_parameters.ElementType == "SphericContPartDEMElement3D" or DEM_parameters.ElementType == "CylinderContPartDEMElement2D"):
    import continuum_sphere_strategy as SolverStrategy

# Import MPI modules if needed
if os.environ.has_key("OMPI_COMM_WORLD_SIZE"):
    # Kratos MPI
    from KratosMultiphysics.MetisApplication import *
    from KratosMultiphysics.MPISearchApplication import *
    from KratosMultiphysics.mpi import *

    # DEM Application MPI
    import DEM_procedures_mpi as DEM_procedures
    import DEM_material_test_script_mpi as DEM_material_test_script
else :
    # DEM Application
    import DEM_procedures
    import DEM_material_test_script

    print("Runing under OpenMP")

##############################################################################
#                                                                            #
#    INITIALIZE                                                              #
#                                                                            #
##############################################################################

# Import utilities from moduels
procedures    = DEM_procedures.Procedures(DEM_parameters)
demio         = DEM_procedures.DEMIo()
report        = DEM_procedures.Report()
parallelutils = DEM_procedures.ParallelUtils()
materialTest  = DEM_procedures.MaterialTest()
 
# Set the print function TODO: do this better...
KRATOSprint   = procedures.KRATOSprint

# Preprocess the model
procedures.PreProcessModel(DEM_parameters)

# Prepare modelparts
dem_model_part = ModelPart("DEM_Part");
fem_model_part = ModelPart("FEM_Part");
mixed_model_part = ModelPart("Mixed_Part");
contact_model_part    = ""

# Add variables
procedures.AddCommonVariables(dem_model_part, DEM_parameters)
SolverStrategy.AddVariables(dem_model_part, DEM_parameters)
procedures.AddMpiVariables(dem_model_part)
procedures.AddCommonVariables(fem_model_part, DEM_parameters)
SolverStrategy.AddFEMVariables(fem_model_part, DEM_parameters)
procedures.AddMpiVariables(fem_model_part)

# Reading the dem_model_part
DEM_model_part_filename = DEM_parameters.problem_name + "_DEM"
model_part_io_DEM = ModelPartIO(DEM_model_part_filename,True)

# Perform the initial partition
[model_part_io_DEM, dem_model_part, MPICommSetup] = parallelutils.PerformInitialPartition(dem_model_part, model_part_io_DEM, DEM_model_part_filename)

model_part_io_DEM.ReadModelPart(dem_model_part)

FEM_model_part_filename = DEM_parameters.problem_name + "_FEM"

model_part_io_FEM = ModelPartIO(FEM_model_part_filename)
model_part_io_FEM.ReadModelPart(fem_model_part)

FEMSolverSettings = FEM_general_variables.SolverSettings #DEMFEM

# import solver file
FEM_solver_constructor = __import__(FEMSolverSettings.solver_type)  #DEMFEM

# construct the solver
FEM_main_step_solver = FEM_solver_constructor.CreateSolver(fem_model_part, FEMSolverSettings) #DEMFEM

FEM_solver_constructor.AddVariables(fem_model_part, FEMSolverSettings) #DEMFEM
#add variables from dem to fem
fem_model_part.AddNodalSolutionStepVariable(ELASTIC_FORCES)
fem_model_part.AddNodalSolutionStepVariable(TOTAL_FORCES)
if(DEM_parameters.PostGroupId):
  fem_model_part.AddNodalSolutionStepVariable(GROUP_ID)

# defining the type, the name and the path of the problem:
problem_type = FEM_general_variables.ProblemType
dem_problem_name = FEM_general_variables.problem_name
fem_problem_name = DEM_parameters.problem_name

if(dem_problem_name != fem_problem_name):
   print("----------------------MIQUEL: only a unique problem name-----------------------")
else:
  problem_name = dem_problem_name
  
problem_path = FEM_general_variables.problem_path

deformed_mesh_flag = WriteDeformedMeshFlag.WriteDeformed
write_conditions = WriteConditionsFlag.WriteConditions

FEM_solver_constructor.AddDofs(fem_model_part, FEMSolverSettings)

# Adding dofs
SolverStrategy.AddDofs(dem_model_part)

# Setting up the buffer size: SHOULD BE DONE AFTER READING!!!

dem_buffer_size = 1
fem_buffer_size = 2

dem_model_part.SetBufferSize(dem_buffer_size)
fem_model_part.SetBufferSize(fem_buffer_size)

# Reading the model_part: binary or ascii, multifile or single --> only binary and single for mpi.

if (DEM_parameters.OutputFileType == "Binary"):
    gid_mode = GiDPostMode.GiD_PostBinary

else:
    gid_mode = GiDPostMode.GiD_PostAscii

if (DEM_parameters.Multifile == "multiple_files"):
    multifile = MultiFileFlag.MultipleFiles

else:
    multifile = MultiFileFlag.SingleFile
gid_io = GidIO(DEM_parameters.problem_name, gid_mode, multifile, deformed_mesh_flag, write_conditions)

# Constructing a creator/destructor object
creator_destructor = ParticleCreatorDestructor()

# Creating a solver object

solver = SolverStrategy.ExplicitStrategy(dem_model_part, fem_model_part, creator_destructor, DEM_parameters)  # here, solver variables initialize as default

#-----------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------
#FEM OPERATIONS

print_lists = FEM_general_variables.PrintLists
output_mode = FEM_general_variables.GidOutputConfiguration.GiDPostMode
list_files = files_utils.ListFilesUtility(problem_path, problem_name, print_lists, output_mode)
list_files.Initialize(FEM_general_variables.file_list)

domain_size = FEM_general_variables.domain_size

problem_restart = restart_utils.RestartUtility(fem_model_part, problem_path, problem_name)
#problem_restart.CleanPreviousFiles()
list_files = files_utils.ListFilesUtility(problem_path, problem_name, print_lists, output_mode)
list_files.Initialize(FEM_general_variables.file_list)
#list_files.RemoveListFiles()
gid_print = gid_utils.GidOutputUtility(problem_name, FEM_general_variables.GidOutputConfiguration)

# set the constitutive law
import constitutive_law_python_utility as constitutive_law_utils

constitutive_law = constitutive_law_utils.ConstitutiveLawUtility(fem_model_part, domain_size);
constitutive_law.Initialize();

# --DEFINE CONDITIONS START--#################
incr_disp = FEM_general_variables.Incremental_Displacement
incr_load = FEM_general_variables.Incremental_Load
rotation_dofs = FEMSolverSettings.RotationDofs
conditions = condition_utils.ConditionsUtility(fem_model_part, domain_size, incr_disp, incr_load, rotation_dofs)


#-----------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------


# Creating necessary directories
main_path = os.getcwd()
[post_path,list_path,data_and_results,graphs_path,MPI_results] = procedures.CreateDirectories(str(main_path),str(DEM_parameters.problem_name))

os.chdir(list_path)

multifiles = (
    DEM_procedures.MultifileList(DEM_parameters.problem_name,1 ),
    DEM_procedures.MultifileList(DEM_parameters.problem_name,5 ),
    DEM_procedures.MultifileList(DEM_parameters.problem_name,10),
    DEM_procedures.MultifileList(DEM_parameters.problem_name,50),
    )

demio.SetMultifileLists(multifiles)
#prev_time = 0.0

os.chdir(main_path)

KRATOSprint("Initializing Problem....")

# Initialize GiD-IO
demio.AddGlobalVariables()
demio.AddBallVariables()
demio.AddContactVariables()
demio.AddMpiVariables()
demio.EnableMpiVariables()

demio.Configure(DEM_parameters.problem_name,
                DEM_parameters.OutputFileType,
                DEM_parameters.Multifile,
                DEM_parameters.ContactMeshOption)

demio.SetOutputName(DEM_parameters.problem_name)

os.chdir(post_path)
demio.InitializeMesh(mixed_model_part,
                     dem_model_part,
                     fem_model_part,
                     contact_model_part)

#initial_pr_time = timer.clock()
#initial_real_time = timer.time()

os.chdir(post_path)
#demio.PrintResults(2)

# Perform a partition to balance the problem
parallelutils.Repart(dem_model_part)
parallelutils.CalculateModelNewIds(dem_model_part)

os.chdir(post_path)
#demio.PrintResults(3)

# Creating a solver object and set the search strategy
solver                 = SolverStrategy.ExplicitStrategy(dem_model_part, fem_model_part,creator_destructor, DEM_parameters)
solver.search_strategy = parallelutils.GetSearchStrategy(solver, dem_model_part)
solver.Initialize()

if ( DEM_parameters.ContactMeshOption =="ON" ) :
    contact_model_part = solver.contact_model_part
  
# constructing a model part for the DEM inlet. it contains the DEM elements to be released during the simulation  
# Initializing the DEM solver must be done before creating the DEM Inlet, because the Inlet configures itself according to some options of the DEM model part
inlet_option                     = 1
dem_inlet_element_type           = "SphericParticle3D"  # "SphericParticle3D", "SphericSwimmingParticle3D"

if (inlet_option):
    max_node_Id = DEM_procedures.FindMaxNodeIdInModelPart(dem_model_part)
    max_FEM_node_Id = DEM_procedures.FindMaxNodeIdInModelPart(fem_model_part)
    if ( max_FEM_node_Id > max_node_Id):
        max_node_Id = max_FEM_node_Id
    creator_destructor.SetMaxNodeId(max_node_Id)
        
    DEM_inlet_model_part = ModelPart("DEMInletPart")
    DEM_Inlet_filename = DEM_parameters.problem_name + "_DEM_Inlet"
    SolverStrategy.AddVariables(DEM_inlet_model_part, DEM_parameters)
    os.chdir(main_path)
    model_part_io_demInlet = ModelPartIO(DEM_Inlet_filename)
    model_part_io_demInlet.ReadModelPart(DEM_inlet_model_part)

    # setting up the buffer size:
    DEM_inlet_model_part.SetBufferSize(1)

    # adding nodal degrees of freedom
    SolverStrategy.AddDofs(DEM_inlet_model_part)
    DEM_inlet_parameters = DEM_inlet_model_part.Properties
    
    for properties in DEM_inlet_model_part.Properties:
            
            DiscontinuumConstitutiveLawString = properties[DEM_DISCONTINUUM_CONSTITUTIVE_LAW_NAME];
            DiscontinuumConstitutiveLaw = globals().get(DiscontinuumConstitutiveLawString)()
            DiscontinuumConstitutiveLaw.SetConstitutiveLawInProperties(properties)             

    # constructing the inlet and intializing it (must be done AFTER the dem_model_part Initialize)    
    DEM_inlet = DEM_Inlet(DEM_inlet_model_part)    
    DEM_inlet.InitializeDEM_Inlet(dem_model_part, creator_destructor, dem_inlet_element_type)
  
#------------------------------------------DEM_PROCEDURES FUNCTIONS & INITIALIZATIONS--------------------------------------------------------
#if (DEM_parameters.PredefinedSkinOption == "ON" ):
   #ProceduresSetPredefinedSkin(dem_model_part)

DEMFEMProcedures = DEM_procedures.DEMFEMProcedures(DEM_parameters, graphs_path, dem_model_part, fem_model_part)

#Procedures.SetCustomSkin(dem_model_part)

materialTest.Initialize(DEM_parameters, procedures, solver, graphs_path, post_path, dem_model_part, fem_model_part)

KRATOSprint("Initialization Complete" + "\n")

step           = 0
time           = 0.0
time_old_print = 0.0

report.Prepare(timer,DEM_parameters.ControlTime)

first_print  = True; index_5 = 1; index_10  = 1; index_50  = 1; control = 0.0
    
## MODEL DATA #~CHARLIE~#:????

if (DEM_parameters.ModelDataInfo == "ON"):
    os.chdir(data_and_results)
    if (DEM_parameters.ContactMeshOption == "ON"):
      (coordination_number) = procedures.ModelData(dem_model_part, contact_model_part, solver)       # calculates the mean number of neighbours the mean radius, etc..
      KRATOSprint ("Coordination Number: " + str(coordination_number) + "\n")
      os.chdir(main_path)
    else:
      KRATOSprint("Activate Contact Mesh for ModelData information")

if(DEM_parameters.Dempack):
#    if(mpi.rank == 0):
    #materialTest.PrintChart();
    materialTest.PrepareDataForGraph()
 
##############################################################################
#                                                                            #
#    MAIN LOOP                                                               #
#                                                                            #
##############################################################################

dt = dem_model_part.ProcessInfo.GetValue(DELTA_TIME) #Calulated in dem initialize
fem_model_part.ProcessInfo[DELTA_TIME] = dt

FEM_main_step_solver.Initialize() #prediction of time step in explicit schemes

dt = fem_model_part.ProcessInfo[DELTA_TIME]
dem_model_part.ProcessInfo[DELTA_TIME] = dt

report.total_steps_expected = int(DEM_parameters.FinalTime / dt)

KRATOSprint(report.BeginReport(timer))

mesh_motion = DEMFEMUtilities()

# creating a Post Utils object that executes several post-related tasks
post_utils = DEM_procedures.PostUtils(DEM_parameters, dem_model_part)

a = 50
step = 0

start_steps = fem_buffer_size - 1
start_time = start_steps * dt;

while ( time < DEM_parameters.FinalTime):

    dt   = dem_model_part.ProcessInfo.GetValue(DELTA_TIME) # Possible modifications of DELTA_TIME
    time = time + dt

    dem_model_part.ProcessInfo[TIME]       = time
    dem_model_part.ProcessInfo[DELTA_TIME] = dt
    dem_model_part.ProcessInfo[TIME_STEPS] = step

    fem_model_part.ProcessInfo[TIME]       = time
    fem_model_part.ProcessInfo[DELTA_TIME] = dt
    fem_model_part.ProcessInfo[TIME_STEPS] = step

    #print("STEP:"+str(step)+"*******************************************")
    #print("*************************************************************")

    # Perform a partition to balance the problem
    #if(not(step%(a-1))):
        #parallelutils.Repart(dem_model_part)
        #parallelutils.CalculateModelNewIds(dem_model_part)
    
    #walls movement:
    mesh_motion.MoveAllMeshes(fem_model_part, time)
    
    #### SOLVE #########################################
    solver.Solve()
    
    if(step > start_steps):
      
      # solve time step non-linear system

      FEM_main_step_solver.Solve()
      
      # plot graphs

      # update previous time step
      fem_model_part.ProcessInfo[PREVIOUS_DELTA_TIME] = dt;

      # incremental load
      incr_steps = 1
      conditions.SetIncrementalLoad(incr_steps, dt);
    #### TIME CONTROL ##################################
    
    # adding DEM elements by the inlet:
    if (inlet_option):
        DEM_inlet.CreateElementsFromInletMesh(dem_model_part, DEM_inlet_model_part, creator_destructor, dem_inlet_element_type)  # After solving, to make sure that neighbours are already set.              

    stepinfo = report.StepiReport(timer,time,step)
    if stepinfo:
        KRATOSprint(stepinfo)
    
    #### PRINTING GRAPHS ####
    os.chdir(graphs_path)
    # measuring mean velocities in a certain control volume (the 'velocity trap')
    if (DEM_parameters.VelocityTrapOption):
        post_utils.ComputeMeanVelocitiesinTrap("Average_Velocity", time)

    #### MATERIAL TEST GRAPHS ############################
    materialTest.MeasureForcesAndPressure()
    materialTest.PrintGraph(step)

    #### GENERAL FORCE GRAPHS ###################################
    DEMFEMProcedures.MeasureForces()
    DEMFEMProcedures.PrintGraph(time)

    #### GiD IO ##########################################
    time_to_print = time - time_old_print

    if ( time_to_print >= DEM_parameters.OutputTimeStep):
        
        
        KRATOSprint("*******************  PRINTING RESULTS FOR GID  ***************************")
        KRATOSprint("                        ("+ str(dem_model_part.NumberOfElements(0)) + " elements)")
        KRATOSprint("")
        sys.stdout.flush()

        os.chdir(data_and_results)

        #properties_list = ProceduresMonitorPhysicalProperties(dem_model_part, physics_calculator, properties_list)

        os.chdir(list_path)
        demio.PrintMultifileLists(time)
        os.chdir(main_path)

        os.chdir(post_path)

        demio.PrintResults(mixed_model_part,dem_model_part,fem_model_part,contact_model_part, time)
        
        if (DEM_parameters.ContactMeshOption == "ON"):
            solver.PrepareContactElementsForPrinting()
            demio.PrintingContactElementsVariables(contact_model_part, time)

        os.chdir(main_path)

        time_old_print = time
  
    step += 1

    #if((step%500) == 0):
      #if (( DEM_parameters.ContactMeshOption =="ON") and (DEM_parameters.TestType!= "None"))  :
          #MaterialTest.OrientationStudy(contact_model_part, step)
    

##############################################################################
#                                                                            #
#    FINALIZATION                                                            #
#                                                                            #
##############################################################################

demio.FinalizeMesh()
materialTest.FinalizeGraphs()
DEMFEMProcedures.FinalizeGraphs()

# Charlie: This didn't exist. I replaced it with the line above
#if((DEM_parameters.TestType == "None")):
#    Procedures.FinalizeGraphs()

demio.CloseMultifiles()

os.chdir(main_path)

# Print tmes and more info
KRATOSprint(report.FinalReport(timer))

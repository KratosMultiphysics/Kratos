from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Python imports
import time as timer
import os
import sys
import math

# Kratos
from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *

# DEM Application
import DEM_explicit_solver_var as DEM_parameters

# TO_DO: Ungly fix. Change it. I don't like this to be in the main...
# Strategy object
if   (DEM_parameters.ElementType == "SphericPartDEMElement3D"     or DEM_parameters.ElementType == "CylinderPartDEMElement2D"):
    import sphere_strategy as SolverStrategy
elif (DEM_parameters.ElementType == "SphericContPartDEMElement3D" or DEM_parameters.ElementType == "CylinderContPartDEMElement2D"):
    import continuum_sphere_strategy as SolverStrategy

# Import MPI modules if needed. This way to do this is only valid when using OpenMPI. For other implementations of MPI it will not work.
if "OMPI_COMM_WORLD_SIZE" in os.environ:
    # Kratos MPI
    from KratosMultiphysics.MetisApplication import *
    from KratosMultiphysics.MPISearchApplication import *
    from KratosMultiphysics.mpi import *

    # DEM Application MPI
    import DEM_procedures_mpi as DEM_procedures
    import DEM_material_test_script_mpi as DEM_material_test_script
else:
    # DEM Application
    import DEM_procedures
    import DEM_material_test_script

    print("Running under OpenMP")

##############################################################################
#                                                                            #
#    INITIALIZE                                                              #
#                                                                            #
##############################################################################

# Import utilities from models
procedures    = DEM_procedures.Procedures(DEM_parameters)
demio         = DEM_procedures.DEMIo()
report        = DEM_procedures.Report()
parallelutils = DEM_procedures.ParallelUtils()
materialTest  = DEM_procedures.MaterialTest()
 
# Set the print function TO_DO: do this better...
KRATOSprint   = procedures.KRATOSprint

# Preprocess the model
procedures.PreProcessModel(DEM_parameters)

# Prepare modelparts
balls_model_part      = ModelPart("SpheresPart");
rigid_face_model_part = ModelPart("RigidFace_Part");  
mixed_model_part      = ModelPart("Mixed_Part");
cluster_model_part    = ModelPart("Cluster_Part");
DEM_inlet_model_part  = ModelPart("DEMInletPart")
contact_model_part    = ""

# Add variables
procedures.AddCommonVariables(balls_model_part, DEM_parameters)
procedures.AddBallsVariables(balls_model_part, DEM_parameters)
procedures.AddMpiVariables(balls_model_part)
SolverStrategy.AddAdditionalVariables(balls_model_part, DEM_parameters)
procedures.AddCommonVariables(rigid_face_model_part, DEM_parameters)
procedures.AddFEMVariables(rigid_face_model_part, DEM_parameters)
procedures.AddMpiVariables(rigid_face_model_part)
procedures.AddCommonVariables(cluster_model_part, DEM_parameters)
procedures.AddClusterVariables(cluster_model_part, DEM_parameters)
procedures.AddMpiVariables(cluster_model_part)
procedures.AddCommonVariables(DEM_inlet_model_part, DEM_parameters)
procedures.AddBallsVariables(DEM_inlet_model_part, DEM_parameters)
SolverStrategy.AddAdditionalVariables(DEM_inlet_model_part, DEM_parameters)  


# Reading the model_part
spheres_mp_filename   = DEM_parameters.problem_name + "DEM"
model_part_io_spheres = ModelPartIO(spheres_mp_filename,True)

# Perform the initial partition
[model_part_io_spheres, balls_model_part, MPICommSetup] = parallelutils.PerformInitialPartition(balls_model_part, model_part_io_spheres, spheres_mp_filename)

model_part_io_spheres.ReadModelPart(balls_model_part)

rigidFace_mp_filename = DEM_parameters.problem_name + "DEM_FEM_boundary"
model_part_io_fem = ModelPartIO(rigidFace_mp_filename)
model_part_io_fem.ReadModelPart(rigid_face_model_part)

clusters_mp_filename = DEM_parameters.problem_name + "DEM_Clusters"
model_part_io_clusters = ModelPartIO(clusters_mp_filename)
model_part_io_clusters.ReadModelPart(cluster_model_part)

DEM_Inlet_filename = DEM_parameters.problem_name + "DEM_Inlet"  
model_part_io_demInlet = ModelPartIO(DEM_Inlet_filename)
model_part_io_demInlet.ReadModelPart(DEM_inlet_model_part)


# Setting up the buffer size
balls_model_part.SetBufferSize(1)
cluster_model_part.SetBufferSize(1)
DEM_inlet_model_part.SetBufferSize(1)
#TODO: what about the buffersize of the rigid_face_model_part??

# Adding dofs
SolverStrategy.AddDofs(balls_model_part)
SolverStrategy.AddDofs(cluster_model_part)
SolverStrategy.AddDofs(DEM_inlet_model_part)

# Constructing a creator/destructor object
creator_destructor = ParticleCreatorDestructor()

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

os.chdir(main_path)

KRATOSprint("Initializing Problem....")

# Initialize GiD-IO
demio.AddGlobalVariables()
demio.AddBallVariables()
demio.AddFEMBoundaryVariables()
demio.AddClusterVariables()
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
                     balls_model_part,
                     rigid_face_model_part,
                     cluster_model_part,
                     contact_model_part)

os.chdir(post_path)

# Perform a partition to balance the problem
parallelutils.Repart(balls_model_part)
parallelutils.CalculateModelNewIds(balls_model_part)

os.chdir(post_path)

#Setting up the BoundingBox
if(DEM_parameters.BoundingBoxOption == "ON"):
    procedures.SetBoundingBox(balls_model_part, cluster_model_part, rigid_face_model_part, creator_destructor)

# Creating a solver object and set the search strategy
solver                 = SolverStrategy.ExplicitStrategy(balls_model_part, rigid_face_model_part, cluster_model_part, DEM_inlet_model_part, creator_destructor, DEM_parameters)
solver.search_strategy = parallelutils.GetSearchStrategy(solver, balls_model_part)

solver.Initialize()

if ( DEM_parameters.ContactMeshOption =="ON" ) :
    contact_model_part = solver.contact_model_part
  
# constructing a model part for the DEM inlet. it contains the DEM elements to be released during the simulation  
# Initializing the DEM solver must be done before creating the DEM Inlet, because the Inlet configures itself according to some options of the DEM model part

if (DEM_parameters.dem_inlet_option):    
    max_node_Id = creator_destructor.FindMaxNodeIdInModelPart(balls_model_part)
    max_FEM_node_Id = creator_destructor.FindMaxNodeIdInModelPart(rigid_face_model_part)

    if ( max_FEM_node_Id > max_node_Id):
        max_node_Id = max_FEM_node_Id
    
    creator_destructor.SetMaxNodeId(max_node_Id)                            
        
    for properties in DEM_inlet_model_part.Properties:
            
            DiscontinuumConstitutiveLawString = properties[DEM_DISCONTINUUM_CONSTITUTIVE_LAW_NAME];
            DiscontinuumConstitutiveLaw = globals().get(DiscontinuumConstitutiveLawString)()
            DiscontinuumConstitutiveLaw.SetConstitutiveLawInProperties(properties)             

    # constructing the inlet and intializing it (must be done AFTER the balls_model_part Initialize)    
    DEM_inlet = DEM_Inlet(DEM_inlet_model_part)    
    DEM_inlet.InitializeDEM_Inlet(balls_model_part, creator_destructor, DEM_parameters.dem_inlet_element_type)
  
#------------------------------------------DEM_PROCEDURES FUNCTIONS & INITIALIZATIONS--------------------------------------------------------
#if (DEM_parameters.PredefinedSkinOption == "ON" ):
   #ProceduresSetPredefinedSkin(balls_model_part)

DEMFEMProcedures = DEM_procedures.DEMFEMProcedures(DEM_parameters, graphs_path, balls_model_part, rigid_face_model_part)

#Procedures.SetCustomSkin(balls_model_part)

materialTest.Initialize(DEM_parameters, procedures, solver, graphs_path, post_path, balls_model_part, rigid_face_model_part)

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
      (coordination_number) = procedures.ModelData(balls_model_part, contact_model_part, solver)       # calculates the mean number of neighbours the mean radius, etc..
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

dt = balls_model_part.ProcessInfo.GetValue(DELTA_TIME) # Possible modifications of DELTA_TIME

report.total_steps_expected = int(DEM_parameters.FinalTime / dt)

KRATOSprint(report.BeginReport(timer))

mesh_motion = DEMFEMUtilities()

# creating a Post Utils object that executes several post-related tasks
post_utils = DEM_procedures.PostUtils(DEM_parameters, balls_model_part)

a = 50
step = 0  
while ( time < DEM_parameters.FinalTime):
    #print("TIME STEP BEGINS. STEP:"+str(step)+"++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
    dt   = balls_model_part.ProcessInfo.GetValue(DELTA_TIME) # Possible modifications of DELTA_TIME
    time = time + dt

    balls_model_part.ProcessInfo[TIME]            = time
    balls_model_part.ProcessInfo[DELTA_TIME]      = dt
    balls_model_part.ProcessInfo[TIME_STEPS]      = step
    
    rigid_face_model_part.ProcessInfo[TIME]       = time
    rigid_face_model_part.ProcessInfo[DELTA_TIME] = dt
    rigid_face_model_part.ProcessInfo[TIME_STEPS] = step

    cluster_model_part.ProcessInfo[TIME]            = time
    cluster_model_part.ProcessInfo[DELTA_TIME]      = dt
    cluster_model_part.ProcessInfo[TIME_STEPS]      = step

    #print("STEP:"+str(step)+"*******************************************")
    #print("*************************************************************")

    # Perform a partition to balance the problem
    #if(not(step%(a-1))):
        #parallelutils.Repart(balls_model_part)
        #parallelutils.CalculateModelNewIds(balls_model_part)
    
    #walls movement:
    mesh_motion.MoveAllMeshes(rigid_face_model_part, time)
    
    #### SOLVE #########################################
    solver.Solve()
    
    #### TIME CONTROL ##################################
    
    # adding DEM elements by the inlet:
    if (DEM_parameters.dem_inlet_option):
        DEM_inlet.CreateElementsFromInletMesh(balls_model_part, DEM_inlet_model_part, creator_destructor, DEM_parameters.dem_inlet_element_type)  # After solving, to make sure that neighbours are already set.              

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
        KRATOSprint("                        ("+ str(balls_model_part.NumberOfElements(0)) + " elements)")
        KRATOSprint("                        ("+ str(balls_model_part.NumberOfNodes(0)) + " nodes)")
        KRATOSprint("")
        sys.stdout.flush()

        os.chdir(data_and_results)

        #properties_list = ProceduresMonitorPhysicalProperties(balls_model_part, physics_calculator, properties_list)

        os.chdir(list_path)
        demio.PrintMultifileLists(time,post_path)
        os.chdir(main_path)

        os.chdir(post_path)

        if (DEM_parameters.ContactMeshOption == "ON"):
            solver.PrepareContactElementsForPrinting()
        
        demio.PrintResults(mixed_model_part, balls_model_part, rigid_face_model_part, cluster_model_part, contact_model_part, time)
                
        os.chdir(main_path)

        time_old_print = time
  
    step += 1

    #if((step%500) == 0):
      #if (( DEM_parameters.ContactMeshOption =="ON") and (DEM_parameters.TestType!= "None"))  :
          #MaterialTest.OrientationStudy(contact_model_part, step)
    
    #print("TIME STEP ENDS ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
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

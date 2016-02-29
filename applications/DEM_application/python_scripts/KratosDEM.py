from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Python imports
import time as timer
import os
import sys
#import math

# Kratos
from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *

sys.path.insert(0,'')
# DEM Application
import DEM_explicit_solver_var as DEM_parameters

# Import MPI modules if needed. This way to do this is only valid when using OpenMPI. For other implementations of MPI it will not work.
if "OMPI_COMM_WORLD_SIZE" in os.environ:
    # Kratos MPI
    from KratosMultiphysics.MetisApplication import *
    from KratosMultiphysics.MPISearchApplication import *
    from KratosMultiphysics.mpi import *

    # DEM Application MPI
    import DEM_procedures_mpi as DEM_procedures
    import DEM_material_test_script_mpi as DEM_material_test_script
    print("Running under MPI...........")
else:
    # DEM Application
    import DEM_procedures
    import DEM_material_test_script
    print("Running under OpenMP........")

# TO_DO: Ugly fix. Change it. I don't like this to be in the main...
# Strategy object
if   (DEM_parameters.ElementType == "SphericPartDEMElement3D"     or DEM_parameters.ElementType == "CylinderPartDEMElement2D"):
    import sphere_strategy as SolverStrategy
elif (DEM_parameters.ElementType == "SphericContPartDEMElement3D" or DEM_parameters.ElementType == "CylinderContPartDEMElement2D"):
    import continuum_sphere_strategy as SolverStrategy
elif (DEM_parameters.ElementType == "ThermalSphericContPartDEMElement3D"):
    import thermal_continuum_sphere_strategy as SolverStrategy    
else:
    KRATOSprint('Error: Strategy unnavailable. Select a different scheme-element')

##############################################################################
#                                                                            #
#    INITIALIZE                                                              #
#                                                                            #
##############################################################################

# Import utilities from models
procedures    = DEM_procedures.Procedures(DEM_parameters)
demio         = DEM_procedures.DEMIo(DEM_parameters)
report        = DEM_procedures.Report()
parallelutils = DEM_procedures.ParallelUtils()
materialTest  = DEM_procedures.MaterialTest()
 
# Set the print function TO_DO: do this better...
KRATOSprint   = procedures.KRATOSprint

# Preprocess the model
procedures.PreProcessModel(DEM_parameters)

# Prepare modelparts
spheres_model_part    = ModelPart("SpheresPart")
rigid_face_model_part = ModelPart("RigidFace_Part")
mixed_model_part      = ModelPart("Mixed_Part")
cluster_model_part    = ModelPart("Cluster_Part")
DEM_inlet_model_part  = ModelPart("DEMInletPart")
mapping_model_part    = ModelPart("Mappingmodel_part")
contact_model_part    = ""

# Constructing a utilities objects
creator_destructor = ParticleCreatorDestructor()
dem_fem_search = DEM_FEM_Search()

# Creating a solver object and set the search strategy

solver = SolverStrategy.ExplicitStrategy(spheres_model_part, rigid_face_model_part, cluster_model_part, DEM_inlet_model_part, creator_destructor, dem_fem_search, DEM_parameters)

#Getting chosen scheme:

if (DEM_parameters.IntegrationScheme == 'Forward_Euler'):
    scheme = ForwardEulerScheme()
elif (DEM_parameters.IntegrationScheme == 'Symplectic_Euler'):
    scheme = SymplecticEulerScheme()
elif (DEM_parameters.IntegrationScheme == 'Taylor_Scheme'):
    scheme = TaylorScheme()
elif (DEM_parameters.IntegrationScheme == 'Newmark_Beta_Method'):
    scheme = NewmarkBetaScheme(0.5, 0.25)
elif (DEM_parameters.IntegrationScheme == 'Verlet_Velocity'):
    scheme = VerletVelocityScheme()
else:
    KRATOSprint('Error: selected scheme not defined. Please select a different scheme')

scheme.SetRotationOption(solver.rotation_option)
solver.time_integration_scheme = scheme

# Add variables
procedures.AddCommonVariables(spheres_model_part, DEM_parameters)
procedures.AddSpheresVariables(spheres_model_part, DEM_parameters)
procedures.AddMpiVariables(spheres_model_part)
solver.AddAdditionalVariables(spheres_model_part, DEM_parameters)
scheme.AddSpheresVariables(spheres_model_part)
procedures.AddCommonVariables(cluster_model_part, DEM_parameters)
procedures.AddClusterVariables(cluster_model_part, DEM_parameters)
procedures.AddMpiVariables(cluster_model_part)
scheme.AddClustersVariables(cluster_model_part)
procedures.AddCommonVariables(DEM_inlet_model_part, DEM_parameters)
procedures.AddSpheresVariables(DEM_inlet_model_part, DEM_parameters)
solver.AddAdditionalVariables(DEM_inlet_model_part, DEM_parameters)  
scheme.AddSpheresVariables(DEM_inlet_model_part)
procedures.AddCommonVariables(rigid_face_model_part, DEM_parameters)
procedures.AddRigidFaceVariables(rigid_face_model_part, DEM_parameters)
procedures.AddMpiVariables(rigid_face_model_part)

# Reading the model_part
spheres_mp_filename   = DEM_parameters.problem_name + "DEM"
model_part_io_spheres = ModelPartIO(spheres_mp_filename)

# Perform the initial partition
[model_part_io_spheres, spheres_model_part, MPICommSetup] = parallelutils.PerformInitialPartition(spheres_model_part, model_part_io_spheres, spheres_mp_filename)

model_part_io_spheres.ReadModelPart(spheres_model_part)

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
spheres_model_part.SetBufferSize(1)
cluster_model_part.SetBufferSize(1)
DEM_inlet_model_part.SetBufferSize(1)
rigid_face_model_part.SetBufferSize(1)

# Adding dofs
solver.AddDofs(spheres_model_part)
solver.AddDofs(cluster_model_part)
solver.AddDofs(DEM_inlet_model_part)

# Creating necessary directories
main_path = os.getcwd()
[post_path,data_and_results,graphs_path,MPI_results] = procedures.CreateDirectories(str(main_path),str(DEM_parameters.problem_name))

os.chdir(main_path)

KRATOSprint("\nInitializing Problem...")

# Initialize GiD-IO
demio.AddGlobalVariables()
demio.AddSpheresVariables()
demio.AddFEMBoundaryVariables()
demio.AddClusterVariables()
demio.AddContactVariables()
#
demio.AddMpiVariables()
demio.EnableMpiVariables()

demio.Configure(DEM_parameters.problem_name,
                DEM_parameters.OutputFileType,
                DEM_parameters.Multifile,
                DEM_parameters.ContactMeshOption)

demio.SetOutputName(DEM_parameters.problem_name)

os.chdir(post_path)

multifiles = (
    DEM_procedures.MultifileList(DEM_parameters.problem_name,1 ),
    DEM_procedures.MultifileList(DEM_parameters.problem_name, 2),
    DEM_procedures.MultifileList(DEM_parameters.problem_name,5 ),
    DEM_procedures.MultifileList(DEM_parameters.problem_name,10),
    DEM_procedures.MultifileList(DEM_parameters.problem_name,20),
    DEM_procedures.MultifileList(DEM_parameters.problem_name,50),
    )

demio.SetMultifileLists(multifiles)
os.chdir(post_path)

demio.InitializeMesh(mixed_model_part,
                     spheres_model_part,
                     rigid_face_model_part,
                     cluster_model_part,
                     contact_model_part,
                     mapping_model_part)

os.chdir(post_path)

# Perform a partition to balance the problem
parallelutils.Repart(spheres_model_part)
parallelutils.CalculateModelNewIds(spheres_model_part)

os.chdir(post_path)

#Setting up the BoundingBox
bounding_box_time_limits = []
if (DEM_parameters.BoundingBoxOption == "ON"):
    procedures.SetBoundingBox(spheres_model_part, cluster_model_part, rigid_face_model_part, creator_destructor)
    bounding_box_time_limits = [solver.bounding_box_start_time, solver.bounding_box_stop_time]

# Creating a solver object and set the search strategy
#solver                 = SolverStrategy.ExplicitStrategy(spheres_model_part, rigid_face_model_part, cluster_model_part, DEM_inlet_model_part, creator_destructor, DEM_parameters)
solver.search_strategy = parallelutils.GetSearchStrategy(solver, spheres_model_part)

dt = DEM_parameters.MaxTimeStep

#Finding the max id of the nodes... (it is necessary for anything that will add spheres to the spheres_model_part, for instance, the INLETS and the CLUSTERS read from mdpa file.
max_node_Id = creator_destructor.FindMaxNodeIdInModelPart(spheres_model_part)
max_FEM_node_Id = creator_destructor.FindMaxNodeIdInModelPart(rigid_face_model_part)

if (max_FEM_node_Id > max_node_Id):
    max_node_Id = max_FEM_node_Id

creator_destructor.SetMaxNodeId(max_node_Id)    

#Strategy Initialization
solver.Initialize()    # Possible modifications of DELTA_TIME, number of elements and number of nodes.

if (DEM_parameters.ContactMeshOption =="ON"):
    contact_model_part = solver.contact_model_part
  
# constructing a model part for the DEM inlet. it contains the DEM elements to be released during the simulation  
# Initializing the DEM solver must be done before creating the DEM Inlet, because the Inlet configures itself according to some options of the DEM model part
if (DEM_parameters.dem_inlet_option):                                        
    # constructing the inlet and initializing it (must be done AFTER the spheres_model_part Initialize)    
    DEM_inlet = DEM_Inlet(DEM_inlet_model_part)    
    DEM_inlet.InitializeDEM_Inlet(spheres_model_part, creator_destructor)
  
#------------------------------------------DEM_PROCEDURES FUNCTIONS & INITIALIZATIONS--------------------------------------------------------
#if (DEM_parameters.PredefinedSkinOption == "ON" ):
   #ProceduresSetPredefinedSkin(spheres_model_part)

DEMFEMProcedures = DEM_procedures.DEMFEMProcedures(DEM_parameters, graphs_path, spheres_model_part, rigid_face_model_part)

#Procedures.SetCustomSkin(spheres_model_part)

materialTest.Initialize(DEM_parameters, procedures, solver, graphs_path, post_path, spheres_model_part, rigid_face_model_part)

KRATOSprint("Initialization Complete" + "\n")

step           = 0
time           = 0.0
time_old_print = 0.0

report.Prepare(timer, DEM_parameters.ControlTime)

first_print = True; index_5 = 1; index_10 = 1; index_50 = 1; control = 0.0
    
## MODEL DATA #~CHARLIE~#:????

if (DEM_parameters.ModelDataInfo == "ON"):
    os.chdir(data_and_results)
    if (DEM_parameters.ContactMeshOption == "ON"):
        (coordination_number) = procedures.ModelData(spheres_model_part, contact_model_part, solver) # Calculates the mean number of neighbours the mean radius, etc..
        KRATOSprint ("Coordination Number: " + str(coordination_number) + "\n")
        os.chdir(main_path)
    else:
        KRATOSprint("Activate Contact Mesh for ModelData information")

materialTest.PrintChart()
materialTest.PrepareDataForGraph()
 
##############################################################################
#                                                                            #
#    MAIN LOOP                                                               #
#                                                                            #
##############################################################################

report.total_steps_expected = int(DEM_parameters.FinalTime / dt)

KRATOSprint(report.BeginReport(timer))

mesh_motion = DEMFEMUtilities()

# creating a Post Utils object that executes several post-related tasks
post_utils = DEM_procedures.PostUtils(DEM_parameters, spheres_model_part)

step = 0

import plot_variables #Related to debugging

nodeplotter = 0 #Related to debugging

list_of_nodes_ids = []

for node in spheres_model_part.Nodes:
    list_of_nodes_ids.append(node.Id)

if nodeplotter: #Related to debugging
    os.chdir(main_path)
    plotter = plot_variables.variable_plotter(spheres_model_part, list_of_nodes_ids) #Related to debugging

while (time < DEM_parameters.FinalTime):
    
    dt    = spheres_model_part.ProcessInfo.GetValue(DELTA_TIME) # Possible modifications of DELTA_TIME
    time  = time + dt
    step += 1

    spheres_model_part.ProcessInfo[TIME]          = time
    spheres_model_part.ProcessInfo[DELTA_TIME]    = dt
    spheres_model_part.ProcessInfo[TIME_STEPS]    = step
    
    rigid_face_model_part.ProcessInfo[TIME]       = time
    rigid_face_model_part.ProcessInfo[DELTA_TIME] = dt
    rigid_face_model_part.ProcessInfo[TIME_STEPS] = step

    cluster_model_part.ProcessInfo[TIME]          = time
    cluster_model_part.ProcessInfo[DELTA_TIME]    = dt
    cluster_model_part.ProcessInfo[TIME_STEPS]    = step

    # Perform a partition to balance the problem
    #if(not(step%(a-1))):
        #parallelutils.Repart(spheres_model_part)
        #parallelutils.CalculateModelNewIds(spheres_model_part)
    
    #### MATERIAL TEST LOAD&UNLOAD  ####################
    #materialTest.ApplyMovementbySteps(time)
    
    #### GENERAL TEST LOAD&UNLOAD  #####################
    #DEMFEMProcedures.ApplyMovementbySteps(time)
    
    #### SOLVE #########################################
    solver.Solve()
    
    # Walls movement:
    mesh_motion.MoveAllMeshes(rigid_face_model_part, time, dt)
    mesh_motion.MoveAllMeshes(spheres_model_part, time, dt)
    mesh_motion.MoveAllMeshes(DEM_inlet_model_part, time, dt)
    
    #### TIME CONTROL ##################################
    
    # adding DEM elements by the inlet:
    if (DEM_parameters.dem_inlet_option):
        DEM_inlet.CreateElementsFromInletMesh(spheres_model_part, cluster_model_part, creator_destructor)  # After solving, to make sure that neighbours are already set.              

    stepinfo = report.StepiReport(timer,time,step)
    if stepinfo:
        KRATOSprint(stepinfo)
    
    #### PRINTING GRAPHS ####
    os.chdir(graphs_path)
    # measuring mean velocities in a certain control volume (the 'velocity trap')
    if (DEM_parameters.VelocityTrapOption):
        compute_flow = False
        post_utils.ComputeMeanVelocitiesinTrap("Average_Velocity.txt", time, compute_flow)

    #### MATERIAL TEST GRAPHS ############################
    materialTest.MeasureForcesAndPressure()
    materialTest.PrintGraph(time)
    
    #### GENERAL FORCE GRAPHS ############################
    #DEMFEMProcedures.MeasureForces()
    DEMFEMProcedures.PrintGraph(time)
    DEMFEMProcedures.PrintBallsGraph(time)

    #### GiD IO ##########################################
    time_to_print = time - time_old_print

    if (DEM_parameters.OutputTimeStep - time_to_print < 1e-2 * dt):

        KRATOSprint("*******************  PRINTING RESULTS FOR GID  ***************************")
        KRATOSprint("                        ("+ str(spheres_model_part.NumberOfElements(0)) + " elements)")
        KRATOSprint("                        ("+ str(spheres_model_part.NumberOfNodes(0)) + " nodes)")
        KRATOSprint("")
        sys.stdout.flush()

        os.chdir(data_and_results)

        #properties_list = ProceduresMonitorPhysicalProperties(spheres_model_part, physics_calculator, properties_list)

        demio.PrintMultifileLists(time, post_path)
        #os.chdir(main_path)

        os.chdir(post_path)

        solver.PrepareElementsForPrinting()
        if (DEM_parameters.ContactMeshOption == "ON"):
            solver.PrepareContactElementsForPrinting()
        
        demio.PrintResults(mixed_model_part, spheres_model_part, rigid_face_model_part, cluster_model_part, contact_model_part, mapping_model_part, creator_destructor, dem_fem_search, time, bounding_box_time_limits)
        os.chdir(main_path)
        
        if nodeplotter:
            plotter.plot_variables(time) #Related to debugging

        time_old_print = time

    #if((step%500) == 0):
      #if (( DEM_parameters.ContactMeshOption =="ON") and (DEM_parameters.TestType!= "None"))  :
          #MaterialTest.OrientationStudy(contact_model_part, step)
    
    #print("TIME STEP ENDS +++++++++++++++++++++++++++++++++++++++++++++++++")
##############################################################################
#                                                                            #
#    FINALIZATION                                                            #
#                                                                            #
##############################################################################

demio.FinalizeMesh()
materialTest.FinalizeGraphs()
DEMFEMProcedures.FinalizeGraphs(rigid_face_model_part)
DEMFEMProcedures.FinalizeBallsGraphs(spheres_model_part)

# Charlie: This didn't exist. I replaced it with the line above
#if((DEM_parameters.TestType == "None")):
#    Procedures.FinalizeGraphs()

demio.CloseMultifiles()

os.chdir(main_path)

# Print times and more info
KRATOSprint(report.FinalReport(timer))

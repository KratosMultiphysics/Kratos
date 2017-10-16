from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Python imports
import time as timer
import os
import sys
import math

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

#EXTRA IMPORTS

# TO_DO: Ungly fix. Change it. I don't like this to be in the main...
# Strategy object
if   (DEM_parameters["ElementType"].GetString() == "SphericPartDEMElement3D"     or DEM_parameters["ElementType"].GetString() == "CylinderPartDEMElement2D"):
    import sphere_strategy as SolverStrategy
elif (DEM_parameters["ElementType"].GetString() == "SphericContPartDEMElement3D" or DEM_parameters["ElementType"].GetString() == "CylinderContPartDEMElement2D"):
    import continuum_sphere_strategy as SolverStrategy
elif (DEM_parameters["ElementType"].GetString() == "ThermalSphericContPartDEMElement3D"):
    import thermal_continuum_sphere_strategy as SolverStrategy    

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
spheres_model_part    = ModelPart("SpheresPart");
rigid_face_model_part = ModelPart("RigidFace_Part");  
mixed_model_part      = ModelPart("Mixed_Part");
cluster_model_part    = ModelPart("Cluster_Part");
DEM_inlet_model_part  = ModelPart("DEMInletPart")
mapping_model_part    = ModelPart("Mappingmodel_part")
contact_model_part    = ""

#EXTRA ModelPart Operations
#
#
#

# Constructing a creator/destructor object
creator_destructor = ParticleCreatorDestructor()
#
#

# Creating a solver object and set the search strategy
solver                 = SolverStrategy.ExplicitStrategy(spheres_model_part, rigid_face_model_part, cluster_model_part, DEM_inlet_model_part, creator_destructor, DEM_parameters)

# Add variables
procedures.AddCommonVariables(spheres_model_part, DEM_parameters)
procedures.AddSpheresVariables(spheres_model_part, DEM_parameters)
procedures.AddMpiVariables(spheres_model_part)
solver.AddAdditionalVariables(spheres_model_part, DEM_parameters)
procedures.AddCommonVariables(cluster_model_part, DEM_parameters)
procedures.AddClusterVariables(cluster_model_part, DEM_parameters)
procedures.AddMpiVariables(cluster_model_part)
procedures.AddCommonVariables(DEM_inlet_model_part, DEM_parameters)
procedures.AddSpheresVariables(DEM_inlet_model_part, DEM_parameters)
solver.AddAdditionalVariables(DEM_inlet_model_part, DEM_parameters)  
procedures.AddCommonVariables(rigid_face_model_part, DEM_parameters)
procedures.AddRigidFaceVariables(rigid_face_model_part, DEM_parameters)
procedures.AddMpiVariables(rigid_face_model_part)

# Reading the model_part
spheres_mp_filename   = DEM_parameters["problem_name"].GetString() + "DEM"
model_part_io_spheres = ModelPartIO(spheres_mp_filename)

if "do_not_perform_initial_partition" in DEM_parameters.keys() and DEM_parameters["do_not_perform_initial_partition"].GetBool():
    pass
else:
    parallelutils.PerformInitialPartition(model_part_io_spheres)

[model_part_io_spheres, spheres_model_part, MPICommSetup] = parallelutils.SetCommunicator(spheres_model_part, model_part_io_spheres, spheres_mp_filename)

model_part_io_spheres.ReadModelPart(spheres_model_part)

rigidFace_mp_filename = DEM_parameters["problem_name"].GetString() + "DEM_FEM_boundary"
model_part_io_fem = ModelPartIO(rigidFace_mp_filename)
model_part_io_fem.ReadModelPart(rigid_face_model_part)

clusters_mp_filename = DEM_parameters["problem_name"].GetString() + "DEM_Clusters"
model_part_io_clusters = ModelPartIO(clusters_mp_filename)
model_part_io_clusters.ReadModelPart(cluster_model_part)

DEM_Inlet_filename = DEM_parameters["problem_name"].GetString() + "DEM_Inlet"  
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
#

#Utilities


# Creating necessary directories
main_path = os.getcwd()
[post_path,list_path,data_and_results,graphs_path,MPI_results] = procedures.CreateDirectories(str(main_path),str(DEM_parameters["problem_name"].GetString()))

os.chdir(main_path)

KRATOSprint("Initializing Problem....")

# Initialize GiD-IO
demio.AddGlobalVariables()
demio.AddSpheresVariables()
demio.AddFEMBoundaryVariables()
demio.AddClusterVariables()
demio.AddContactVariables()
#
demio.AddMpiVariables()

demio.Configure(DEM_parameters["problem_name"].GetString(),
                DEM_parameters["OutputFileType"].GetString(),
                DEM_parameters["Multifile"].GetString(),
                DEM_parameters["ContactMeshOption"].GetBool())

demio.SetOutputName(DEM_parameters["problem_name"].GetString())

os.chdir(post_path)

multifiles = (
    DEM_procedures.MultifileList(DEM_parameters["problem_name"].GetString(),1 ),
    DEM_procedures.MultifileList(DEM_parameters["problem_name"].GetString(), 2),
    DEM_procedures.MultifileList(DEM_parameters["problem_name"].GetString(),5 ),
    DEM_procedures.MultifileList(DEM_parameters["problem_name"].GetString(),10),
    DEM_procedures.MultifileList(DEM_parameters["problem_name"].GetString(),20),
    DEM_procedures.MultifileList(DEM_parameters["problem_name"].GetString(),50),
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
if DEM_parameters["BoundingBoxOption"].GetBool():
    procedures.SetBoundingBox(spheres_model_part, cluster_model_part, rigid_face_model_part, creator_destructor)

# Creating a solver object and set the search strategy
#solver                 = SolverStrategy.ExplicitStrategy(spheres_model_part, rigid_face_model_part, cluster_model_part, DEM_inlet_model_part, creator_destructor, DEM_parameters)
solver.search_strategy = parallelutils.GetSearchStrategy(solver, spheres_model_part)

dt = DEM_parameters["MaxTimeStep"].GetDouble()
solver.Initialize()    # Possible modifications of DELTA_TIME
##

if DEM_parameters["ContactMeshOption"].GetBool():
    contact_model_part = solver.contact_model_part
  
# constructing a model part for the DEM inlet. it contains the DEM elements to be released during the simulation  
# Initializing the DEM solver must be done before creating the DEM Inlet, because the Inlet configures itself according to some options of the DEM model part

if DEM_parameters["dem_inlet_option"].GetBool():    
    max_node_Id = creator_destructor.FindMaxNodeIdInModelPart(spheres_model_part)
    max_FEM_node_Id = creator_destructor.FindMaxNodeIdInModelPart(rigid_face_model_part)

    if ( max_FEM_node_Id > max_node_Id):
        max_node_Id = max_FEM_node_Id
    
    creator_destructor.SetMaxNodeId(max_node_Id)                            

    # constructing the inlet and intializing it (must be done AFTER the spheres_model_part Initialize)    
    DEM_inlet = DEM_Inlet(DEM_inlet_model_part)    
    DEM_inlet.InitializeDEM_Inlet(spheres_model_part, creator_destructor)

DEMFEMProcedures = DEM_procedures.DEMFEMProcedures(DEM_parameters, graphs_path, spheres_model_part, rigid_face_model_part)

#Procedures.SetCustomSkin(spheres_model_part)

materialTest.Initialize(DEM_parameters, procedures, solver, graphs_path, post_path, spheres_model_part, rigid_face_model_part)

KRATOSprint("Initialization Complete" + "\n")

step           = 0
time           = 0.0
time_old_print = 0.0

report.Prepare(timer,DEM_parameters["ControlTime"].GetDouble())

first_print  = True; index_5 = 1; index_10  = 1; index_50  = 1; control = 0.0
    
## MODEL DATA #~CHARLIE~#:????

if DEM_parameters["ModelDataInfo"].GetBool():
    os.chdir(data_and_results)
    if DEM_parameters["ContactMeshOption"].GetBool():
      (coordination_number) = procedures.ModelData(spheres_model_part, contact_model_part, solver)       # calculates the mean number of neighbours the mean radius, etc..
      KRATOSprint ("Coordination Number: " + str(coordination_number) + "\n")
      os.chdir(main_path)
    else:
      KRATOSprint("Activate Contact Mesh for ModelData information")

if DEM_parameters["Dempack"].GetBool():
#    if(mpi.rank == 0):
    materialTest.PrintChart();
    materialTest.PrepareDataForGraph()
 
##############################################################################
#                                                                            #
#    MAIN LOOP                                                               #
#                                                                            #
##############################################################################

report.total_steps_expected = int(DEM_parameters["FinalTime"].GetDouble() / dt)

KRATOSprint(report.BeginReport(timer))

mesh_motion = DEMFEMUtilities()

# creating a Post Utils object that executes several post-related tasks
post_utils = DEM_procedures.PostUtils(DEM_parameters, spheres_model_part)

step = 0

#G
import spreader
import spreader_var as spp

numer_of_particles_rate = DEM_inlet_model_part.GetProperties()[1][INLET_NUMBER_OF_PARTICLES]
start_time =  DEM_inlet_model_part.GetProperties()[1][INLET_START_TIME]
stop_time =  DEM_inlet_model_part.GetProperties()[1][INLET_STOP_TIME]
total_balls_to_be_added = numer_of_particles_rate * (stop_time - start_time)
total_balls = total_balls_to_be_added + spheres_model_part.NumberOfNodes(0) + rigid_face_model_part.NumberOfNodes(0)
maximum_expected_particle_id = int(total_balls * spp.n_balls_security_factor)

spreader_scanner = spreader.scanner(spheres_model_part, maximum_expected_particle_id, spp.outermost_disc_radius, spp.cone_angle, spp.number_of_vanes)   
        
flags = Flags()
#Z

#G
# creating mpda
creator = DEM_procedures.MdpaCreator(main_path, DEM_parameters)
#Z

while time < DEM_parameters["FinalTime"].GetDouble():
    dt   = spheres_model_part.ProcessInfo.GetValue(DELTA_TIME) # Possible modifications of DELTA_TIME
    time = time + dt
    step += 1

    spheres_model_part.ProcessInfo[TIME]            = time
    spheres_model_part.ProcessInfo[DELTA_TIME]      = dt
    spheres_model_part.ProcessInfo[TIME_STEPS]      = step
    
    rigid_face_model_part.ProcessInfo[TIME]       = time
    rigid_face_model_part.ProcessInfo[DELTA_TIME] = dt
    rigid_face_model_part.ProcessInfo[TIME_STEPS] = step

    cluster_model_part.ProcessInfo[TIME]            = time
    cluster_model_part.ProcessInfo[DELTA_TIME]      = dt
    cluster_model_part.ProcessInfo[TIME_STEPS]      = step
    
    #walls movement:
    mesh_motion.MoveAllMeshes(rigid_face_model_part, time, dt)
    #########mesh_motion.MoveSphereMeshes(spheres_model_part, time, dt)
    
    #### SOLVE #########################################
    solver.Solve()
    
#G
    if step % spp.steps_between_exit_control == 0:
        spreader_scanner.UpdateData(time)               
#Z      
    
    #### TIME CONTROL ##################################
    
    # adding DEM elements by the inlet:
    if DEM_parameters["dem_inlet_option"].GetBool(): 
        DEM_inlet.CreateElementsFromInletMesh(spheres_model_part, cluster_model_part, creator_destructor)  # After solving, to make sure that neighbours are already set.              

    stepinfo = report.StepiReport(timer,time,step)
    if stepinfo:
        KRATOSprint(stepinfo)
    
    #### PRINTING GRAPHS ####
    os.chdir(graphs_path)
    # measuring mean velocities in a certain control volume (the 'velocity trap')
    if DEM_parameters["VelocityTrapOption"].GetBool():
        compute_flow = False
        post_utils.ComputeMeanVelocitiesinTrap("Average_Velocity.txt", time)

    #### MATERIAL TEST GRAPHS ############################
    materialTest.MeasureForcesAndPressure()
    materialTest.PrintGraph(time)
    
    #### GENERAL FORCE GRAPHS ############################
    DEMFEMProcedures.PrintGraph(time)
    #DEMFEMProcedures.PrintBallsGraph(time)  

    #### GiD IO ##########################################
    time_to_print = time - time_old_print

    if ( DEM_parameters["OutputTimeStep"].GetDouble() - time_to_print < 1e-2*dt  ):            
        
        KRATOSprint("*******************  PRINTING RESULTS FOR GID  ***************************")
        KRATOSprint("                        ("+ str(spheres_model_part.NumberOfElements(0)) + " elements)")
        KRATOSprint("                        ("+ str(spheres_model_part.NumberOfNodes(0)) + " nodes)")
        KRATOSprint("")

        os.chdir(data_and_results)

        #properties_list = ProceduresMonitorPhysicalProperties(spheres_model_part, physics_calculator, properties_list)

        os.chdir(list_path)
        demio.PrintMultifileLists(time, post_path)
        os.chdir(main_path)
        os.chdir(post_path)

        if DEM_parameters["ContactMeshOption"].GetBool():
            solver.PrepareContactElementsForPrinting()
        
        demio.PrintResults(mixed_model_part, spheres_model_part, rigid_face_model_part, cluster_model_part, contact_model_part, mapping_model_part, time)
                
        os.chdir(main_path)

        time_old_print = time
        
    #if((step%500) == 0):
      #if (( DEM_parameters["ContactMeshOption"].GetBool()) and (DEM_parameters["TestType"].GetString() != "None"))  :
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
#if((DEM_parameters["TestType"].GetString()  == "None")):
#    Procedures.FinalizeGraphs()

demio.CloseMultifiles()

os.chdir(main_path)
#G
spreader_out_path = os.getcwd() + '/out_data.txt'
spreader_in_path  = os.getcwd() +  '/in_data.txt'
spreader_scanner.PrintData(spreader_in_path, spreader_out_path)
#Z
# Print times and more info
KRATOSprint(report.FinalReport(timer))

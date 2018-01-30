from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Python imports
import time as timer
import os
import sys

# Kratos
from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *
from KratosMultiphysics.SolidMechanicsApplication import *                   #DEMFEM
sys.path.insert(0,'')
# DEM Application
import DEM_explicit_solver_var as DEM_parameters

# Import MPI modules if needed. This way to do this is only valid when using OpenMPI. For other implementations of MPI it will not work.
if "OMPI_COMM_WORLD_SIZE" in os.environ or "I_MPI_INFO_NUMA_NODE_NUM" in os.environ:
    print("Running under MPI...........")
    from KratosMultiphysics.MetisApplication import *
    from KratosMultiphysics.MPISearchApplication import *
    from KratosMultiphysics.mpi import *
    import DEM_procedures_mpi as DEM_procedures
    import DEM_material_test_script_mpi as DEM_material_test_script    
    def model_part_reader(modelpart, nodeid=0, elemid=0, condid=0):
        return ReorderConsecutiveFromGivenIdsModelPartIO(modelpart, nodeid, elemid, condid)
         
else:
    print("Running under OpenMP........")
    import DEM_procedures
    import DEM_material_test_script
    #def model_part_reader(modelpart, a=0, b=0, c=0):
    #    return ModelPartIO(modelpart)
    def model_part_reader(modelpart, nodeid=0, elemid=0, condid=0):
        return ReorderConsecutiveFromGivenIdsModelPartIO(modelpart, nodeid, elemid, condid)
    
#EXTRA IMPORTS

from KratosMultiphysics.ExternalSolversApplication import *                  #DEMFEM
import ProjectParameters as FEM_general_variables                            #DEMFEM
import constitutive_law_python_utility as constitutive_law_utils             #DEMFEM
import conditions_python_utility as condition_utils                          #DEMFEM

# TODO: Ugly fix. Change it. I don't like this to be in the main...
# Strategy object
if (DEM_parameters["ElementType"].GetString() == "SphericPartDEMElement3D" or DEM_parameters["ElementType"].GetString() == "CylinderPartDEMElement2D"):
    import sphere_strategy as SolverStrategy
elif (DEM_parameters["ElementType"].GetString() == "SphericContPartDEMElement3D" or DEM_parameters["ElementType"].GetString() == "CylinderContPartDEMElement2D"):
    import continuum_sphere_strategy as SolverStrategy
elif (DEM_parameters["ElementType"].GetString() == "ThermalSphericContPartDEMElement3D"):
    import thermal_continuum_sphere_strategy as SolverStrategy
elif (DEM_parameters["ElementType"].GetString() == "ThermalSphericPartDEMElement3D"):
    import thermal_sphere_strategy as SolverStrategy  
elif (DEM_parameters["ElementType"].GetString() == "SinteringSphericConPartDEMElement3D"):
    import thermal_continuum_sphere_strategy as SolverStrategy     
else:
    KRATOSprint('Error: Strategy unavailable. Select a different scheme-element')
        
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
#procedures.PreProcessModel(DEM_parameters)

# Prepare modelparts
spheres_model_part    = ModelPart("SpheresPart")
rigid_face_model_part = ModelPart("RigidFace_Part")
cluster_model_part    = ModelPart("Cluster_Part")
DEM_inlet_model_part  = ModelPart("DEMInletPart")
mapping_model_part    = ModelPart("Mappingmodel_part")
contact_model_part    = ""

#EXTRA ModelPart Operations
FEMSolverSettings = FEM_general_variables.SolverSettings                           #DEMFEM
FEM_solver_constructor = __import__(FEMSolverSettings.solver_type)                 #DEMFEM
FEM_main_step_solver = FEM_solver_constructor.CreateSolver(rigid_face_model_part, FEMSolverSettings) #DEMFEM

# Constructing a utilities objects
creator_destructor = ParticleCreatorDestructor()
dem_fem_search = DEM_FEM_Search()

#Getting chosen scheme:
if (DEM_parameters["IntegrationScheme"].GetString() == 'Forward_Euler'):
    scheme = ForwardEulerScheme()
elif (DEM_parameters["IntegrationScheme"].GetString() == 'Symplectic_Euler'):
    scheme = SymplecticEulerScheme()
elif (DEM_parameters["IntegrationScheme"].GetString() == 'Taylor_Scheme'):
    scheme = TaylorScheme()
elif (DEM_parameters["IntegrationScheme"].GetString() == 'Newmark_Beta_Method'):
    scheme = NewmarkBetaScheme(0.5, 0.25)
elif (DEM_parameters["IntegrationScheme"].GetString() == 'Verlet_Velocity'):
    scheme = VerletVelocityScheme()
else:
    KRATOSprint('Error: selected scheme not defined. Please select a different scheme')


# Creating a solver object and set the search strategy
solver = SolverStrategy.ExplicitStrategy(spheres_model_part, rigid_face_model_part, cluster_model_part, DEM_inlet_model_part, creator_destructor, dem_fem_search, scheme, DEM_parameters, procedures)

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
FEM_solver_constructor.AddVariables(rigid_face_model_part, FEMSolverSettings)       #DEMFEM
#procedures.AddCommonVariables(rigid_face_model_part, DEM_parameters)               #DEMFEM
procedures.AddElasticFaceVariables(rigid_face_model_part, DEM_parameters)           #DEMFEM
procedures.AddMpiVariables(rigid_face_model_part)

# Reading the model_part
spheres_mp_filename   = DEM_parameters["problem_name"].GetString() + "DEM"
model_part_io_spheres = model_part_reader(spheres_mp_filename)

if "do_not_perform_initial_partition" in DEM_parameters.keys() and DEM_parameters["do_not_perform_initial_partition"].GetBool():
    pass
else:
    parallelutils.PerformInitialPartition(model_part_io_spheres)

[model_part_io_spheres, spheres_model_part, MPICommSetup] = parallelutils.SetCommunicator(spheres_model_part, model_part_io_spheres, spheres_mp_filename)


model_part_io_spheres.ReadModelPart(spheres_model_part)


max_node_Id = creator_destructor.FindMaxNodeIdInModelPart(spheres_model_part)
max_elem_Id = creator_destructor.FindMaxElementIdInModelPart(spheres_model_part)
old_max_elem_Id_spheres = max_elem_Id
max_cond_Id = creator_destructor.FindMaxElementIdInModelPart(spheres_model_part)
rigidFace_mp_filename = DEM_parameters["problem_name"].GetString() + "DEM_FEM_boundary"
model_part_io_fem = model_part_reader(rigidFace_mp_filename,max_node_Id+1, max_elem_Id+1, max_cond_Id+1)
model_part_io_fem.ReadModelPart(rigid_face_model_part)


max_node_Id = creator_destructor.FindMaxNodeIdInModelPart(rigid_face_model_part)
max_elem_Id = creator_destructor.FindMaxElementIdInModelPart(rigid_face_model_part)
max_cond_Id = creator_destructor.FindMaxElementIdInModelPart(rigid_face_model_part)
clusters_mp_filename = DEM_parameters["problem_name"].GetString() + "DEM_Clusters"
model_part_io_clusters = model_part_reader(clusters_mp_filename,max_node_Id+1, max_elem_Id+1, max_cond_Id+1)
model_part_io_clusters.ReadModelPart(cluster_model_part)
max_elem_Id = creator_destructor.FindMaxElementIdInModelPart(spheres_model_part)
if(max_elem_Id != old_max_elem_Id_spheres):
    creator_destructor.RenumberElementIdsFromGivenValue(cluster_model_part, max_elem_Id)

max_node_Id = creator_destructor.FindMaxNodeIdInModelPart(cluster_model_part)
max_elem_Id = creator_destructor.FindMaxElementIdInModelPart(cluster_model_part)
max_cond_Id = creator_destructor.FindMaxElementIdInModelPart(cluster_model_part)
DEM_Inlet_filename = DEM_parameters["problem_name"].GetString() + "DEM_Inlet"  
model_part_io_demInlet = model_part_reader(DEM_Inlet_filename,max_node_Id+1, max_elem_Id+1, max_cond_Id+1)
model_part_io_demInlet.ReadModelPart(DEM_inlet_model_part)


# Setting up the buffer size
spheres_model_part.SetBufferSize(1)
cluster_model_part.SetBufferSize(1)
DEM_inlet_model_part.SetBufferSize(1)
rigid_face_model_part.SetBufferSize(FEM_main_step_solver.buffer_size)         #DEMFFEM

# Adding dofs
solver.AddDofs(spheres_model_part)
solver.AddDofs(cluster_model_part)
solver.AddDofs(DEM_inlet_model_part)
FEM_solver_constructor.AddDofs(rigid_face_model_part, FEMSolverSettings)        

#Utilities

# set the constitutive law
constitutive_law = constitutive_law_utils.ConstitutiveLawUtility(rigid_face_model_part, DEM_parameters["Dimension"].GetInt());
constitutive_law.Initialize();
conditions    = condition_utils.ConditionsUtility(rigid_face_model_part, DEM_parameters["Dimension"].GetInt(), FEM_general_variables.Incremental_Displacement, FEM_general_variables.Incremental_Load, FEMSolverSettings.RotationDofs)

# Creating necessary directories
main_path = os.getcwd()
[post_path, data_and_results, graphs_path, MPI_results] = procedures.CreateDirectories(str(main_path), str(DEM_parameters["problem_name"].GetString()))

os.chdir(main_path)

KRATOSprint("\nInitializing Problem...")

# Initialize GiD-IO
demio.AddGlobalVariables()
demio.AddSpheresVariables()
demio.AddSpheresAndClustersVariables()
demio.AddSpheresNotInClusterAndClustersVariables()
demio.AddFEMBoundaryVariables()
demio.AddClusterVariables()
demio.AddContactVariables()
# MPI
demio.AddMpiVariables()

demio.Configure(DEM_parameters["problem_name"].GetString(),
                DEM_parameters["OutputFileType"].GetString(),
                DEM_parameters["Multifile"].GetString(),
                DEM_parameters["ContactMeshOption"].GetBool())
demio.SetOutputName(DEM_parameters["problem_name"].GetString())

os.chdir(post_path)

multifiles = (
    DEM_procedures.MultifileList(DEM_parameters["problem_name"].GetString(), 1),
    DEM_procedures.MultifileList(DEM_parameters["problem_name"].GetString(), 2),
    DEM_procedures.MultifileList(DEM_parameters["problem_name"].GetString(), 5),
    DEM_procedures.MultifileList(DEM_parameters["problem_name"].GetString(),10),
    DEM_procedures.MultifileList(DEM_parameters["problem_name"].GetString(),20),
    DEM_procedures.MultifileList(DEM_parameters["problem_name"].GetString(),50),
    )

demio.SetMultifileLists(multifiles)
os.chdir(post_path)

demio.InitializeMesh(spheres_model_part,
                     rigid_face_model_part,
                     cluster_model_part,
                     solver.contact_model_part,
                     mapping_model_part)

# Perform a partition to balance the problem
solver.search_strategy = parallelutils.GetSearchStrategy(solver, spheres_model_part)
solver.CreateCPlusPlusStrategy()
solver.RebuildListOfDiscontinuumSphericParticles()
solver.SetNormalRadiiOnAllParticles()
solver.SetSearchRadiiOnAllParticles()
parallelutils.Repart(spheres_model_part)
#parallelutils.CalculateModelNewIds(spheres_model_part)

os.chdir(post_path)

#Setting up the BoundingBox
bounding_box_time_limits = []
if DEM_parameters["BoundingBoxOption"].GetBool():
    procedures.SetBoundingBox(spheres_model_part, cluster_model_part, rigid_face_model_part, creator_destructor)
    bounding_box_time_limits = [solver.bounding_box_start_time, solver.bounding_box_stop_time]

#Creating a solver object and set the search strategy
#solver                 = SolverStrategy.ExplicitStrategy(spheres_model_part, rigid_face_model_part, cluster_model_part, DEM_inlet_model_part, creator_destructor, DEM_parameters)

dt = DEM_parameters["MaxTimeStep"].GetDouble()

#Finding the max id of the nodes... (it is necessary for anything that will add spheres to the spheres_model_part, for instance, the INLETS and the CLUSTERS read from mdpa file.
max_node_Id = creator_destructor.FindMaxNodeIdInModelPart(spheres_model_part)
max_elem_Id = creator_destructor.FindMaxElementIdInModelPart(spheres_model_part)
max_FEM_node_Id = creator_destructor.FindMaxNodeIdInModelPart(rigid_face_model_part)
max_cluster_node_Id = creator_destructor.FindMaxNodeIdInModelPart(cluster_model_part)

max_Id = max(max_FEM_node_Id, max_node_Id, max_elem_Id, max_cluster_node_Id)
creator_destructor.SetMaxNodeId(max_Id)    

#Strategy Initialization
os.chdir(main_path)
solver.Initialize() # Possible modifications of number of elements and number of nodes

dt_dem = min(DEM_parameters["MaxTimeStep"].GetDouble(), spheres_model_part.ProcessInfo.GetValue(DELTA_TIME)) # Possible modifications of DELTA_TIME
print("dem calculated" +str(dt_dem))
rigid_face_model_part.ProcessInfo[DELTA_TIME] = dt #DEMFEM
FEM_main_step_solver.Initialize() 
FEM_main_step_solver.SetRestart(False)
dt_fem = rigid_face_model_part.ProcessInfo[DELTA_TIME]   #DEMFEM
print("fem calculated" +str(dt_fem))
dt = min(dt_fem,dt_dem)
print("using this one " +str(dt))

rigid_face_model_part.ProcessInfo[DELTA_TIME] = dt #DEMFEM
spheres_model_part.ProcessInfo[DELTA_TIME] = dt  #DEMFEM
rigid_face_model_part.CloneTimeStep(0.0)
conditions.Initialize(dt)

#Constructing a model part for the DEM inlet. It contains the DEM elements to be released during the simulation  
#Initializing the DEM solver must be done before creating the DEM Inlet, because the Inlet configures itself according to some options of the DEM model part
if DEM_parameters["dem_inlet_option"].GetBool(): 
    #Constructing the inlet and initializing it (must be done AFTER the spheres_model_part Initialize)    
    DEM_inlet = DEM_Inlet(DEM_inlet_model_part)    
    DEM_inlet.InitializeDEM_Inlet(spheres_model_part, creator_destructor, solver.continuum_type)
  
DEMFEMProcedures = DEM_procedures.DEMFEMProcedures(DEM_parameters, graphs_path, spheres_model_part, rigid_face_model_part)

materialTest.Initialize(DEM_parameters, procedures, solver, graphs_path, post_path, spheres_model_part, rigid_face_model_part)

KRATOSprint("Initialization Complete" + "\n")

step           = 0
time           = 0.0
time_old_print = 0.0

report.Prepare(timer, DEM_parameters["ControlTime"].GetDouble())

first_print = True; index_5 = 1; index_10 = 1; index_50 = 1; control = 0.0

coordination_number = procedures.ModelData(spheres_model_part, solver)
KRATOSprint ("Coordination Number: " + str(coordination_number) + "\n")

materialTest.PrintChart()
materialTest.PrepareDataForGraph()

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


##############################################################################
#                                                                            #
#    MAIN LOOP                                                               #
#                                                                            #
##############################################################################
report.total_steps_expected = int(DEM_parameters["FinalTime"].GetDouble() / dt)
KRATOSprint(report.BeginReport(timer))

while time < DEM_parameters["FinalTime"].GetDouble():
    
    rigid_face_model_part.ProcessInfo[PREVIOUS_DELTA_TIME] = dt
    time  = time + dt
    step += 1

    spheres_model_part.ProcessInfo[TIME]          = time
    spheres_model_part.ProcessInfo[DELTA_TIME]    = dt
    spheres_model_part.ProcessInfo[TIME_STEPS]    = step
    
    rigid_face_model_part.CloneTimeStep(time)
    rigid_face_model_part.ProcessInfo[TIME]       = time
    rigid_face_model_part.ProcessInfo[DELTA_TIME] = dt
    rigid_face_model_part.ProcessInfo[TIME_STEPS] = step

    cluster_model_part.ProcessInfo[TIME]          = time
    cluster_model_part.ProcessInfo[DELTA_TIME]    = dt
    cluster_model_part.ProcessInfo[TIME_STEPS]    = step
    
    
    #### SOLVE #########################################
    solver.Solve()
    FEM_main_step_solver.Solve()                            #DEMFEM
    conditions.SetIncrementalLoad(step, dt);
    ####################################################
    
    # Walls movement:
    mesh_motion.MoveAllMeshes(rigid_face_model_part, time, dt)
    mesh_motion.MoveAllMeshes(spheres_model_part, time, dt)
    mesh_motion.MoveAllMeshes(DEM_inlet_model_part, time, dt)

    
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
    DEMFEMProcedures.PrintBallsGraph(time)

    #### GiD IO ##########################################
    time_to_print = time - time_old_print

    if (DEM_parameters["OutputTimeStep"].GetDouble() - time_to_print < 1e-2 * dt):
        
        if solver.poisson_ratio_option:
            DEMFEMProcedures.PrintPoisson(spheres_model_part, DEM_parameters, "Poisson_ratio.txt", time)
            
        if DEM_parameters["PostEulerAngles"].GetBool():
            post_utils.PrintEulerAngles(spheres_model_part)

        KRATOSprint("*******************  PRINTING RESULTS FOR GID  ***************************")
        KRATOSprint("                        ("+ str(spheres_model_part.NumberOfElements(0)) + " elements)")
        KRATOSprint("                        ("+ str(spheres_model_part.NumberOfNodes(0)) + " nodes)")
        KRATOSprint("")

        os.chdir(data_and_results)

        #properties_list = ProceduresMonitorPhysicalProperties(spheres_model_part, physics_calculator, properties_list)

        demio.PrintMultifileLists(time, post_path)
        #os.chdir(main_path)

        os.chdir(post_path)

        solver.PrepareElementsForPrinting()
        if DEM_parameters["ContactMeshOption"].GetBool():
            solver.PrepareContactElementsForPrinting()
        
        demio.PrintResults(all_model_parts, creator_destructor, dem_fem_search, time, bounding_box_time_limits)
        os.chdir(main_path)
        
        if nodeplotter:
            plotter.plot_variables(time) #Related to debugging

        time_old_print = time
    
    #AFTER PRINTING OPERATIONS
    conditions.RestartImposedDisp()

    #if((step%500) == 0):
        #if (( DEM_parameters["ContactMeshOption"].GetBool()) and (DEM_parameters["TestType"].GetString() != "None"))  :
            #MaterialTest.OrientationStudy(solver.contact_model_part, step)
    
    #print("TIME STEP ENDS +++++++++++++++++++++++++++++++++++++++++++++++++")
##############################################################################
#                                                                            #
#    FINALIZATION                                                            #
#                                                                            #
##############################################################################

# Print times and more info
KRATOSprint(report.FinalReport(timer))

demio.FinalizeMesh()
materialTest.FinalizeGraphs()
DEMFEMProcedures.FinalizeGraphs(rigid_face_model_part)
DEMFEMProcedures.FinalizeBallsGraphs(spheres_model_part)

demio.CloseMultifiles()

os.chdir(main_path)

# Freeing up memory
objects_to_destroy = [demio, procedures, creator_destructor, dem_fem_search, solver, DEMFEMProcedures, post_utils, 
                      cluster_model_part, rigid_face_model_part, spheres_model_part, DEM_inlet_model_part, mapping_model_part]

if DEM_parameters["dem_inlet_option"].GetBool(): 
    objects_to_destroy.append(DEM_inlet)

for obj in objects_to_destroy:
    del obj

procedures.DeleteFiles()


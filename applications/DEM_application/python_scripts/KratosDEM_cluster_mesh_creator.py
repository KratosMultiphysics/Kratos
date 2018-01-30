from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Python imports
import time as timer
import os
import sys
import cluster_mdpa_creator

# Kratos
from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *

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
        #return ModelPartIO(modelpart)
        return ReorderConsecutiveFromGivenIdsModelPartIO(modelpart, nodeid, elemid, condid)
    

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
elif (DEM_parameters["ElementType"].GetString() == "IceContPartDEMElement3D"):
    import ice_continuum_sphere_strategy as SolverStrategy
else:
    KRATOSprint('Error: Strategy unavailable. Select a different scheme-element')
        
##############################################################################
#                                                                            #
#    INITIALIZE                                                              #
#                                                                            #
##############################################################################

procedures    = DEM_procedures.Procedures(DEM_parameters)

# Creating necessary directories
main_path = os.getcwd()
[post_path, data_and_results, graphs_path, MPI_results] = procedures.CreateDirectories(str(main_path), str(DEM_parameters["problem_name"].GetString()))

demio         = DEM_procedures.DEMIo(DEM_parameters, post_path)
report        = DEM_procedures.Report()
parallelutils = DEM_procedures.ParallelUtils()
materialTest  = DEM_procedures.MaterialTest()
 
# Set the print function TO_DO: do this better...
KRATOSprint   = procedures.KRATOSprint

# Prepare modelparts
spheres_model_part    = ModelPart("SpheresPart")
rigid_face_model_part = ModelPart("RigidFacePart")
cluster_model_part    = ModelPart("ClusterPart")
DEM_inlet_model_part  = ModelPart("DEMInletPart")
mapping_model_part    = ModelPart("MappingPart")
contact_model_part    = ModelPart("ContactPart")
mp_list = []
mp_list.append(spheres_model_part)
mp_list.append(rigid_face_model_part)
mp_list.append(cluster_model_part)
mp_list.append(DEM_inlet_model_part)
mp_list.append(mapping_model_part)
mp_list.append(contact_model_part)

all_model_parts = DEM_procedures.SetOfModelParts(mp_list)

# Constructing a utilities objects
creator_destructor = ParticleCreatorDestructor()
dem_fem_search = DEM_FEM_Search()

scheme = procedures.SetScheme()
#solver = SolverStrategy.ExplicitStrategy(spheres_model_part, rigid_face_model_part, cluster_model_part, DEM_inlet_model_part, contact_model_part, creator_destructor, dem_fem_search, scheme, DEM_parameters, procedures)
solver = SolverStrategy.ExplicitStrategy(all_model_parts, creator_destructor, dem_fem_search, scheme, DEM_parameters, procedures)

procedures.AddAllVariablesInAllModelParts(solver, scheme, all_model_parts, DEM_parameters)

os.chdir(main_path)
# Reading the model_part
spheres_mp_filename   = DEM_parameters["problem_name"].GetString() + "DEM"
model_part_io_spheres = model_part_reader(spheres_mp_filename)

if "do_not_perform_initial_partition" in DEM_parameters.keys() and DEM_parameters["do_not_perform_initial_partition"].GetBool():
    pass
else:
    parallelutils.PerformInitialPartition(model_part_io_spheres)

os.chdir(main_path)
[model_part_io_spheres, spheres_model_part, MPICommSetup] = parallelutils.SetCommunicator(spheres_model_part, model_part_io_spheres, spheres_mp_filename)

model_part_io_spheres.ReadModelPart(spheres_model_part)

max_node_Id = creator_destructor.FindMaxNodeIdInModelPart(spheres_model_part)
max_elem_Id = creator_destructor.FindMaxElementIdInModelPart(spheres_model_part)
old_max_elem_Id_spheres = max_elem_Id
max_cond_Id = creator_destructor.FindMaxConditionIdInModelPart(spheres_model_part)
rigidFace_mp_filename = DEM_parameters["problem_name"].GetString() + "DEM_FEM_boundary"
model_part_io_fem = model_part_reader(rigidFace_mp_filename,max_node_Id+1, max_elem_Id+1, max_cond_Id+1)
model_part_io_fem.ReadModelPart(rigid_face_model_part)

max_node_Id = creator_destructor.FindMaxNodeIdInModelPart(rigid_face_model_part)
max_elem_Id = creator_destructor.FindMaxElementIdInModelPart(rigid_face_model_part)
max_cond_Id = creator_destructor.FindMaxConditionIdInModelPart(rigid_face_model_part)
clusters_mp_filename = DEM_parameters["problem_name"].GetString() + "DEM_Clusters"
cluster_mdpa_creator.WriteClusterMdpa()
model_part_io_clusters = model_part_reader(clusters_mp_filename,max_node_Id+1, max_elem_Id+1, max_cond_Id+1)
model_part_io_clusters.ReadModelPart(cluster_model_part)
max_elem_Id = creator_destructor.FindMaxElementIdInModelPart(spheres_model_part)
if (max_elem_Id != old_max_elem_Id_spheres):
    creator_destructor.RenumberElementIdsFromGivenValue(cluster_model_part, max_elem_Id)

max_node_Id = creator_destructor.FindMaxNodeIdInModelPart(cluster_model_part)
max_elem_Id = creator_destructor.FindMaxElementIdInModelPart(cluster_model_part)
max_cond_Id = creator_destructor.FindMaxConditionIdInModelPart(cluster_model_part)
DEM_Inlet_filename = DEM_parameters["problem_name"].GetString() + "DEM_Inlet"  
model_part_io_demInlet = model_part_reader(DEM_Inlet_filename,max_node_Id+1, max_elem_Id+1, max_cond_Id+1)
model_part_io_demInlet.ReadModelPart(DEM_inlet_model_part)

# Setting up the buffer size
procedures.SetUpBufferSizeInAllModelParts(spheres_model_part, 1, cluster_model_part, 1, DEM_inlet_model_part, 1, rigid_face_model_part, 1)
# Adding dofs
solver.AddDofs(spheres_model_part)
solver.AddDofs(cluster_model_part)
solver.AddDofs(DEM_inlet_model_part)

os.chdir(main_path)

KRATOSprint("\nInitializing Problem...")

demio.Initialize(DEM_parameters)

os.chdir(post_path)
demio.InitializeMesh(all_model_parts)

# Perform a partition to balance the problem
solver.search_strategy = parallelutils.GetSearchStrategy(solver, spheres_model_part)
solver.BeforeInitialize()
parallelutils.Repart(spheres_model_part)

#Setting up the BoundingBox
bounding_box_time_limits = procedures.SetBoundingBoxLimits(all_model_parts, creator_destructor)

dt = DEM_parameters["MaxTimeStep"].GetDouble()

#Finding the max id of the nodes... (it is necessary for anything that will add spheres to the spheres_model_part, for instance, the INLETS and the CLUSTERS read from mdpa file.z
max_Id = procedures.FindMaxNodeIdAccrossModelParts(creator_destructor, all_model_parts)
creator_destructor.SetMaxNodeId(max_Id)    

#Strategy Initialization
os.chdir(main_path)
solver.Initialize() # Possible modifications of number of elements and number of nodes
#dt = min(DEM_parameters["MaxTimeStep"].GetDouble(), spheres_model_part.ProcessInfo.GetValue(DELTA_TIME)) # under revision. linked to automatic timestep? Possible modifications of DELTA_TIME
dt = DEM_parameters["MaxTimeStep"].GetDouble()
#Constructing a model part for the DEM inlet. It contains the DEM elements to be released during the simulation  
#Initializing the DEM solver must be done before creating the DEM Inlet, because the Inlet configures itself according to some options of the DEM model part
if DEM_parameters["dem_inlet_option"].GetBool(): 
    #Constructing the inlet and initializing it (must be done AFTER the spheres_model_part Initialize)    
    DEM_inlet = DEM_Inlet(DEM_inlet_model_part)    
    DEM_inlet.InitializeDEM_Inlet(spheres_model_part, creator_destructor, solver.continuum_type)
  
DEMFEMProcedures = DEM_procedures.DEMFEMProcedures(DEM_parameters, graphs_path, spheres_model_part, rigid_face_model_part)

os.chdir(graphs_path)
DEMEnergyCalculator = DEM_procedures.DEMEnergyCalculator(DEM_parameters, spheres_model_part, cluster_model_part, "EnergyPlot.grf")

materialTest.Initialize(DEM_parameters, procedures, solver, graphs_path, post_path, spheres_model_part, rigid_face_model_part)

KRATOSprint("Initialization Complete" + "\n")

step           = 0
time           = 0.0
time_old_print = 0.0

report.Prepare(timer, DEM_parameters["ControlTime"].GetDouble())

procedures.ModelData(spheres_model_part, solver)

materialTest.PrintChart()
materialTest.PrepareDataForGraph()

post_utils = DEM_procedures.PostUtils(DEM_parameters, spheres_model_part)

step = 0

##############################################################################
#    MAIN LOOP                                                               #
##############################################################################
report.total_steps_expected = int(DEM_parameters["FinalTime"].GetDouble() / dt)
KRATOSprint(report.BeginReport(timer))

while time < DEM_parameters["FinalTime"].GetDouble():
    
    dt    = spheres_model_part.ProcessInfo.GetValue(DELTA_TIME) # Possible modifications of DELTA_TIME
    time  = time + dt
    step += 1

    DEMFEMProcedures.UpdateTimeInModelParts(all_model_parts, time,dt,step) 
    
    #### SOLVE #########################################
    solver.Solve()
    ####################################################
    
    DEMFEMProcedures.MoveAllMeshes(all_model_parts, time, dt)
       
    ##### adding DEM elements by the inlet ######
    if DEM_parameters["dem_inlet_option"].GetBool(): 
        DEM_inlet.CreateElementsFromInletMesh(spheres_model_part, cluster_model_part, creator_destructor)  # After solving, to make sure that neighbours are already set.              

    stepinfo = report.StepiReport(timer,time,step)
    if stepinfo:
        KRATOSprint(stepinfo)
    
    #### PRINTING GRAPHS ####
    os.chdir(graphs_path)
    post_utils.ComputeMeanVelocitiesinTrap("Average_Velocity.txt", time)

    materialTest.MeasureForcesAndPressure()
    materialTest.PrintGraph(time)
    
    DEMFEMProcedures.PrintGraph(time)
    DEMFEMProcedures.PrintBallsGraph(time)

    DEMEnergyCalculator.CalculateEnergyAndPlot(time)

    #### GiD IO ##########################################
    time_to_print = time - time_old_print

    if (DEM_parameters["OutputTimeStep"].GetDouble() - time_to_print < 1e-2 * dt):
        
        if solver.poisson_ratio_option:
            DEMFEMProcedures.PrintPoisson(spheres_model_part, DEM_parameters, "Poisson_ratio.txt", time)
            
        if DEM_parameters["PostEulerAngles"].GetBool():
            post_utils.PrintEulerAngles(spheres_model_part, cluster_model_part)

        demio.ShowPrintingResultsOnScreen(all_model_parts)

        os.chdir(data_and_results)
        demio.PrintMultifileLists(time, post_path)
        os.chdir(post_path)

        solver.PrepareElementsForPrinting()
        if DEM_parameters["ContactMeshOption"].GetBool():
            solver.PrepareContactElementsForPrinting()
        
        demio.PrintResults(all_model_parts, creator_destructor, dem_fem_search, time, bounding_box_time_limits)
        os.chdir(main_path)

        time_old_print = time
   
KRATOSprint("Finalizing execution...")

demio.FinalizeMesh()
materialTest.FinalizeGraphs()
DEMFEMProcedures.FinalizeGraphs(rigid_face_model_part)
DEMFEMProcedures.FinalizeBallsGraphs(spheres_model_part)

DEMEnergyCalculator.FinalizeEnergyPlot()

demio.CloseMultifiles()

os.chdir(main_path)

objects_to_destroy = [demio, procedures, creator_destructor, dem_fem_search, solver, DEMFEMProcedures, post_utils, 
                      cluster_model_part, rigid_face_model_part, spheres_model_part, DEM_inlet_model_part, mapping_model_part]

if DEM_parameters["dem_inlet_option"].GetBool(): 
    objects_to_destroy.append(DEM_inlet)

for obj in objects_to_destroy:
    del obj

procedures.DeleteFiles()

KRATOSprint(report.FinalReport(timer))

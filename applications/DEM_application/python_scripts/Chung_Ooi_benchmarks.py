from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import time as timer
import os
import sys
from math import sin, cos
import shutil

# Kratos
from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *

import plot_variables           # Related to benchmarks in Chung, Ooi
import Chung_Ooi_class as COC   # Related to benchmarks in Chung, Ooi
from glob import glob

sys.path.insert(0,'')
# DEM Application

import DEM_explicit_solver_var as DEM_parameters

# Strategy object
if   (DEM_parameters["ElementType"].GetString() == "SphericPartDEMElement3D"     or DEM_parameters["ElementType"].GetString() == "CylinderPartDEMElement2D"):
    import sphere_strategy as SolverStrategy
elif (DEM_parameters["ElementType"].GetString() == "SphericContPartDEMElement3D" or DEM_parameters["ElementType"].GetString() == "CylinderContPartDEMElement2D"):
    import continuum_sphere_strategy as SolverStrategy
elif (DEM_parameters["ElementType"].GetString() == "ThermalSphericContPartDEMElement3D"):
    import thermal_continuum_sphere_strategy as SolverStrategy

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

########################################################## CHUNG, OOI BENCHMARKS
### BENCHMARK NUMBER

benchmark_number = int(sys.argv[1])

if benchmark_number==1 or benchmark_number==2:
    nodeplotter = 1
else:
    nodeplotter = 0
    
####################
final_time, dt, output_time_step, number_of_points_in_the_graphic = COC.initialize_time_parameters(benchmark_number)
benchmark_class_name = 'Benchmark' + str(benchmark_number)
benchmark_class = getattr(COC, benchmark_class_name)
benchmark = benchmark_class()
#
for iteration in range(1, number_of_points_in_the_graphic + 1):
#
########################################################## CHUNG, OOI BENCHMARKS

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

    # Constructing a creator/destructor object
    creator_destructor = ParticleCreatorDestructor()
    # Creating a solver object and set the search strategy
    solver             = SolverStrategy.ExplicitStrategy(spheres_model_part, rigid_face_model_part, cluster_model_part, DEM_inlet_model_part, creator_destructor, DEM_parameters)
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
    #
    procedures.AddCommonVariables(rigid_face_model_part, DEM_parameters)
    procedures.AddRigidFaceVariables(rigid_face_model_part, DEM_parameters)
    procedures.AddMpiVariables(rigid_face_model_part)
    
    problem_name                = 'benchmark' + str(benchmark_number)
    spheres_mp_filename         = problem_name + "DEM"
    DEM_parameters["problem_name"] = problem_name
    
    model_part_io_spheres = ModelPartIO(spheres_mp_filename)
    
    if "do_not_perform_initial_partition" in DEM_parameters.keys() and DEM_parameters["do_not_perform_initial_partition"].GetBool(): #TODO: add subIndex _option to this variable? Or the type is enough? Discuss with the team
        pass
    else:
        parallelutils.PerformInitialPartition(model_part_io_spheres)

    [model_part_io_spheres, spheres_model_part, MPICommSetup] = parallelutils.SetCommunicator(spheres_model_part, model_part_io_spheres, spheres_mp_filename)
        
    model_part_io_spheres.ReadModelPart(spheres_model_part)
        
    ###################################################### CHUNG, OOI BENCHMARKS
    
    benchmark.get_initial_data(spheres_model_part, iteration, number_of_points_in_the_graphic)
        
    ###################################################### CHUNG, OOI BENCHMARKS
            
    rigidFace_mp_filename = DEM_parameters["problem_name"].GetString() + "DEM_FEM_boundary"
    model_part_io_fem = ModelPartIO(rigidFace_mp_filename)
    model_part_io_fem.ReadModelPart(rigid_face_model_part)

    clusters_mp_filename = 'benchmark' + "DEM_Clusters"
    model_part_io_clusters = ModelPartIO(clusters_mp_filename)
    model_part_io_clusters.ReadModelPart(cluster_model_part)

    DEM_Inlet_filename = 'benchmark' + "DEM_Inlet"  
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
    [post_path, data_and_results, graphs_path, MPI_results] = procedures.CreateDirectories(str(main_path), str(DEM_parameters["problem_name"].GetString()))

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
        DEM_procedures.MultifileList(DEM_parameters["problem_name"].GetString(),2 ),
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

    solver.Initialize()    # Possible modifications of DELTA_TIME
        
    if (DEM_parameters["ContactMeshOption"].GetBool()):
        contact_model_part = solver.contact_model_part

    # constructing a model part for the DEM inlet. it contains the DEM elements to be released during the simulation  
    # Initializing the DEM solver must be done before creating the DEM Inlet, because the Inlet configures itself according to some options of the DEM model part

    if DEM_parameters["dem_inlet_option"].GetBool():    
        max_node_Id = creator_destructor.FindMaxNodeIdInModelPart(spheres_model_part)
        max_FEM_node_Id = creator_destructor.FindMaxNodeIdInModelPart(rigid_face_model_part)

        if ( max_FEM_node_Id > max_node_Id):
            max_node_Id = max_FEM_node_Id

        creator_destructor.SetMaxNodeId(max_node_Id)                            

        for properties in DEM_inlet_model_part.Properties:

                DiscontinuumConstitutiveLawString = properties[DEM_DISCONTINUUM_CONSTITUTIVE_LAW_NAME];
                DiscontinuumConstitutiveLaw = globals().get(DiscontinuumConstitutiveLawString)()
                DiscontinuumConstitutiveLaw.SetConstitutiveLawInProperties(properties, True)             

        # constructing the inlet and intializing it (must be done AFTER the spheres_model_part Initialize)    
        DEM_inlet = DEM_Inlet(DEM_inlet_model_part)    
        DEM_inlet.InitializeDEM_Inlet(spheres_model_part, creator_destructor)

    DEMFEMProcedures = DEM_procedures.DEMFEMProcedures(DEM_parameters, graphs_path, spheres_model_part, rigid_face_model_part)

    materialTest.Initialize(DEM_parameters, procedures, solver, graphs_path, post_path, spheres_model_part, rigid_face_model_part)

    KRATOSprint("Initialization Complete" + "\n")

    step           = 0
    time           = 0.0
    time_old_print = 0.0

    report.Prepare(timer,DEM_parameters["ControlTime"].GetDouble())

    first_print  = True; index_5 = 1; index_10  = 1; index_50  = 1; control = 0.0

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

    report.total_steps_expected = int(final_time / dt)

    KRATOSprint(report.BeginReport(timer))

    print("Computing points in the curve...", 1 + number_of_points_in_the_graphic - iteration, "point(s) left to finish....",'\n')

    mesh_motion = DEMFEMUtilities()

    # creating a Post Utils object that executes several post-related tasks
    post_utils = DEM_procedures.PostUtils(DEM_parameters, spheres_model_part)
    
    list_of_nodes_ids = [1]
       
    if nodeplotter:
        os.chdir(main_path)
        plotter = plot_variables.variable_plotter(spheres_model_part, list_of_nodes_ids) #Related to the benchmark in Chung, Ooi
    
    step = 0  

    while (time < final_time):
        #
        #SLS dt   = spheres_model_part.ProcessInfo.GetValue(DELTA_TIME) # Possible modifications of DELTA_TIME
        time = time + dt
        step += 1

        spheres_model_part.ProcessInfo[TIME]            = time
        spheres_model_part.ProcessInfo[DELTA_TIME]      = dt
        spheres_model_part.ProcessInfo[TIME_STEPS]      = step

        rigid_face_model_part.ProcessInfo[TIME]         = time
        rigid_face_model_part.ProcessInfo[DELTA_TIME]   = dt
        rigid_face_model_part.ProcessInfo[TIME_STEPS]   = step

        cluster_model_part.ProcessInfo[TIME]            = time
        cluster_model_part.ProcessInfo[DELTA_TIME]      = dt
        cluster_model_part.ProcessInfo[TIME_STEPS]      = step

        #walls movement:
        mesh_motion.MoveAllMeshes(rigid_face_model_part, time, dt)
        #########mesh_motion.MoveSphereMeshes(spheres_model_part, time, dt)

        #### SOLVE #########################################
        
        solver.Solve()
                
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
            post_utils.ComputeMeanVelocitiesinTrap("Average_Velocity", time)

        #### MATERIAL TEST GRAPHS ############################
        materialTest.MeasureForcesAndPressure()
        materialTest.PrintGraph(time)

        #### GENERAL FORCE GRAPHS ############################
        DEMFEMProcedures.PrintGraph(time)
        DEMFEMProcedures.PrintBallsGraph(time)  

        #### GiD IO ##########################################
        time_to_print = time - time_old_print

        if (output_time_step - time_to_print < 1e-2*dt):            

            '''KRATOSprint("*******************  PRINTING RESULTS FOR GID  ***************************")
            KRATOSprint("                        ("+ str(spheres_model_part.NumberOfElements(0)) + " elements)")
            KRATOSprint("                        ("+ str(spheres_model_part.NumberOfNodes(0)) + " nodes)")
            KRATOSprint("")'''

            os.chdir(data_and_results)

            demio.PrintMultifileLists(time, post_path)
            os.chdir(main_path)

            os.chdir(post_path)

            if DEM_parameters["ContactMeshOption"].GetBool():
                solver.PrepareContactElementsForPrinting()

            demio.PrintResults(mixed_model_part, spheres_model_part, rigid_face_model_part, cluster_model_part, contact_model_part, mapping_model_part, time)

            os.chdir(main_path)

            time_old_print = time
                
        if nodeplotter:
            os.chdir(main_path)
            plotter.plot_variables(time) #Related to the benchmark in Chung, Ooi

    ############################################################################
    #                                                                          #
    #    FINALIZATION                                                          #
    #                                                                          #
    ############################################################################
    #
    #
    ###################################################### CHUNG, OOI BENCHMARKS
                               
    benchmark.get_final_data(spheres_model_part)
        
    ###################################################### CHUNG, OOI BENCHMARKS
    #
        
    demio.FinalizeMesh()
    materialTest.FinalizeGraphs()
    DEMFEMProcedures.FinalizeGraphs(rigid_face_model_part)
    DEMFEMProcedures.FinalizeBallsGraphs(spheres_model_part)

    demio.CloseMultifiles()

    if nodeplotter:
        os.chdir(main_path)
        plotter.close_files() #Related to the benchmark in Chung, Ooi

    os.chdir(main_path)
                
    # Print times and more info
      
    KRATOSprint(report.FinalReport(timer))
    
benchmark.print_results(number_of_points_in_the_graphic, dt)

COC.delete_archives(nodeplotter) #.......Removing some unuseful files 


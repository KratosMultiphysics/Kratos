from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import time as timer
import os
import sys
#from math import sin, cos
#import shutil

# Kratos
from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *

import plot_variables                # Related to benchmarks in Chung, Ooi
import DEM_benchmarks_class as DBC   # Related to benchmarks in Chung, Ooi
#from glob import glob

sys.path.insert(0,'')
# DEM Application

benchmark_number = int(sys.argv[1])

list = list(range(1,12))

if benchmark_number in list:
    import DEM_explicit_solver_var as DEM_parameters
elif benchmark_number ==12:
    import DEM_explicit_solver_var_12 as DEM_parameters
else:
    import DEM_explicit_solver_var2 as DEM_parameters


# Strategy object
if   (DEM_parameters.ElementType == "SphericPartDEMElement3D"     or DEM_parameters.ElementType == "CylinderPartDEMElement2D"):
    import sphere_strategy as SolverStrategy
elif (DEM_parameters.ElementType == "SphericContPartDEMElement3D" or DEM_parameters.ElementType == "CylinderContPartDEMElement2D"):
    import continuum_sphere_strategy as SolverStrategy
elif (DEM_parameters.ElementType == "ThermalSphericContPartDEMElement3D"):
    import thermal_continuum_sphere_strategy as SolverStrategy

# Import MPI modules if needed. This way to do this is only valid when using OpenMPI. For other implementations of MPI it will not work.
if "OMPI_COMM_WORLD_SIZE" in os.environ:
    # Kratos MPI
    from KratosMultiphysics.MetisApplication import *
    from KratosMultiphysics.MPISearchApplication import *
    from KratosMultiphysics.mpi import *

    # DEM Application MPI
    import DEM_procedures_mpi as DEM_procedures
    #import DEM_material_test_script_mpi as DEM_material_test_script
    print("Running under MPI...........")
else:
    # DEM Application
    import DEM_procedures
    #import DEM_material_test_script

    print("Running under OpenMP........")

########################################################## CHUNG, OOI BENCHMARKS
    
nodeplotter = 0  
  
if benchmark_number in list:
    nodeplotter = 1  
####################
final_time, dt, output_time_step, number_of_points_in_the_graphic, number_of_coeffs_of_restitution = DBC.initialize_time_parameters(benchmark_number)
benchmark_class_name = 'Benchmark' + str(benchmark_number)
benchmark_class = getattr(DBC, benchmark_class_name)
benchmark = benchmark_class()
#
for coeff_of_restitution_iteration in range(1, number_of_coeffs_of_restitution + 1):
    
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

        # Constructing utilities objects
        creator_destructor = ParticleCreatorDestructor()
        dem_fem_search = DEM_FEM_Search()
        
        # Creating a solver object and set the search strategy
        solver = SolverStrategy.ExplicitStrategy(spheres_model_part, rigid_face_model_part, cluster_model_part, DEM_inlet_model_part, creator_destructor, dem_fem_search, DEM_parameters)
        
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
        DEM_parameters.problem_name = problem_name

        model_part_io_spheres = ModelPartIO(spheres_mp_filename)

        # Perform the initial partition
        [model_part_io_spheres, spheres_model_part, MPICommSetup] = parallelutils.PerformInitialPartition(spheres_model_part, model_part_io_spheres, spheres_mp_filename)

        model_part_io_spheres.ReadModelPart(spheres_model_part)

        ###################################################### CHUNG, OOI BENCHMARKS

        benchmark.set_initial_data(spheres_model_part, iteration, number_of_points_in_the_graphic, coeff_of_restitution_iteration)

        ###################################################### CHUNG, OOI BENCHMARKS

        rigidFace_mp_filename = DEM_parameters.problem_name + "DEM_FEM_boundary"
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
        [post_path, data_and_results, graphs_path, MPI_results] = procedures.CreateDirectories(str(main_path), str(DEM_parameters.problem_name))

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
        demio.EnableMpiVariables()

        demio.Configure(DEM_parameters.problem_name,
                        DEM_parameters.OutputFileType,
                        DEM_parameters.Multifile,
                        DEM_parameters.ContactMeshOption)

        demio.SetOutputName(DEM_parameters.problem_name)

        os.chdir(post_path)

        multifiles = (
            DEM_procedures.MultifileList(DEM_parameters.problem_name,1 ),
            DEM_procedures.MultifileList(DEM_parameters.problem_name,2 ),
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

        # Setting up the BoundingBox
        bounding_box_time_limits = []
        if (DEM_parameters.BoundingBoxOption == "ON"):
            procedures.SetBoundingBox(spheres_model_part, cluster_model_part, rigid_face_model_part, creator_destructor)
            bounding_box_time_limits = [solver.bounding_box_start_time, solver.bounding_box_stop_time]
        
        # Creating a solver object and set the search strategy
        # solver                 = SolverStrategy.ExplicitStrategy(spheres_model_part, rigid_face_model_part, cluster_model_part, DEM_inlet_model_part, creator_destructor, DEM_parameters)
        solver.search_strategy = parallelutils.GetSearchStrategy(solver, spheres_model_part)

        solver.Initialize()    # Possible modifications of DELTA_TIME

        if (DEM_parameters.ContactMeshOption =="ON"):
            contact_model_part = solver.contact_model_part

        # constructing a model part for the DEM inlet. it contains the DEM elements to be released during the simulation  
        # Initializing the DEM solver must be done before creating the DEM Inlet, because the Inlet configures itself according to some options of the DEM model part

        if (DEM_parameters.dem_inlet_option):    
            max_node_Id = creator_destructor.FindMaxNodeIdInModelPart(spheres_model_part)
            max_FEM_node_Id = creator_destructor.FindMaxNodeIdInModelPart(rigid_face_model_part)

            if ( max_FEM_node_Id > max_node_Id):
                max_node_Id = max_FEM_node_Id

            creator_destructor.SetMaxNodeId(max_node_Id)                            

            for properties in DEM_inlet_model_part.Properties:

                    DiscontinuumConstitutiveLawString = properties[DEM_DISCONTINUUM_CONSTITUTIVE_LAW_NAME];
                    DiscontinuumConstitutiveLaw = globals().get(DiscontinuumConstitutiveLawString)()
                    DiscontinuumConstitutiveLaw.SetConstitutiveLawInProperties(properties)             

            # constructing the inlet and intializing it (must be done AFTER the spheres_model_part Initialize)    
            DEM_inlet = DEM_Inlet(DEM_inlet_model_part)    
            DEM_inlet.InitializeDEM_Inlet(spheres_model_part, creator_destructor)

        DEMFEMProcedures = DEM_procedures.DEMFEMProcedures(DEM_parameters, graphs_path, spheres_model_part, rigid_face_model_part)

        materialTest.Initialize(DEM_parameters, procedures, solver, graphs_path, post_path, spheres_model_part, rigid_face_model_part)

        KRATOSprint("Initialization Complete" + "\n")

        step           = 0
        time           = 0.0
        time_old_print = 0.0

        report.Prepare(timer,DEM_parameters.ControlTime)

        first_print  = True; index_5 = 1; index_10  = 1; index_50  = 1; control = 0.0

        if (DEM_parameters.ModelDataInfo == "ON"):
            os.chdir(data_and_results)
            if (DEM_parameters.ContactMeshOption == "ON"):
              (coordination_number) = procedures.ModelData(spheres_model_part, contact_model_part, solver)       # calculates the mean number of neighbours the mean radius, etc..
              KRATOSprint ("Coordination Number: " + str(coordination_number) + "\n")
              os.chdir(main_path)
            else:
              KRATOSprint("Activate Contact Mesh for ModelData information")

        if (DEM_parameters.Dempack):
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
            plotter = plot_variables.variable_plotter(spheres_model_part, list_of_nodes_ids) #Related to Benchmarks
            tang_plotter = plot_variables.tangential_force_plotter(spheres_model_part, list_of_nodes_ids, iteration) #Related to Benchmarks

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


            #### GENERAL TEST LOAD&UNLOAD  ######################
            benchmark.ApplyNodalRotation(time, spheres_model_part)

            #walls movement:
            mesh_motion.MoveAllMeshes(rigid_face_model_part, time, dt)
            #########mesh_motion.MoveSphereMeshes(spheres_model_part, time, dt)

            #### SOLVE #########################################

            solver.Solve()

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
                post_utils.ComputeMeanVelocitiesinTrap("Average_Velocity", time)

            #### MATERIAL TEST GRAPHS ############################
            materialTest.MeasureForcesAndPressure()
            materialTest.PrintGraph(time)

            #### GENERAL FORCE GRAPHS ############################
            #DEMFEMProcedures.MeasureForces()
            DEMFEMProcedures.PrintGraph(time)
            benchmark.generate_graph_points(spheres_model_part, time, output_time_step, dt)

            #### GiD IO ##########################################
            time_to_print = time - time_old_print

            if (output_time_step - time_to_print < 1e-2*dt):            

                '''KRATOSprint("*******************  PRINTING RESULTS FOR GID  ***************************")
                KRATOSprint("                        ("+ str(spheres_model_part.NumberOfElements(0)) + " elements)")
                KRATOSprint("                        ("+ str(spheres_model_part.NumberOfNodes(0)) + " nodes)")
                KRATOSprint("")
                sys.stdout.flush()'''

                os.chdir(data_and_results)

                demio.PrintMultifileLists(time, post_path)
                os.chdir(main_path)

                os.chdir(post_path)

                if (DEM_parameters.ContactMeshOption == "ON"):
                    solver.PrepareContactElementsForPrinting()

                demio.PrintResults(mixed_model_part, spheres_model_part, rigid_face_model_part, cluster_model_part, contact_model_part, mapping_model_part, creator_destructor, dem_fem_search, time, bounding_box_time_limits)

                os.chdir(main_path)

                time_old_print = time

            if nodeplotter:
                os.chdir(main_path)
                plotter.plot_variables(time) #Related to the benchmark in Chung, Ooi
                tang_plotter.plot_tangential_force(time)

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
            tang_plotter.close_files()

        os.chdir(main_path)

        # Print times and more info

        KRATOSprint(report.FinalReport(timer))
    
    benchmark.print_results(number_of_points_in_the_graphic, dt)

DBC.delete_archives() #.......Removing some unuseful files 


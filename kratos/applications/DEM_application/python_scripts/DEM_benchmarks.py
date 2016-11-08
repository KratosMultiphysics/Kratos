from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import time as timer
import os
import sys

# Kratos
from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *

import plot_variables                # Related to benchmarks in Chung, Ooi
import DEM_benchmarks_class as DBC   # Related to benchmarks in Chung, Ooi
#from glob import glob

sys.path.insert(0,'')
# DEM Application

benchmark_number = int(sys.argv[1])

listDISCONT = list(range(1,12))
listROLLFR  = list(range(12,13))
listDEMFEM  = list(range(13,18))
listCONT    = list(range(20,27))

'''
listDISCONT = list(range(3,4)) #list(range(1,12))
listROLLFR  = [] #list(range(12,13))
listDEMFEM  = [] #list(range(13,18))
listCONT    = [] #list(range(20,27))
'''

if benchmark_number in listDISCONT:
    import DEM_explicit_solver_var as DEM_parameters
elif benchmark_number in listROLLFR:
    import DEM_explicit_solver_var_ROLLFR as DEM_parameters
elif benchmark_number in listDEMFEM:
    import DEM_explicit_solver_var_DEMFEM as DEM_parameters
elif benchmark_number in listCONT:
    import DEM_explicit_solver_var_CONT as DEM_parameters
elif benchmark_number == 27:
    import DEM_explicit_solver_var_UCS as DEM_parameters
elif benchmark_number == 28:
    import DEM_explicit_solver_var_PENDULO3D as DEM_parameters
else:
    print('Benchmark number does not exist')    
    sys.exit()

# Strategy object
if   (DEM_parameters.ElementType == "SphericPartDEMElement3D"     or DEM_parameters.ElementType == "CylinderPartDEMElement2D"):
    import sphere_strategy as SolverStrategy
elif (DEM_parameters.ElementType == "SphericContPartDEMElement3D" or DEM_parameters.ElementType == "CylinderContPartDEMElement2D"):
    import continuum_sphere_strategy as SolverStrategy
elif (DEM_parameters.ElementType == "ThermalSphericContPartDEMElement3D"):
    import thermal_continuum_sphere_strategy as SolverStrategy
else:
    KRATOSprint('Error: Strategy unavailable. Select a different scheme element')    

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
        return ModelPartIO(modelpart)
    

########################################################## CHUNG, OOI BENCHMARKS
    
nodeplotter = 0  
  
if benchmark_number in listDISCONT:
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
        # Creating necessary directories
        main_path = os.getcwd()
        [post_path, data_and_results, graphs_path, MPI_results] = procedures.CreateDirectories(str(main_path), str(DEM_parameters.problem_name))
        demio         = DEM_procedures.DEMIo(DEM_parameters, post_path)
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

        # Constructing utilities objects
        creator_destructor = ParticleCreatorDestructor()
        dem_fem_search = DEM_FEM_Search()
        
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


        # Creating a solver object and set the search strategy
        solver = SolverStrategy.ExplicitStrategy(spheres_model_part, rigid_face_model_part, cluster_model_part, DEM_inlet_model_part, creator_destructor, dem_fem_search, scheme, DEM_parameters, procedures)
        
        procedures.AddAllVariablesInAllModelParts(solver, scheme, spheres_model_part, cluster_model_part, DEM_inlet_model_part, rigid_face_model_part, DEM_parameters)

        problem_name                = 'benchmark' + str(benchmark_number)
        spheres_mp_filename         = problem_name + "DEM"
        DEM_parameters.problem_name = problem_name

        os.chdir(main_path)
        model_part_io_spheres = model_part_reader(spheres_mp_filename)

        if (hasattr(DEM_parameters, "do_not_perform_initial_partition") and DEM_parameters.do_not_perform_initial_partition == 1):
            pass
        else:
            parallelutils.PerformInitialPartition(model_part_io_spheres)

        os.chdir(main_path)
        [model_part_io_spheres, spheres_model_part, MPICommSetup] = parallelutils.SetCommunicator(spheres_model_part, model_part_io_spheres, spheres_mp_filename)

        model_part_io_spheres.ReadModelPart(spheres_model_part) #########################3??????????????????????

        ###################################################### CHUNG, OOI BENCHMARKS

        benchmark.set_initial_data(spheres_model_part, rigid_face_model_part, iteration, number_of_points_in_the_graphic, coeff_of_restitution_iteration)
        
        ###################################################### CHUNG, OOI BENCHMARKS

        max_node_Id = creator_destructor.FindMaxNodeIdInModelPart(spheres_model_part)
        max_elem_Id = creator_destructor.FindMaxElementIdInModelPart(spheres_model_part)
        old_max_elem_Id_spheres = max_elem_Id
        max_cond_Id = creator_destructor.FindMaxElementIdInModelPart(spheres_model_part)        
        rigidFace_mp_filename = DEM_parameters.problem_name + "DEM_FEM_boundary"
        model_part_io_fem = model_part_reader(rigidFace_mp_filename,max_node_Id+1, max_elem_Id+1, max_cond_Id+1)
        model_part_io_fem.ReadModelPart(rigid_face_model_part)

        max_node_Id = creator_destructor.FindMaxNodeIdInModelPart(rigid_face_model_part)
        max_elem_Id = creator_destructor.FindMaxElementIdInModelPart(rigid_face_model_part)
        max_cond_Id = creator_destructor.FindMaxElementIdInModelPart(rigid_face_model_part)        
        clusters_mp_filename = 'benchmark' + "DEM_Clusters"
        model_part_io_clusters = model_part_reader(clusters_mp_filename,max_node_Id+1, max_elem_Id+1, max_cond_Id+1)
        model_part_io_clusters.ReadModelPart(cluster_model_part)

        max_node_Id = creator_destructor.FindMaxNodeIdInModelPart(cluster_model_part)
        max_elem_Id = creator_destructor.FindMaxElementIdInModelPart(cluster_model_part)
        max_cond_Id = creator_destructor.FindMaxElementIdInModelPart(cluster_model_part)        
        DEM_Inlet_filename = 'benchmark' + "DEM_Inlet"  
        model_part_io_demInlet = model_part_reader(DEM_Inlet_filename,max_node_Id+1, max_elem_Id+1, max_cond_Id+1)
        model_part_io_demInlet.ReadModelPart(DEM_inlet_model_part)

        # Setting up the buffer size
        procedures.SetUpBufferSizeInAllModelParts(spheres_model_part, 1, cluster_model_part, 1, DEM_inlet_model_part, 1, rigid_face_model_part, 1)
        # Adding dofs
        solver.AddDofs(spheres_model_part)
        solver.AddDofs(cluster_model_part)
        solver.AddDofs(DEM_inlet_model_part)

        os.chdir(main_path)

        KRATOSprint("Initializing Problem....")

        demio.Initialize(DEM_parameters)

        os.chdir(post_path)
        demio.InitializeMesh(spheres_model_part,rigid_face_model_part, cluster_model_part, solver.contact_model_part, mapping_model_part)

        # Perform a partition to balance the problem
        solver.search_strategy = parallelutils.GetSearchStrategy(solver, spheres_model_part)
        solver.BeforeInitialize()
        parallelutils.Repart(spheres_model_part)

        # Setting up the BoundingBox
        bounding_box_time_limits = []
        if (DEM_parameters.BoundingBoxOption == "ON"):
            procedures.SetBoundingBox(spheres_model_part, cluster_model_part, rigid_face_model_part, creator_destructor)
            bounding_box_time_limits = [solver.bounding_box_start_time, solver.bounding_box_stop_time]
        
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

        #Constructing a model part for the DEM inlet. It contains the DEM elements to be released during the simulation  
        #Initializing the DEM solver must be done before creating the DEM Inlet, because the Inlet configures itself according to some options of the DEM model part
        if (DEM_parameters.dem_inlet_option):
            #Constructing the inlet and initializing it (must be done AFTER the spheres_model_part Initialize)    
            DEM_inlet = DEM_Inlet(DEM_inlet_model_part)    
            DEM_inlet.InitializeDEM_Inlet(spheres_model_part, creator_destructor, solver.continuum_type)

        DEMFEMProcedures = DEM_procedures.DEMFEMProcedures(DEM_parameters, graphs_path, spheres_model_part, rigid_face_model_part)

        if (hasattr(DEM_parameters, "EnergyCalculationOption")):
            if (DEM_parameters.EnergyCalculationOption):
                os.chdir(graphs_path)
                DEMEnergyCalculator = DEM_procedures.DEMEnergyCalculator(DEM_parameters, spheres_model_part, cluster_model_part, "EnergyPlot.grf")

        materialTest.Initialize(DEM_parameters, procedures, solver, graphs_path, post_path, spheres_model_part, rigid_face_model_part)

        KRATOSprint("Initialization Complete" + "\n")

        step           = 0
        time           = 0.0
        time_old_print = 0.0

        report.Prepare(timer,DEM_parameters.ControlTime)

        first_print  = True; index_5 = 1; index_10  = 1; index_50  = 1; control = 0.0

        coordination_number = procedures.ModelData(spheres_model_part, solver.contact_model_part, solver)
        KRATOSprint ("Coordination Number: " + str(coordination_number) + "\n")

        materialTest.PrintChart()
        materialTest.PrepareDataForGraph()

        post_utils = DEM_procedures.PostUtils(DEM_parameters, spheres_model_part)

        step = 0

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

            DEMFEMProcedures.UpdateTimeInModelParts(spheres_model_part, rigid_face_model_part, cluster_model_part, time,dt,step) 

            benchmark.ApplyNodalRotation(time, dt, spheres_model_part)            
            #### SOLVE #########################################
            
            #benchmark.ApplyNodalRotation(time, dt, spheres_model_part)
            solver.Solve()
            #benchmark.ApplyNodalRotation(time, dt, spheres_model_part)
            
            DEMFEMProcedures.MoveAllMeshes(rigid_face_model_part, spheres_model_part, DEM_inlet_model_part, time, dt)

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
            DEMFEMProcedures.PrintGraph(time)
            benchmark.generate_graph_points(spheres_model_part, rigid_face_model_part, time, output_time_step, dt)

            #### GiD IO ##########################################
            time_to_print = time - time_old_print

            if (output_time_step - time_to_print < 1e-2*dt):            

                '''KRATOSprint("*******************  PRINTING RESULTS FOR GID  ***************************")
                KRATOSprint("                        ("+ str(spheres_model_part.NumberOfElements(0)) + " elements)")
                KRATOSprint("                        ("+ str(spheres_model_part.NumberOfNodes(0)) + " nodes)")
                KRATOSprint("")'''

                os.chdir(data_and_results)

                demio.PrintMultifileLists(time, post_path)

                os.chdir(post_path)

                solver.PrepareElementsForPrinting()
                if (DEM_parameters.ContactMeshOption == "ON"):
                    solver.PrepareContactElementsForPrinting()

                demio.PrintResults(spheres_model_part, rigid_face_model_part, cluster_model_part, solver.contact_model_part, mapping_model_part, creator_destructor, dem_fem_search, time, bounding_box_time_limits)
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

        benchmark.get_final_data(spheres_model_part, rigid_face_model_part)

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

# Freeing up memory
objects_to_destroy = [demio, procedures, creator_destructor, dem_fem_search, solver, DEMFEMProcedures, post_utils, 
                      cluster_model_part, rigid_face_model_part, spheres_model_part, DEM_inlet_model_part, mapping_model_part]

if (DEM_parameters.dem_inlet_option):
    objects_to_destroy.append(DEM_inlet)

del objects_to_destroy[:]

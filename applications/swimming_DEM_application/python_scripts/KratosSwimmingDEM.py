# This script contains an algorithm that models fluid-particle interaction.
# It combines two parts: a FEM model for the fluid and a DEM model for the particles.
# It has been conceived by adding the DEM part and the interaction on top of an original fluid-only script (see kratos/applications/FluidDynamicswApplication)
# Some parts of the original fluid script have been kept practically untouched and are clearly marked.
# Whenever a minor modification has been made on one of these parts, the corresponding line is indicated with a comment: # MOD.

from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Python imports
import os
import sys

class Logger(object):
    def __init__(self):
        self.terminal = sys.stdout
        self.path_to_console_out_file = 'console_output.txt'
        self.log = open(self.path_to_console_out_file, "a")

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)

    def flush(self):
        #this flush method is needed for python 3 compatibility.
        #this handles the flush command by doing nothing.
        #you might want to specify some extra behavior here.
        pass

sys.stdout = Logger()

# Kratos
from KratosMultiphysics import *
from KratosMultiphysics.ExternalSolversApplication import *
from KratosMultiphysics.MeshingApplication import *
from KratosMultiphysics.DEMApplication import *
from KratosMultiphysics.FluidDynamicsApplication import *
from KratosMultiphysics.SolidMechanicsApplication import *
from KratosMultiphysics.IncompressibleFluidApplication import *
from KratosMultiphysics.SwimmingDEMApplication import *

sys.path.insert(0,'')
# import the configuration data as read from the GiD
import define_output

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

class Solution:
    import swimming_DEM_algorithm
    def __init__(self, algorithm = swimming_DEM_algorithm, varying_parameters = dict()):
        import DEM_explicit_solver_var as DEM_parameters
        import ProjectParameters as pp
        self.main_path = os.getcwd()
        self.pp = pp
        self.pp.main_path = os.getcwd()
        self.pp.CFD_DEM = DEM_parameters
        self.alg = algorithm.Algorithm(self.pp)
        self.alg.CreateParts()
        self.alg.SetCustomBetaParamters(varying_parameters)

    def Run(self):
        import math
        import swimming_DEM_procedures as swim_proc
        import CFD_DEM_coupling
        import embedded
        import swimming_DEM_algorithm
        # import the configuration data as read from the GiD
        import define_output

        run_code = self.alg.GetRunCode()

        # Moving to the recently created folder
        os.chdir(self.main_path)
        [post_path, data_and_results, graphs_path, MPI_results] = self.alg.procedures.CreateDirectories(str(self.main_path), str(self.pp.CFD_DEM.problem_name))
        swim_proc.CopyInputFilesIntoFolder(self.main_path, post_path)

        self.alg.AddExtraVariables()

        [post_path, data_and_results, graphs_path, MPI_results] = self.alg.procedures.CreateDirectories(str(self.main_path), str(self.pp.CFD_DEM.problem_name))

        os.chdir(self.main_path)

        # Initialize GiD-IO
        self.demio = DEM_procedures.DEMIo(self.pp.CFD_DEM, post_path)
        self.demio.AddGlobalVariables()
        self.demio.AddSpheresVariables()
        self.demio.AddFEMBoundaryVariables()
        self.demio.AddClusterVariables()
        self.demio.AddContactVariables()
        # MPI
        self.demio.AddMpiVariables()

        self.demio.Configure(self.pp.CFD_DEM.problem_name,
                             self.pp.CFD_DEM.OutputFileType,
                             self.pp.CFD_DEM.Multifile,
                             self.pp.CFD_DEM.ContactMeshOption)

        self.demio.SetOutputName(self.pp.CFD_DEM.problem_name)
        os.chdir(post_path)
        self.demio.InitializeMesh(self.alg.all_model_parts)

        #Setting up the BoundingBox
        bounding_box_time_limits = []
        if (self.pp.CFD_DEM.BoundingBoxOption == "ON"):
            self.alg.procedures.SetBoundingBox(self.alg.spheres_model_part, self.alg.cluster_model_part, self.alg.rigid_face_model_part, self.alg.creator_destructor)
            bounding_box_time_limits = [self.alg.solver.bounding_box_start_time, self.alg.solver.bounding_box_stop_time]

        # Creating the fluid solver
        SolverSettings = self.pp.FluidSolverConfiguration
        # solver_module = import_solver(SolverSettings)
        fluid_model_part = self.alg.all_model_parts.Get('FluidPart')
        mixed_model_part = self.alg.all_model_parts.Get('MixedPart')
        self.alg.fluid_solver = self.alg.solver_module.CreateSolver(fluid_model_part, SolverSettings)

        Dt_DEM = self.pp.CFD_DEM.MaxTimeStep

        # reading the fluid part
        self.alg.Initialize()

        # activate turbulence model
        self.alg.ActivateTurbulenceModel()

        # constructing a model part for the DEM inlet. it contains the DEM elements to be released during the simulation
        # Initializing the DEM solver must be done before creating the DEM Inlet, because the Inlet configures itself according to some options of the DEM model part

        if not self.pp.VolumeOutput:
            cut_list = define_output.DefineCutPlanes()
            gid_io.define_cuts(fluid_model_part, cut_list)

        # gid_io.initialize_results(fluid_model_part) # MOD.

        import swimming_DEM_gid_output
        swimming_DEM_gid_io = swimming_DEM_gid_output.SwimmingDEMGiDOutput(self.pp.problem_name,
                                                                           self.pp.VolumeOutput,
                                                                           self.pp.GiDPostMode,
                                                                           self.pp.GiDMultiFileFlag,
                                                                           self.pp.GiDWriteMeshFlag,
                                                                           self.pp.GiDWriteConditionsFlag)

        swimming_DEM_gid_io.initialize_swimming_DEM_results(self.alg.spheres_model_part, self.alg.cluster_model_part, self.alg.rigid_face_model_part, mixed_model_part)

        # define the drag computation list
        drag_list = define_output.DefineDragList()
        drag_file_output_list = []
        for it in drag_list:
            f = open(it[1], 'w')
            drag_file_output_list.append(f)
            tmp = "#Drag for group " + it[1] + "\n"
            f.write(tmp)
            tmp = "time RX RY RZ"
            f.write(tmp)
            f.flush()

        print(drag_file_output_list)


        def PrintDrag(drag_list, drag_file_output_list, fluid_model_part, time):
            i = 0
            for it in drag_list:
                print(it[0])
                nodes = fluid_model_part.GetNodes(it[0])
                drag = Vector(3)
                drag[0] = 0.0
                drag[1] = 0.0
                drag[2] = 0.0
                for node in nodes:
                    reaction = node.GetSolutionStepValue(REACTION, 0)
                    drag[0] += reaction[0]
                    drag[1] += reaction[1]
                    drag[2] += reaction[2]

                output = str(time) + " " + str(drag[0]) + " " + str(drag[1]) + " " + str(drag[2]) + "\n"
                # print drag_file_output_list[i]
                # print output
                drag_file_output_list[i].write(output)
                drag_file_output_list[i].flush()
                i = i + 1


        # preparing output of point graphs
        #import point_graph_printer

        #output_nodes_list = define_output.DefineOutputPoints()
        #graph_printer = point_graph_printer.PrintGraphPrinter(
            #output_nodes_list,
            #fluid_model_part,
            #variables_dictionary,
            #domain_size)

        # setting fluid's body force to the same as DEM's
        if self.pp.CFD_DEM.body_force_on_fluid_option:

            for node in fluid_model_part.Nodes:
                node.SetSolutionStepValue(BODY_FORCE_X, 0, self.pp.CFD_DEM.GravityX)
                node.SetSolutionStepValue(BODY_FORCE_Y, 0, self.pp.CFD_DEM.GravityY)
                node.SetSolutionStepValue(BODY_FORCE_Z, 0, self.pp.CFD_DEM.GravityZ)

        # coarse-graining: applying changes to the physical properties of the model to adjust for
        # the similarity transformation if required (fluid effects only).
        swim_proc.ApplySimilarityTransformations(fluid_model_part, self.pp.CFD_DEM.similarity_transformation_type , self.pp.CFD_DEM.model_over_real_diameter_factor )

        # creating a Post Utils object that executes several post-related tasks
        post_utils = swim_proc.PostUtils(swimming_DEM_gid_io,
                                                       self.pp,
                                                       fluid_model_part,
                                                       self.alg.spheres_model_part,
                                                       self.alg.cluster_model_part,
                                                       self.alg.rigid_face_model_part,
                                                       mixed_model_part)

        # creating an IOTools object to perform other printing tasks
        io_tools = swim_proc.IOTools(self.pp)

        # creating a projection module for the fluid-DEM coupling
        h_min = 0.01
        n_balls = 1
        fluid_volume = 10
        self.pp.CFD_DEM.n_particles_in_depth = int(math.sqrt(n_balls / fluid_volume)) # only relevant in 2D problems
        # creating a physical calculations module to analyse the DEM model_part
        dem_physics_calculator = SphericElementGlobalPhysicsCalculator(self.alg.spheres_model_part)

        if self.pp.CFD_DEM.coupling_level_type:

            if self.pp.CFD_DEM.meso_scale_length <= 0.0 and self.alg.spheres_model_part.NumberOfElements(0) > 0:
                biggest_size = 2 * dem_physics_calculator.CalculateMaxNodalVariable(self.alg.spheres_model_part, RADIUS)
                self.pp.CFD_DEM.meso_scale_length  = 20 * biggest_size

            elif self.alg.spheres_model_part.NumberOfElements(0) == 0:
                self.pp.CFD_DEM.meso_scale_length  = 1.0

            field_utility = self.alg.GetFieldUtility()

            self.alg.projection_module = CFD_DEM_coupling.ProjectionModule(fluid_model_part, self.alg.spheres_model_part, self.alg.rigid_face_model_part, self.pp.domain_size, self.pp, field_utility)
            self.alg.projection_module.UpdateDatabase(h_min)

        # creating a custom functions calculator for the implementation of additional custom functions
        custom_functions_tool = swim_proc.FunctionsCalculator(self.pp)

        # creating a derivative recovery tool to calculate the necessary derivatives from the fluid solution (gradient, laplacian, material acceleration...)
        derivative_recovery_tool = DerivativeRecoveryTool3D(fluid_model_part)

        # creating a basset_force tool to perform the operations associated with the calculation of this force along the path of each particle
        if self.pp.CFD_DEM.basset_force_type > 0:
            basset_force_tool = swim_proc.BassetForceTools()

        # creating a stationarity assessment tool
        stationarity_tool = swim_proc.StationarityAssessmentTool(self.pp.CFD_DEM.max_pressure_variation_rate_tol , custom_functions_tool)

        # creating a debug tool
        dem_volume_tool = swim_proc.ProjectionDebugUtils(self.pp.CFD_DEM.fluid_domain_volume, fluid_model_part, self.alg.spheres_model_part, custom_functions_tool)

        # creating a distance calculation process for the embedded technology
        # (used to calculate elemental distances defining the structure embedded in the fluid mesh)
        if self.pp.CFD_DEM.embedded_option:
            calculate_distance_process = CalculateSignedDistanceTo3DSkinProcess(self.alg.rigid_face_model_part, fluid_model_part)
            calculate_distance_process.Execute()

        self.alg.KRATOSprint("Initialization Complete" + "\n")

        step           = 0
        time           = self.pp.Start_time
        Dt             = self.pp.Dt
        out            = Dt
        Nsteps         = self.pp.nsteps
        final_time     = self.pp.CFD_DEM.FinalTime
        output_time    = self.pp.CFD_DEM.OutputTimeStep

        self.alg.report.Prepare(self.alg.timer, self.pp.CFD_DEM.ControlTime)

        first_print = True; index_5 = 1; index_10 = 1; index_50 = 1; control = 0.0

        if (self.pp.CFD_DEM.ModelDataInfo == "ON"):
            os.chdir(data_and_results)
            if (self.pp.CFD_DEM.ContactMeshOption == "ON"):
                (coordination_number) = self.alg.procedures.ModelData(self.alg.spheres_model_part, self.alg.solver) # Calculates the mean number of neighbours the mean radius, etc..
                self.alg.KRATOSprint ("Coordination Number: " + str(coordination_number) + "\n")
                os.chdir(self.main_path)
            else:
                self.alg.KRATOSprint("Activate Contact Mesh for ModelData information")

        if self.pp.CFD_DEM.flow_in_porous_medium_option:
            fluid_frac_util = swim_proc.FluidFractionFieldUtility(fluid_model_part, self.pp.CFD_DEM.min_fluid_fraction )

            for field in self.pp.fluid_fraction_fields:
                fluid_frac_util.AppendLinearField(field)

            fluid_frac_util.AddFluidFractionField()

        if self.pp.CFD_DEM.flow_in_porous_DEM_medium_option:
            swim_proc.FixModelPart(self.alg.spheres_model_part)

        # choosing the directory in which we want to work (print to)

        os.chdir(post_path)

        def yield_DEM_time(current_time, current_time_plus_increment, delta_time):
            current_time += delta_time

            tolerance = 0.0001
            while current_time < (current_time_plus_increment - tolerance * delta_time):
                yield current_time
                current_time += delta_time

            current_time = current_time_plus_increment
            yield current_time

        ######################################################################################################################################

        #                      I N I T I A L I Z I N G    T I M E    L O O P     ...   ( M I X E D    F L U I D / D E M    B L O C K )

        ######################################################################################################################################

        # setting up loop counters: Counter(steps_per_tick_step, initial_step, active_or_inactive_boolean, dead_or_not)
        fluid_solve_counter          = self.alg.GetFluidSolveCounter()
        embedded_counter             = self.alg.GetEmbeddedCounter()
        DEM_to_fluid_counter         = self.alg.GetBackwardCouplingCounter()
        derivative_recovery_counter  = self.alg.GetRecoveryCounter()
        stationarity_counter         = self.alg.GetStationarityCounter()
        print_counter                = self.alg.GetPrintCounter()
        debug_info_counter           = self.alg.GetDebugInfo()
        particles_results_counter    = self.alg.GetParticlesResultsCounter()
        quadrature_counter           = self.alg.HistoryForceQuadratureCounter()
        mat_deriv_averager           = swim_proc.Averager(1, 3)
        laplacian_averager           = swim_proc.Averager(1, 3)

        ##############################################################################
        #                                                                            #
        #    MAIN LOOP                                                               #
        #                                                                            #
        ##############################################################################

        DEM_step     = 0      # necessary to get a good random insertion of particles   # relevant to the stationarity assessment tool
        time_dem     = 0.0
        Dt_DEM       = self.alg.spheres_model_part.ProcessInfo.GetValue(DELTA_TIME)
        self.alg.rigid_face_model_part.ProcessInfo[DELTA_TIME] = Dt_DEM
        self.alg.cluster_model_part.ProcessInfo[DELTA_TIME] = Dt_DEM
        stationarity = False

        self.alg.report.total_steps_expected = int(self.pp.CFD_DEM.FinalTime / Dt_DEM)

        self.alg.KRATOSprint(self.alg.report.BeginReport(self.alg.timer))

        mesh_motion = DEMFEMUtilities()

        # creating a Post Utils object that executes several post-related tasks
        post_utils_DEM = DEM_procedures.PostUtils(self.pp.CFD_DEM, self.alg.spheres_model_part)

        swim_proc.InitializeVariablesWithNonZeroValues(fluid_model_part, self.alg.spheres_model_part, self.pp) # otherwise variables are set to 0 by default

        self.alg.SetUpResultsDatabase()

        # ANALYTICS BEGIN
        self.pp.CFD_DEM.perform_analytics_option = False

        if self.pp.CFD_DEM.perform_analytics_option:
            import analytics
            variables_to_measure = [PRESSURE]
            steps_between_measurements = 100
            gauge = analytics.Gauge(fluid_model_part, Dt, final_time, variables_to_measure, steps_between_measurements)
            point_coors = [0.0, 0.0, 0.01]
            target_node = swim_proc.FindClosestNode(fluid_model_part, point_coors)
            target_id = target_node.Id
            print(target_node.X, target_node.Y, target_node.Z)
            print(target_id)
            def condition(node):
                return node.Id == target_id

            gauge.ConstructArrayOfNodes(condition)
            print(gauge.variables)
            print_analytics_counter = swim_proc.Counter( 5 * steps_between_measurements, 1, 1)
        # ANALYTICS END

        import derivative_recovery.derivative_recovery_strategy as derivative_recoverer
        recovery = derivative_recoverer.DerivativeRecoveryStrategy(self.pp, fluid_model_part, derivative_recovery_tool, custom_functions_tool)

        N_steps = int(final_time / Dt_DEM) + 20

        if self.pp.CFD_DEM.basset_force_type > 0:
            basset_force_tool.FillDaitcheVectors(N_steps, self.pp.CFD_DEM.quadrature_order, self.pp.CFD_DEM.time_steps_per_quadrature_step)
        if self.pp.CFD_DEM.basset_force_type >= 3 or self.pp.CFD_DEM.basset_force_type == 1:
            basset_force_tool.FillHinsbergVectors(self.alg.spheres_model_part, self.pp.CFD_DEM.number_of_exponentials, self.pp.CFD_DEM.number_of_quadrature_steps_in_window)

        self.alg.PerformZeroStepInitializations()

        post_utils.Writeresults(time)

        while (time <= final_time):

            time = time + Dt
            step += 1
            fluid_model_part.CloneTimeStep(time)
            self.alg.TellTime(time)


            if self.pp.CFD_DEM.coupling_scheme_type  == "UpdatedDEM":
                time_final_DEM_substepping = time + Dt

            else:
                time_final_DEM_substepping = time

            # calculating elemental distances defining the structure embedded in the fluid mesh
            if self.pp.CFD_DEM.embedded_option:
                calculate_distance_process.Execute()

            if embedded_counter.Tick():
                embedded.ApplyEmbeddedBCsToFluid(fluid_model_part)
                embedded.ApplyEmbeddedBCsToBalls(self.alg.spheres_model_part, self.pp.CFD_DEM)

            # solving the fluid part

            if step >= 3 and not stationarity:
                print("Solving Fluid... (", fluid_model_part.NumberOfElements(0), "elements )")
                sys.stdout.flush()

                if fluid_solve_counter.Tick():
                    self.alg.FluidSolve(time)

            # assessing stationarity

                if stationarity_counter.Tick():
                    print("Assessing Stationarity...")
                    stationarity = stationarity_tool.Assess(fluid_model_part)
                    sys.stdout.flush()

            # printing if required

            if particles_results_counter.Tick():
                # eliminating remote balls

                #if self.pp.dem.BoundingBoxOption == "ON":
                #    self.alg.creator_destructor.DestroyParticlesOutsideBoundingBox(self.alg.spheres_model_part)

                io_tools.PrintParticlesResults(self.alg.pp.variables_to_print_in_file, time, self.alg.spheres_model_part)
                graph_printer.PrintGraphs(time)
                PrintDrag(drag_list, drag_file_output_list, fluid_model_part, time)

            if output_time <= out and self.alg.pp.CFD_DEM.coupling_scheme_type == "UpdatedDEM":

                if self.pp.CFD_DEM.coupling_level_type > 0:
                    self.alg.projection_module.ComputePostProcessResults(self.alg.spheres_model_part.ProcessInfo)

                post_utils.Writeresults(time)
                out = 0

            # solving the DEM part
            derivative_recovery_counter.Switch(time > self.pp.CFD_DEM.interaction_start_time)

            if derivative_recovery_counter.Tick():
                recovery.Recover()

            print("Solving DEM... (", self.alg.spheres_model_part.NumberOfElements(0), "elements )")
            sys.stdout.flush()
            first_dem_iter = True

            for time_dem in yield_DEM_time(time_dem, time_final_DEM_substepping, Dt_DEM):
                DEM_step += 1   # this variable is necessary to get a good random insertion of particles
                self.alg.spheres_model_part.ProcessInfo[TIME_STEPS]    = DEM_step
                self.alg.rigid_face_model_part.ProcessInfo[TIME_STEPS] = DEM_step
                self.alg.cluster_model_part.ProcessInfo[TIME_STEPS]    = DEM_step

                self.alg.PerformInitialDEMStepOperations(time_dem)

                if time >= self.pp.CFD_DEM.interaction_start_time and self.pp.CFD_DEM.coupling_level_type and (self.pp.CFD_DEM.project_at_every_substep_option or first_dem_iter):

                    if self.pp.CFD_DEM.coupling_scheme_type == "UpdatedDEM":
                        self.alg.ApplyForwardCoupling()

                    else:
                        self.alg.ApplyForwardCoupling((time_final_DEM_substepping - time_dem) / Dt)

                        if self.alg.pp.CFD_DEM.IntegrationScheme == 'Hybrid_Bashforth':
                            self.alg.solver.Solve() # only advance in space
                            self.alg.ApplyForwardCouplingOfVelocityOnly(time_dem)
                        else:
                            if self.alg.pp.CFD_DEM.basset_force_type > 0:
                                node.SetSolutionStepValue(SLIP_VELOCITY_X, vx)
                                node.SetSolutionStepValue(SLIP_VELOCITY_Y, vy)

                        if quadrature_counter.Tick():
                            if self.pp.CFD_DEM.basset_force_type == 1 or self.pp.CFD_DEM.basset_force_type >= 3:
                                basset_force_tool.AppendIntegrandsWindow(self.alg.spheres_model_part)
                            elif self.pp.CFD_DEM.basset_force_type == 2:
                                basset_force_tool.AppendIntegrands(self.alg.spheres_model_part)

                # performing the time integration of the DEM part

                self.alg.spheres_model_part.ProcessInfo[TIME]    = time_dem
                self.alg.rigid_face_model_part.ProcessInfo[TIME] = time_dem
                self.alg.cluster_model_part.ProcessInfo[TIME]    = time_dem

                if self.alg.pp.do_solve_dem:
                    self.alg.DEMSolve(time_dem)

                # Walls movement:
                mesh_motion.MoveAllMeshes(self.alg.rigid_face_model_part, time, Dt)
                mesh_motion.MoveAllMeshes(self.alg.spheres_model_part, time, Dt)
                mesh_motion.MoveAllMeshes(self.alg.DEM_inlet_model_part, time, Dt)

                #### TIME CONTROL ##################################

                # adding DEM elements by the inlet:
                if self.pp.CFD_DEM.dem_inlet_option:
                    self.alg.DEM_inlet.CreateElementsFromInletMesh(self.alg.spheres_model_part, self.alg.cluster_model_part, self.alg.creator_destructor)  # After solving, to make sure that neighbours are already set.

                if output_time <= out and self.pp.CFD_DEM.coupling_scheme_type == "UpdatedFluid":

                    if self.pp.CFD_DEM.coupling_level_type:
                        self.alg.projection_module.ComputePostProcessResults(self.alg.spheres_model_part.ProcessInfo)

                    post_utils.Writeresults(time_dem)
                    out = 0

                out = out + Dt_DEM
                first_dem_iter = False

            #### PRINTING GRAPHS ####
            os.chdir(graphs_path)
            # measuring mean velocities in a certain control volume (the 'velocity trap')
            if self.pp.CFD_DEM.VelocityTrapOption:
                post_utils_DEM.ComputeMeanVelocitiesinTrap("Average_Velocity.txt", time)

            os.chdir(post_path)

            # applying DEM-to-fluid coupling
            if DEM_to_fluid_counter.Tick() and time >= self.pp.CFD_DEM.interaction_start_time:
                self.alg.projection_module.ProjectFromParticles()

            # coupling checks (debugging)
            if debug_info_counter.Tick():
                dem_volume_tool.UpdateDataAndPrint(self.pp.CFD_DEM.fluid_domain_volume)

            # printing if required

            if particles_results_counter.Tick():
                io_tools.PrintParticlesResults(self.pp.variables_to_print_in_file, time, self.alg.spheres_model_part)
                graph_printer.PrintGraphs(time)
                PrintDrag(drag_list, drag_file_output_list, fluid_model_part, time)

        swimming_DEM_gid_io.finalize_results()

        self.alg.PerformFinalOperations(time_dem)

        for i in drag_file_output_list:
            i.close()

        self.alg.TellFinalSummary(step, time, DEM_step)

        return self.alg.GetReturnValue()

if __name__=="__main__":
    Solution().Run()

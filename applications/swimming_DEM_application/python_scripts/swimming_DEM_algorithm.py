from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

#TODO: test DEM bounding box

import os
import sys
import math
import time as timer
import weakref

from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *
from KratosMultiphysics.SwimmingDEMApplication import *

from DEM_procedures import KratosPrint as Say
import CFD_DEM_coupling
import swimming_DEM_procedures as SDP
import swimming_DEM_gid_output
import embedded
import variables_management as vars_man

try:
    import define_output  # MA: some GUI write this file, some others not!
except ImportError:
    pass

# Import MPI modules if needed. This way to do this is only valid when using OpenMPI. For other implementations of MPI it will not work.
if "OMPI_COMM_WORLD_SIZE" in os.environ:
    # Kratos MPI
    from KratosMultiphysics.MetisApplication import *
    from KratosMultiphysics.MPISearchApplication import *
    from KratosMultiphysics.mpi import *

    # DEM Application MPI
    import DEM_procedures_mpi as DEM_procedures
    import DEM_material_test_script_mpi as DEM_material_test_script
    Say('Running under MPI...........\n')
else:
    # DEM Application
    import DEM_procedures
    import DEM_material_test_script
    Say('Running under OpenMP........\n')

sys.path.insert(0,'')

class Logger(object):
    def __init__(self):
        self.terminal = sys.stdout
        self.console_output_file_name = 'console_output.txt'
        self.path_to_console_out_file = os.getcwd()
        self.path_to_console_out_file += '/' + self.console_output_file_name
        self.log = open(self.path_to_console_out_file, "a")

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)

    def flush(self):
        #this flush method is needed for python 3 compatibility.
        #this handles the flush command by doing nothing.
        #you might want to specify some extra behavior here.
        pass

class Algorithm(object):
    def __enter__ (self):
        # sys.stdout = Logger()
        return self

    def __exit__(self, exception_type, exception_value, traceback):
        pass

    def __init__(self, varying_parameters = Parameters("{}")):
        sys.stdout = Logger()
        self.StartTimer()
        self.main_path = os.getcwd()

        self.SetFluidAlgorithm()
        self.fluid_solution.coupling_algorithm = weakref.proxy(self)

        self.pp = self.fluid_solution.pp

        self.SetCouplingParameters(varying_parameters)

        self.SetDispersePhaseAlgorithm()
        self.disperse_phase_solution.coupling_algorithm = weakref.proxy(self)

        self.procedures = weakref.proxy(self.disperse_phase_solution.procedures)
        self.report = DEM_procedures.Report()

        # creating a basset_force tool to perform the operations associated with the calculation of this force along the path of each particle
        self.GetBassetForceTools()

    def SetFluidAlgorithm(self):
        import eulerian_fluid_ready_for_coupling
        self.fluid_solution = eulerian_fluid_ready_for_coupling.Solution()
        self.fluid_solution.main_path = self.main_path

    def SetDispersePhaseAlgorithm(self):
        import dem_main_script_ready_for_coupling as DEM_algorithm
        self.disperse_phase_solution = DEM_algorithm.Solution(self.pp)

    def ReadDispersePhaseAndCouplingParameters(self):

        with open(self.main_path + '/ProjectParametersDEM.json', 'r') as parameters_file:
            self.pp.CFD_DEM = Parameters(parameters_file.read())

        import dem_default_input_parameters
        dem_defaults = dem_default_input_parameters.GetDefaultInputParameters()

        import swimming_dem_default_input_parameters
        only_swimming_defaults = swimming_dem_default_input_parameters.GetDefaultInputParameters()

        for key in only_swimming_defaults.keys():
            dem_defaults.AddValue(key,only_swimming_defaults[key])

        self.pp.CFD_DEM.ValidateAndAssignDefaults(dem_defaults)

    def SetCouplingParameters(self, varying_parameters):

        # First, read the parameters generated from the interface
        self.ReadDispersePhaseAndCouplingParameters()

        # Second, set the default 'beta' parameters (candidates to be moved to the interface)
        self.SetBetaParameters()

        # Third, set the parameters fed to the particular case that you are running
        self.SetCustomBetaParameters(varying_parameters)

        # Finally adjust some of the paramters for consistency
        #   This function should be reduced to a minimum since,
        #   in principle, there should be no need to change the parameters
        self.SetDerivedParameters()

    def SetDerivedParameters(self):
        self.SetDoSolveDEMVariable()

    def SetAllModelParts(self):
        self.all_model_parts = weakref.proxy(self.disperse_phase_solution.all_model_parts)

        # defining a fluid model
        self.all_model_parts.Add(self.fluid_solution.fluid_model_part)

        self.fluid_model_part = self.fluid_solution.fluid_model_part

        # defining a model part for the mixed part
        self.all_model_parts.Add(ModelPart("MixedPart"))

        self.mixed_model_part = self.all_model_parts.Get('MixedPart')


    def StartTimer(self):
        self.timer = timer
        self.simulation_start_time = timer.time()

    def SetBetaParameters(self): # These are input parameters that have not yet been transferred to the interface
        # import the configuration data as read from the GiD
        ##############################################################################
        #                                                                            #
        #    INITIALIZE                                                              #
        #                                                                            #
        ##############################################################################

        #G
        self.pp.CFD_DEM.AddEmptyValue("do_solve_dem").SetBool(True)
        self.pp.CFD_DEM.AddEmptyValue("fluid_already_calculated").SetBool(False)
        self.pp.CFD_DEM.AddEmptyValue("recovery_echo_level").SetInt(1)
        self.pp.CFD_DEM.AddEmptyValue("gradient_calculation_type").SetInt(1)
        self.pp.CFD_DEM.AddEmptyValue("pressure_grad_recovery_type").SetInt(0)
        self.pp.CFD_DEM.AddEmptyValue("fluid_fraction_grad_type").SetInt(0)
        self.pp.CFD_DEM.AddEmptyValue("store_full_gradient_option").SetBool(False)
        self.pp.CFD_DEM.AddEmptyValue("store_fluid_pressure_option").SetBool(False)
        self.pp.CFD_DEM.AddEmptyValue("laplacian_calculation_type").SetInt(0)
        self.pp.CFD_DEM.AddEmptyValue("do_search_neighbours").SetBool(True)
        self.pp.CFD_DEM.AddEmptyValue("faxen_terms_type").SetInt(0)
        self.pp.CFD_DEM.AddEmptyValue("material_acceleration_calculation_type").SetInt(1)
        self.pp.CFD_DEM.AddEmptyValue("faxen_force_type").SetInt(0)
        self.pp.CFD_DEM.AddEmptyValue("vorticity_calculation_type").SetInt(5)
        self.pp.CFD_DEM.AddEmptyValue("print_FLUID_VEL_PROJECTED_RATE_option").SetBool(False)
        self.pp.CFD_DEM.AddEmptyValue("print_MATERIAL_FLUID_ACCEL_PROJECTED_option").SetBool(True)
        self.pp.CFD_DEM.AddEmptyValue("basset_force_type").SetInt(0)
        self.pp.CFD_DEM.AddEmptyValue("print_BASSET_FORCE_option").SetBool(True)
        self.pp.CFD_DEM.AddEmptyValue("basset_force_integration_type").SetInt(2)
        self.pp.CFD_DEM.AddEmptyValue("n_init_basset_steps").SetInt(0)
        self.pp.CFD_DEM.AddEmptyValue("time_steps_per_quadrature_step").SetInt(1)
        self.pp.CFD_DEM.AddEmptyValue("delta_time_quadrature").SetDouble(
            self.pp.CFD_DEM["time_steps_per_quadrature_step"].GetInt() * self.pp.CFD_DEM["MaxTimeStep"].GetDouble())
        self.pp.CFD_DEM.AddEmptyValue("quadrature_order").SetInt(2)
        self.pp.CFD_DEM.AddEmptyValue("time_window").SetDouble(0.8)
        self.pp.CFD_DEM.AddEmptyValue("number_of_exponentials").SetInt(2)
        self.pp.CFD_DEM.AddEmptyValue("number_of_quadrature_steps_in_window").SetInt(
            int(self.pp.CFD_DEM["time_window"].GetDouble() / self.pp.CFD_DEM["delta_time_quadrature"].GetDouble()))
        self.pp.CFD_DEM.AddEmptyValue("do_impose_flow_from_field_option").SetBool(False)
        self.pp.CFD_DEM.AddEmptyValue("print_MATERIAL_ACCELERATION_option").SetBool(True)
        self.pp.CFD_DEM.AddEmptyValue("print_FLUID_ACCEL_FOLLOWING_PARTICLE_PROJECTED_option").SetBool(False)
        self.pp.CFD_DEM.AddEmptyValue("print_VORTICITY_option").SetBool(True)
        self.pp.CFD_DEM.AddEmptyValue("print_MATERIAL_ACCELERATION_option").SetBool(True)
        self.pp.CFD_DEM.AddEmptyValue("print_VISCOSITY_option").SetBool(False)
        self.pp.CFD_DEM.AddEmptyValue("print_VELOCITY_GRADIENT_option").SetBool(True)
        self.pp.CFD_DEM.AddEmptyValue("print_DISPERSE_FRACTION_option").SetBool(False)
        self.pp.CFD_DEM.AddEmptyValue("print_FLUID_FRACTION_GRADIENT_option").SetBool(False)
        self.pp.CFD_DEM.AddEmptyValue("print_FLUID_FRACTION_GRADIENT_PROJECTED_option").SetBool(False)
        self.pp.CFD_DEM.AddEmptyValue("calculate_diffusivity_option").SetBool(False)
        self.pp.CFD_DEM.AddEmptyValue("print_CONDUCTIVITY_option").SetBool(False)
        self.pp.CFD_DEM.AddEmptyValue("filter_velocity_option").SetBool(False)
        self.pp.CFD_DEM.AddEmptyValue("print_PARTICLE_VEL_option").SetBool(False)
        self.pp.CFD_DEM.AddEmptyValue("apply_time_filter_to_fluid_fraction_option").SetBool(False)
        self.pp.CFD_DEM.AddEmptyValue("full_particle_history_watcher").SetString("Empty")
        self.pp.CFD_DEM.AddEmptyValue("prerun_fluid_file_name").SetString("")
        self.pp.CFD_DEM.AddEmptyValue("frame_of_reference_type").SetInt(0)
        self.pp.CFD_DEM.AddEmptyValue("angular_velocity_of_frame_X").SetDouble(0.0)
        self.pp.CFD_DEM.AddEmptyValue("angular_velocity_of_frame_Y").SetDouble(0.0)
        self.pp.CFD_DEM.AddEmptyValue("angular_velocity_of_frame_Z").SetDouble(0.0)
        self.pp.CFD_DEM.AddEmptyValue("angular_velocity_of_frame_old_X").SetDouble(0.0)
        self.pp.CFD_DEM.AddEmptyValue("angular_velocity_of_frame_old_Y").SetDouble(0.0)
        self.pp.CFD_DEM.AddEmptyValue("angular_velocity_of_frame_old_Z").SetDouble(0.0)
        self.pp.CFD_DEM.AddEmptyValue("acceleration_of_frame_origin_X").SetDouble(0.0)
        self.pp.CFD_DEM.AddEmptyValue("acceleration_of_frame_origin_Y").SetDouble(0.0)
        self.pp.CFD_DEM.AddEmptyValue("acceleration_of_frame_origin_Z").SetDouble(0.0)
        self.pp.CFD_DEM.AddEmptyValue("angular_acceleration_of_frame_X").SetDouble(0.0)
        self.pp.CFD_DEM.AddEmptyValue("angular_acceleration_of_frame_Y").SetDouble(0.0)
        self.pp.CFD_DEM.AddEmptyValue("angular_acceleration_of_frame_Z").SetDouble(0.0)
        self.pp.CFD_DEM.AddEmptyValue("ALE_option").SetBool(False)
        self.pp.CFD_DEM.AddEmptyValue("frame_rotation_axis_initial_point").SetVector(Vector([0., 0., 0.]))
        self.pp.CFD_DEM.AddEmptyValue("frame_rotation_axis_final_point").SetVector(Vector([0., 0., 1.]))
        self.pp.CFD_DEM.AddEmptyValue("angular_velocity_magnitude").SetDouble(1.0)
        self.pp.CFD_DEM.print_DISPERSE_FRACTION_option = False
        self.pp.CFD_DEM.print_steps_per_plot_step = 1
        self.pp.CFD_DEM.PostCationConcentration = False
        # Making the fluid step an exact multiple of the DEM step
        self.pp.Dt = int(self.pp.Dt / self.pp.CFD_DEM["MaxTimeStep"].GetDouble()) * self.pp.CFD_DEM["MaxTimeStep"].GetDouble()
        self.output_time = int(self.pp.CFD_DEM["OutputTimeStep"].GetDouble() / self.pp.CFD_DEM["MaxTimeStep"].GetDouble()) * self.pp.CFD_DEM["MaxTimeStep"].GetDouble()
        Say('Dt_DEM', self.pp.CFD_DEM["MaxTimeStep"].GetDouble())
        Say('self.pp.Dt', self.pp.Dt)
        Say('self.output_time', self.output_time)
        self.pp.viscosity_modification_type = 0.0
        self.domain_size = 3
        self.pp.type_of_inlet = 'VelocityImposed' # 'VelocityImposed' or 'ForceImposed'
        self.pp.force = Vector(3)
        self.pp.force[0] = 0
        self.pp.force[1] = 0
        self.pp.force[2] = 1e-10

        # defining and adding imposed porosity fields
        self.pp.fluid_fraction_fields = []
        field1 = SDP.FluidFractionFieldUtility.LinearField(0.0,
                                                          [0.0, 0.0, 0.0],
                                                          [-1.0, -1.0, 0.15],
                                                          [1.0, 1.0, 0.3])
        from math import pi
        self.pp.CFD_DEM.AddEmptyValue("fluid_domain_volume").SetDouble(0.5 ** 2 * 2 * pi) # write down the volume you know it has

        self.pp.fluid_fraction_fields.append(field1)

    def SetDoSolveDEMVariable(self):
        self.do_solve_dem = self.pp.CFD_DEM["do_solve_dem"].GetBool()

        if self.pp.CFD_DEM["flow_in_porous_DEM_medium_option"].GetBool():
            self.do_solve_dem = False

    def SetCustomBetaParameters(self, custom_parameters): # this method is ugly. The way to go is to have all input parameters as a dictionary
        custom_parameters.ValidateAndAssignDefaults(self.pp.CFD_DEM)
        self.pp.CFD_DEM = custom_parameters
        # TO DO: remove next lines as info is taken from Parameters object everywhere
        # var_names = [k for k in dictionary.keys()]
        # var_values = [k for k in dictionary.values()]
        # for name, value in zip(var_names, var_values):
        #     self.pp.CFD_DEM.__setitem__(name, value)

    def Run(self):
        self.Initialize()
        self.RunMainTemporalLoop()
        self.Finalize()

        return self.GetReturnValue()

    def SetUpResultsDatabase(self):
        pass

    def ReadDispersePhaseModelParts(self, starting_node_Id = 0, starting_elem_Id = 0, starting_cond_Id = 0):
        fluid_mp = self.fluid_model_part
        max_node_Id = self.disperse_phase_solution.creator_destructor.FindMaxNodeIdInModelPart(fluid_mp)
        max_elem_Id = self.disperse_phase_solution.creator_destructor.FindMaxElementIdInModelPart(fluid_mp)
        max_cond_Id = self.disperse_phase_solution.creator_destructor.FindMaxConditionIdInModelPart(fluid_mp)
        self.disperse_phase_solution.BaseReadModelParts(max_node_Id, max_elem_Id, max_cond_Id)

    def Initialize(self):

        Say('\nInitializing Problem...\n')

        self.run_code = self.GetRunCode()

        # Moving to the recently created folder
        os.chdir(self.main_path)
        [self.post_path, data_and_results, self.graphs_path, MPI_results] = self.procedures.CreateDirectories(str(self.main_path), str(self.pp.CFD_DEM["problem_name"].GetString()), self.run_code)
        SDP.CopyInputFilesIntoFolder(self.main_path, self.post_path)

        #self.mixed_model_part = self.all_model_parts.Get('MixedPart')

        vars_man.ConstructListsOfVariables(self.pp)

        self.FluidInitialize()
        self.DispersePhaseInitialize()

        self.SetAllModelParts()

        self.SetCutsOutput()

        self.swimming_DEM_gid_io = swimming_DEM_gid_output.SwimmingDEMGiDOutput(self.pp.problem_name,
                                                                                self.pp.VolumeOutput,
                                                                                self.pp.GiDPostMode,
                                                                                self.pp.GiDMultiFileFlag,
                                                                                self.pp.GiDWriteMeshFlag,
                                                                                self.pp.GiDWriteConditionsFlag)

        self.swimming_DEM_gid_io.initialize_swimming_DEM_results(self.disperse_phase_solution.spheres_model_part,
                                                                 self.disperse_phase_solution.cluster_model_part,
                                                                 self.disperse_phase_solution.rigid_face_model_part,
                                                                 self.mixed_model_part)

        self.SetDragOutput()

        self.SetPointGraphPrinter()

        self.TransferGravityFromDisperseToFluid()

        # coarse-graining: applying changes to the physical properties of the model to adjust for
        # the similarity transformation if required (fluid effects only).
        SDP.ApplySimilarityTransformations(self.fluid_model_part,
                                           self.pp.CFD_DEM["similarity_transformation_type"].GetInt(),
                                           self.pp.CFD_DEM["model_over_real_diameter_factor"].GetDouble())

        self.SetPostUtils()

        # creating an IOTools object to perform other printing tasks
        self.io_tools = SDP.IOTools(self.pp)

        # creating a projection module for the fluid-DEM coupling
        self.h_min = 0.01
        n_balls = 1
        fluid_volume = 10
        self.pp.CFD_DEM.n_particles_in_depth = int(math.sqrt(n_balls / fluid_volume)) # only relevant in 2D problems
        # creating a physical calculations module to analyse the DEM model_part
        dem_physics_calculator = SphericElementGlobalPhysicsCalculator(self.disperse_phase_solution.spheres_model_part)

        if self.pp.CFD_DEM["coupling_level_type"].GetInt():

            if self.pp.CFD_DEM["meso_scale_length"].GetDouble() <= 0.0 and self.disperse_phase_solution.spheres_model_part.NumberOfElements(0) > 0:
                biggest_size = 2 * dem_physics_calculator.CalculateMaxNodalVariable(self.disperse_phase_solution.spheres_model_part, RADIUS)
                self.pp.CFD_DEM.meso_scale_length  = 20 * biggest_size

            elif self.disperse_phase_solution.spheres_model_part.NumberOfElements(0) == 0:
                self.pp.CFD_DEM.meso_scale_length  = 1.0

            self.projection_module = CFD_DEM_coupling.ProjectionModule(self.fluid_model_part,
                                                                       self.disperse_phase_solution.spheres_model_part,
                                                                       self.disperse_phase_solution.rigid_face_model_part,
                                                                       self.pp,
                                                                       flow_field = self.GetFieldUtility())
            self.projection_module.UpdateDatabase(self.h_min)

        # creating a custom functions calculator for the implementation of additional custom functions
        self.custom_functions_tool = SDP.FunctionsCalculator(self.pp)

        # creating a stationarity assessment tool
        self.stationarity_tool = SDP.StationarityAssessmentTool(self.pp.CFD_DEM["max_pressure_variation_rate_tol"].GetDouble() , self.custom_functions_tool)

        # creating a debug tool
        self.dem_volume_tool = self.GetVolumeDebugTool()

        #self.SetEmbeddedTools()

        Say('Initialization Complete\n')

        self.step           = 0
        self.time           = self.pp.Start_time
        self.Dt             = self.pp.Dt
        self.final_time     = self.pp.CFD_DEM["FinalTime"].GetDouble()
        self.report.Prepare(self.timer, self.pp.CFD_DEM["ControlTime"].GetDouble())

        #first_print = True; index_5 = 1; index_10 = 1; index_50 = 1; control = 0.0

        if self.pp.CFD_DEM["ModelDataInfo"].GetBool():
            os.chdir(data_and_results)
            if self.pp.CFD_DEM.ContactMeshOption == "ON":
                (coordination_number) = self.procedures.ModelData(self.disperse_phase_solution.spheres_model_part, self.solver) # Calculates the mean number of neighbours the mean radius, etc..
                Say('Coordination Number: ' + str(coordination_number) + '\n')
                os.chdir(self.main_path)
            else:
                Say('Activate Contact Mesh for ModelData information\n')

        if self.pp.CFD_DEM["flow_in_porous_medium_option"].GetBool():
            fluid_frac_util = SDP.FluidFractionFieldUtility(self.fluid_model_part, self.pp.CFD_DEM.min_fluid_fraction )

            for field in self.pp.fluid_fraction_fields:
                fluid_frac_util.AppendLinearField(field)

            fluid_frac_util.AddFluidFractionField()

        if self.pp.CFD_DEM["flow_in_porous_DEM_medium_option"].GetBool():
            SDP.FixModelPart(self.disperse_phase_solution.spheres_model_part)

        # choosing the directory in which we want to work (print to)

        os.chdir(self.post_path)



        ######################################################################################################################################

        #                      I N I T I A L I Z I N G    T I M E    L O O P     ...   ( M I X E D    F L U I D / D E M    B L O C K )

        ######################################################################################################################################

        self.DEM_step     = 0      # necessary to get a good random insertion of particles   # relevant to the stationarity assessment tool
        self.time_dem     = 0.0
        self.Dt_DEM       = self.disperse_phase_solution.spheres_model_part.ProcessInfo.GetValue(DELTA_TIME)
        self.disperse_phase_solution.rigid_face_model_part.ProcessInfo[DELTA_TIME] = self.Dt_DEM
        self.disperse_phase_solution.cluster_model_part.ProcessInfo[DELTA_TIME] = self.Dt_DEM
        self.stationarity = False

        # setting up loop counters: Counter(steps_per_tick_step, initial_step, active_or_inactive_boolean, dead_or_not)
        self.fluid_solve_counter          = self.GetFluidSolveCounter()
        #self.embedded_counter             = self.GetEmbeddedCounter()
        self.DEM_to_fluid_counter         = self.GetBackwardCouplingCounter()
        self.derivative_recovery_counter  = self.GetRecoveryCounter()
        self.stationarity_counter         = self.GetStationarityCounter()
        self.print_counter_updated_DEM    = self.GetPrintCounterUpdatedDEM()
        self.print_counter_updated_fluid  = self.GetPrintCounterUpdatedFluid()
        self.debug_info_counter           = self.GetDebugInfo()
        self.particles_results_counter    = self.GetParticlesResultsCounter()
        self.quadrature_counter           = self.GetHistoryForceQuadratureCounter()
        self.mat_deriv_averager           = SDP.Averager(1, 3)
        self.laplacian_averager           = SDP.Averager(1, 3)

        self.report.total_steps_expected = int(self.pp.CFD_DEM["FinalTime"].GetDouble() / self.Dt_DEM)

        Say(self.report.BeginReport(self.timer))

        # creating a Post Utils object that executes several post-related tasks
        self.post_utils_DEM = DEM_procedures.PostUtils(self.pp.CFD_DEM, self.disperse_phase_solution.spheres_model_part)

        SDP.InitializeVariablesWithNonZeroValues(self.fluid_model_part, self.disperse_phase_solution.spheres_model_part, self.pp) # otherwise variables are set to 0 by default

        self.SetUpResultsDatabase()

        # ANALYTICS BEGIN
        self.pp.CFD_DEM.perform_analytics_option = False

        if self.pp.CFD_DEM.perform_analytics_option:
            import analytics
            variables_to_measure = [PRESSURE]
            steps_between_measurements = 100
            gauge = analytics.Gauge(self.fluid_model_part, self.Dt, self.final_time, variables_to_measure, steps_between_measurements)
            point_coors = [0.0, 0.0, 0.01]
            target_node = SDP.FindClosestNode(self.fluid_model_part, point_coors)
            target_id = target_node.Id
            Say(target_node.X, target_node.Y, target_node.Z)
            Say(target_id)
            def condition(node):
                return node.Id == target_id

            gauge.ConstructArrayOfNodes(condition)
            Say(gauge.variables)
            #print_analytics_counter = SDP.Counter( 5 * steps_between_measurements, 1, 1) # MA: not used anywhere?
        # ANALYTICS END

        import derivative_recovery.derivative_recovery_strategy as derivative_recoverer

        self.recovery = derivative_recoverer.DerivativeRecoveryStrategy(self.pp, self.fluid_model_part, self.custom_functions_tool)

        self.FillHistoryForcePrecalculatedVectors()

        self.PerformZeroStepInitializations()

        self.post_utils.Writeresults(self.time)

    def AddExtraProcessInfoVariablesToFluid(self):
        vars_man.AddExtraProcessInfoVariablesToFluidModelPart(self.pp, self.fluid_model_part)

    def FluidInitialize(self):
        self.fluid_model_part = self.fluid_solution.fluid_model_part
        self.fluid_solution.vars_man = vars_man

        self.fluid_solution.SetFluidSolverModule()
        self.fluid_solution.AddFluidVariables()
        self.AddExtraProcessInfoVariablesToFluid()
        self.ReadFluidModelParts()
        self.fluid_solution.SetFluidBufferSizeAndAddDofs()
        SDP.AddExtraDofs(self.pp, self.fluid_model_part, self.disperse_phase_solution.spheres_model_part, self.disperse_phase_solution.cluster_model_part, self.disperse_phase_solution.DEM_inlet_model_part)
        self.fluid_solution.SetFluidSolver()
        self.fluid_solution.fluid_solver.Initialize()
        self.fluid_solution.ActivateTurbulenceModel()

    def ReadFluidModelParts(self):
        self.fluid_solution.ReadFluidModelPart()

    def DispersePhaseInitialize(self):
        self.spheres_model_part = self.disperse_phase_solution.spheres_model_part
        vars_man.AddNodalVariables(self.disperse_phase_solution.spheres_model_part, self.pp.dem_vars)
        vars_man.AddNodalVariables(self.disperse_phase_solution.rigid_face_model_part, self.pp.rigid_faces_vars)
        vars_man.AddNodalVariables(self.disperse_phase_solution.DEM_inlet_model_part, self.pp.inlet_vars)
        vars_man.AddExtraProcessInfoVariablesToDispersePhaseModelPart(self.pp, self.disperse_phase_solution.spheres_model_part)

        self.disperse_phase_solution.Initialize()

    def SetPostUtils(self):
          # creating a Post Utils object that executes several post-related tasks
        self.post_utils = SDP.PostUtils(self.swimming_DEM_gid_io,
                                        self.pp,
                                        self.fluid_model_part,
                                        self.disperse_phase_solution.spheres_model_part,
                                        self.disperse_phase_solution.cluster_model_part,
                                        self.disperse_phase_solution.rigid_face_model_part,
                                        self.mixed_model_part)

    def SetEmbeddedTools(self):
    # creating a distance calculation process for the embedded technology
        # (used to calculate elemental distances defining the structure embedded in the fluid mesh)
        if self.pp.CFD_DEM["embedded_option"].GetBool():
            self.calculate_distance_process = CalculateSignedDistanceTo3DSkinProcess(self.disperse_phase_solution.rigid_face_model_part, self.fluid_model_part)
            self.calculate_distance_process.Execute()

    def TheSimulationMustGoOn(self):
        return self.time <= self.final_time

    def RunMainTemporalLoop(self):

        while self.TheSimulationMustGoOn():

            self.time = self.time + self.Dt
            self.step += 1
            self.CloneTimeStep()
            self.TellTime(self.time)

            if self.pp.CFD_DEM["coupling_scheme_type"].GetString() == "UpdatedDEM":
                time_final_DEM_substepping = self.time + self.Dt

            else:
                time_final_DEM_substepping = self.time

            #self.PerformEmbeddedOperations() #TODO: it's crashing

            self.UpdateALEMeshMovement(self.time)

            # solving the fluid part
            if self.step >= self.GetFirstStepForFluidComputation():
                self.FluidSolve(self.time, solve_system = self.fluid_solve_counter.Tick() and not self.stationarity)

            # assessing stationarity

                if self.stationarity_counter.Tick():
                    Say("Assessing Stationarity...\n")
                    self.stationarity = self.stationarity_tool.Assess(self.fluid_model_part)
                    self.stationarity_counter.Deactivate(self.stationarity)

            # printing if required

            if self.particles_results_counter.Tick():
                self.io_tools.PrintParticlesResults(self.pp.variables_to_print_in_file, self.time, self.disperse_phase_solution.spheres_model_part)
                self.PrintDrag(self.drag_list, self.drag_file_output_list, self.fluid_model_part, self.time)

            if self.print_counter_updated_DEM.Tick():

                if coupling_level_type:
                    self.projection_module.ComputePostProcessResults(self.disperse_phase_solution.spheres_model_part.ProcessInfo)

                self.post_utils.Writeresults(self.time_dem)

            # solving the DEM part
            interaction_start_time = self.pp.CFD_DEM["interaction_start_time"].GetDouble()

            self.derivative_recovery_counter.Activate(self.time > interaction_start_time)

            if self.derivative_recovery_counter.Tick():
                self.recovery.Recover()

            Say('Solving DEM... (', self.disperse_phase_solution.spheres_model_part.NumberOfElements(0), 'elements )')
            first_dem_iter = True

            coupling_level_type = self.pp.CFD_DEM["coupling_level_type"].GetInt()
            project_at_every_substep_option = self.pp.CFD_DEM["project_at_every_substep_option"].GetBool()
            coupling_scheme_type = self.pp.CFD_DEM["coupling_scheme_type"].GetString()
            integration_scheme = self.pp.CFD_DEM["TranslationalIntegrationScheme"].GetString()
            basset_force_type = self.pp.CFD_DEM["basset_force_type"].GetInt()
            dem_inlet_option = self.pp.CFD_DEM["dem_inlet_option"].GetBool()

            for self.time_dem in self.yield_DEM_time(self.time_dem, time_final_DEM_substepping, self.Dt_DEM):
                self.DEM_step += 1   # this variable is necessary to get a good random insertion of particles
                self.disperse_phase_solution.spheres_model_part.ProcessInfo[TIME_STEPS]    = self.DEM_step
                self.disperse_phase_solution.rigid_face_model_part.ProcessInfo[TIME_STEPS] = self.DEM_step
                self.disperse_phase_solution.cluster_model_part.ProcessInfo[TIME_STEPS]    = self.DEM_step

                self.PerformInitialDEMStepOperations(self.time_dem)

                if self.time >= interaction_start_time and coupling_level_type and (project_at_every_substep_option or first_dem_iter):

                    if coupling_scheme_type == "UpdatedDEM":
                        self.ApplyForwardCoupling()

                    else:
                        self.ApplyForwardCoupling(alpha = 1.0 - (time_final_DEM_substepping - self.time_dem) / self.Dt)

                        if self.quadrature_counter.Tick():
                            self.AppendValuesForTheHistoryForce()

                        if integration_scheme in {'Hybrid_Bashforth', 'TerminalVelocityScheme'}:
                            # Advance in space only
                            self.DEMSolve(self.time_dem)
                            self.ApplyForwardCouplingOfVelocityToSlipVelocityOnly(self.time_dem)

                # performing the time integration of the DEM part

                self.disperse_phase_solution.spheres_model_part.ProcessInfo[TIME]    = self.time_dem
                self.disperse_phase_solution.rigid_face_model_part.ProcessInfo[TIME] = self.time_dem
                self.disperse_phase_solution.cluster_model_part.ProcessInfo[TIME]    = self.time_dem

                if self.do_solve_dem:
                    self.DEMSolve(self.time_dem)

                self.disperse_phase_solution.DEMFEMProcedures.MoveAllMeshes(self.all_model_parts, self.time_dem, self.Dt_DEM)

                #### TIME CONTROL ##################################

                # adding DEM elements by the inlet:
                if dem_inlet_option:
                    self.disperse_phase_solution.DEM_inlet.CreateElementsFromInletMesh(self.disperse_phase_solution.spheres_model_part, self.disperse_phase_solution.cluster_model_part, self.disperse_phase_solution.creator_destructor)  # After solving, to make sure that neighbours are already set.

                if self.print_counter_updated_fluid.Tick():

                    if coupling_level_type:
                        self.projection_module.ComputePostProcessResults(self.disperse_phase_solution.spheres_model_part.ProcessInfo)

                    self.post_utils.Writeresults(self.time_dem)

                first_dem_iter = False

                # applying DEM-to-fluid coupling

                if self.DEM_to_fluid_counter.Tick() and self.time >= interaction_start_time:
                    self.projection_module.ProjectFromParticles()

            #### PRINTING GRAPHS ####
            os.chdir(self.graphs_path)
            # measuring mean velocities in a certain control volume (the 'velocity trap')
            if self.pp.CFD_DEM["VelocityTrapOption"].GetBool():
                self.post_utils_DEM.ComputeMeanVelocitiesinTrap("Average_Velocity.txt", self.time)

            os.chdir(self.post_path)

            # coupling checks (debugging)
            if self.debug_info_counter.Tick():
                self.dem_volume_tool.UpdateDataAndPrint(self.pp.CFD_DEM["fluid_domain_volume"].GetDouble())

            # printing if required

            if self.particles_results_counter.Tick():
                self.io_tools.PrintParticlesResults(self.pp.variables_to_print_in_file, self.time, self.disperse_phase_solution.spheres_model_part)
                #self.graph_printer.PrintGraphs(self.time) #MA: commented out because the constructor was already commented out
                self.PrintDrag(self.drag_list, self.drag_file_output_list, self.fluid_model_part, self.time)



    def GetFirstStepForFluidComputation(self):
        return 3;

    def CloneTimeStep(self):
        self.fluid_model_part.CloneTimeStep(self.time)

    def DEMSolve(self, time = 'None'): # time is passed in case it is needed
        self.disperse_phase_solution.solver.Solve()

    def UpdateALEMeshMovement(self, time):
        pass

    def FluidSolve(self, time = 'None', solve_system = True):
        Say('Solving Fluid... (', self.fluid_model_part.NumberOfElements(0), 'elements )\n')

        if solve_system:
            self.fluid_solution.fluid_solver.Solve()
        else:
            Say("Skipping solving system...\n")

    def PerformZeroStepInitializations(self):
        pass

    def PerformInitialDEMStepOperations(self, time = None):
        pass

    def PerformEmbeddedOperations(self):
        # calculating elemental distances defining the structure embedded in the fluid mesh
        if self.pp.CFD_DEM["embedded_option"].GetBool():
            self.calculate_distance_process.Execute()

        if self.embedded_counter.Tick():
            embedded.ApplyEmbeddedBCsToFluid(self.fluid_model_part)
            embedded.ApplyEmbeddedBCsToBalls(self.disperse_phase_solution.spheres_model_part, self.pp.CFD_DEM)

    def SetInlet(self):
        if self.pp.CFD_DEM["dem_inlet_option"].GetBool():
            #Constructing the inlet and initializing it (must be done AFTER the self.disperse_phase_solution.spheres_model_part Initialize)
            # Note that right now only inlets of a single type are possible. This should be generalized.
            if self.pp.type_of_inlet == 'VelocityImposed':
                self.DEM_inlet = DEM_Inlet(self.DEM_inlet_model_part)
            elif self.pp.type_of_inlet == 'ForceImposed':
                self.DEM_inlet = DEM_Force_Based_Inlet(self.DEM_inlet_model_part, self.pp.force)

            self.DEM_inlet.InitializeDEM_Inlet(self.disperse_phase_solution.spheres_model_part, self.creator_destructor)

    def SetAnalyticFaceWatcher(self):
        from analytic_tools import analytic_data_procedures
        self.watcher = AnalyticFaceWatcher()
        self.watcher_analyser = analytic_data_procedures.FaceWatcherAnalyzer(analytic_face_watcher = self.watcher, path = self.main_path)

    def SetAnalyticParticleWatcher(self):
        from analytic_tools import analytic_data_procedures
        self.particle_watcher = AnalyticParticleWatcher()
        self.particle_watcher_analyser = analytic_data_procedures.ParticleWatcherAnalyzer(analytic_particle_watcher = self.particle_watcher, path = self.main_path)

    def SetInletWatcher(self):
        self.watcher_analyser.SetInlet(self.DEM_inlet)

    def TellTime(self, time):
        Say('\nTIME = ', time)
        Say('ELAPSED TIME = ', self.timer.time() - self.simulation_start_time, '\n')

    def TellFinalSummary(self, time, step, DEM_step):
        Say('*************************************************************')
        Say('CALCULATIONS FINISHED. THE SIMULATION ENDED SUCCESSFULLY.')
        simulation_elapsed_time = self.timer.time() - self.simulation_start_time
        Say('Elapsed time: ' + '%.5f'%(simulation_elapsed_time) + ' s ')
        Say('per fluid time step: ' + '%.5f'%(simulation_elapsed_time / step) + ' s ')
        Say('per DEM time step: ' + '%.5f'%(simulation_elapsed_time / DEM_step) + ' s ')
        Say('*************************************************************\n')

    def GetFluidSolveCounter(self):
        return SDP.Counter(is_dead = (self.pp.CFD_DEM["drag_force_type"].GetInt() == 9))

    def GetEmbeddedCounter(self):
        return SDP.Counter(1, 3, self.pp.CFD_DEM["embedded_option"].GetBool())  # MA: because I think DISTANCE,1 (from previous time step) is not calculated correctly for step=1

    def GetBackwardCouplingCounter(self):
        return SDP.Counter(1, 1, self.pp.CFD_DEM["coupling_level_type"].GetInt() > 1)

    def GetRecoveryCounter(self):
        return SDP.Counter(1, 1, self.pp.CFD_DEM["coupling_level_type"].GetInt() or self.pp.CFD_DEM["print_PRESSURE_GRADIENT_option"].GetBool())

    def GetStationarityCounter(self):
        return SDP.Counter(steps_in_cycle = self.pp.CFD_DEM["time_steps_per_stationarity_step"].GetInt(),
                           beginning_step = 1,
                           is_active = self.pp.CFD_DEM["stationary_problem_option"].GetBool())

    def GetPrintCounterUpdatedDEM(self):
        counter = SDP.Counter(steps_in_cycle = int(self.output_time / self.Dt_DEM + 0.5),
                                     beginning_step = int(self.Dt / self.Dt_DEM))

        if 'UpdatedDEM' != self.pp.CFD_DEM["coupling_scheme_type"].GetString():
            counter.Kill()
        return counter

    def GetPrintCounterUpdatedFluid(self):
        counter = SDP.Counter(steps_in_cycle = int(self.output_time / self.Dt_DEM + 0.5),
                           beginning_step = int(self.Dt / self.Dt_DEM))

        if 'UpdatedFluid' != self.pp.CFD_DEM["coupling_scheme_type"].GetString():
            counter.Kill()
        return counter

    def GetDebugInfo(self):
        return SDP.Counter(self.pp.CFD_DEM["debug_tool_cycle"].GetInt(), 1, self.pp.CFD_DEM["print_debug_info_option"].GetBool())

    def GetParticlesResultsCounter(self):
        return SDP.Counter(self.pp.CFD_DEM["print_particles_results_cycle"].GetInt(), 1, self.pp.CFD_DEM["print_particles_results_option"].GetBool())

    def GetHistoryForceQuadratureCounter(self):
        return SDP.Counter(self.pp.CFD_DEM["time_steps_per_quadrature_step"].GetInt(), 1, self.pp.CFD_DEM["basset_force_type"].GetInt())

    def GetVolumeDebugTool(self):
        return SDP.ProjectionDebugUtils(self.pp.CFD_DEM["fluid_domain_volume"].GetDouble(), self.fluid_model_part, self.disperse_phase_solution.spheres_model_part, self.custom_functions_tool)

    def GetRunCode(self):
        return SDP.CreateRunCode(self.pp)

    def FillHistoryForcePrecalculatedVectors(self):
        # Warning: this estimation is based on a constant time step for DEM. This is usually the case, but could not be so. A more robust implementation is needed!
        N_steps = int(self.pp.CFD_DEM["FinalTime"].GetDouble() / self.pp.CFD_DEM["MaxTimeStep"].GetDouble()) + 20
        spheres_model_part = self.all_model_parts.Get('SpheresPart')
        if self.pp.CFD_DEM["basset_force_type"].GetInt() > 0:
            self.basset_force_tool.FillDaitcheVectors(N_steps, self.pp.CFD_DEM["quadrature_order"].GetInt(), self.pp.CFD_DEM["time_steps_per_quadrature_step"].GetInt())
        if self.pp.CFD_DEM["basset_force_type"].GetInt() >= 3 or self.pp.CFD_DEM["basset_force_type"].GetInt() == 1:
            self.basset_force_tool.FillHinsbergVectors(spheres_model_part, self.pp.CFD_DEM["number_of_exponentials"].GetInt(), self.pp.CFD_DEM["number_of_quadrature_steps_in_window"].GetInt())

    def AppendValuesForTheHistoryForce(self):
        spheres_model_part = self.all_model_parts.Get('SpheresPart')

        if self.pp.CFD_DEM["basset_force_type"].GetInt() == 1 or self.pp.CFD_DEM["basset_force_type"].GetInt() >= 3:
            self.basset_force_tool.AppendIntegrandsWindow(spheres_model_part)
        elif self.pp.CFD_DEM["basset_force_type"].GetInt() == 2:
            self.basset_force_tool.AppendIntegrands(spheres_model_part)

    def GetBassetForceTools(self):
        self.basset_force_tool = SDP.BassetForceTools()

    def GetFieldUtility(self):
        return None

    def GetResultsCreator(self):
        return None

    def ApplyForwardCoupling(self, alpha = 'None'):
        self.projection_module.ApplyForwardCoupling(alpha)

    def ApplyForwardCouplingOfVelocityToSlipVelocityOnly(self, time = None):
        self.projection_module.ApplyForwardCouplingOfVelocityToSlipVelocityOnly()

    def PerformFinalOperations(self, time = None):
        os.chdir(self.main_path)
        del self.post_utils
        self.ModifyResultsFolderName(time)

    def ModifyResultsFolderName(self, time):
        pass

    def Finalize(self):

        self.swimming_DEM_gid_io.finalize_results()

        self.PerformFinalOperations(self.time_dem)

        self.FinalizeDragOutput()

        self.TellFinalSummary(self.step, self.time, self.DEM_step)

    def FinalizeDragOutput(self):
        for i in self.drag_file_output_list:
            i.close()

    def SetCutsOutput(self):
        if not self.pp.VolumeOutput:
            cut_list = define_output.DefineCutPlanes()
            self.swimming_DEM_gid_io.define_cuts(self.fluid_model_part, cut_list)

    def SetDragOutput(self):
        # define the drag computation list
        self.drag_list = define_output.DefineDragList()
        self.drag_file_output_list = []
        for it in self.drag_list:
            f = open(it[1], 'w')
            self.drag_file_output_list.append(f)
            tmp = "#Drag for group " + it[1] + "\n"
            f.write(tmp)
            tmp = "time RX RY RZ"
            f.write(tmp)
            f.flush()

        Say(self.drag_file_output_list)

    def SetPointGraphPrinter(self):
        pass
         # preparing output of point graphs
        #import point_graph_printer

        #output_nodes_list = define_output.DefineOutputPoints()
        #self.graph_printer = point_graph_printer.PrintGraphPrinter(
            #output_nodes_list,
            #fluid_model_part,
            #variables_dictionary,
            #domain_size)

    def TransferGravityFromDisperseToFluid(self):
        # setting fluid's body force to the same as DEM's
        if self.pp.CFD_DEM["body_force_on_fluid_option"].GetBool():

            for node in self.fluid_model_part.Nodes:
                node.SetSolutionStepValue(BODY_FORCE_X, 0, self.pp.CFD_DEM["GravityX"].GetDouble())
                node.SetSolutionStepValue(BODY_FORCE_Y, 0, self.pp.CFD_DEM["GravityY"].GetDouble())
                node.SetSolutionStepValue(BODY_FORCE_Z, 0, self.pp.CFD_DEM["GravityZ"].GetDouble())

    def yield_DEM_time(self, current_time, current_time_plus_increment, delta_time):
            current_time += delta_time

            tolerance = 0.0001
            while current_time < (current_time_plus_increment - tolerance * delta_time):
                yield current_time
                current_time += delta_time

            current_time = current_time_plus_increment
            yield current_time

    def PrintDrag(self, drag_list, drag_file_output_list, fluid_model_part, time):
        i = 0
        for it in drag_list:
            nodes = self.fluid_model_part.GetNodes(it[0])
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
            # print self.drag_file_output_list[i]
            # print output
            self.drag_file_output_list[i].write(output)
            self.drag_file_output_list[i].flush()
            i = i + 1

    def GetReturnValue(self):
        return 0.0

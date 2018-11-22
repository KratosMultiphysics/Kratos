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

from analysis_stage import AnalysisStage

import CFD_DEM_coupling
import swimming_DEM_procedures as SDP
import swimming_DEM_gid_output
import embedded
import variables_management as vars_man

def Say(*args):
    Logger.PrintInfo("DEM-FLUID", *args)
    Logger.Flush()


try:
    import define_output  # MA: some GUI write this file, some others not!
except ImportError:
    pass

# Import MPI modules if needed. This way to do this is only valid when using OpenMPI.
# For other implementations of MPI it will not work.
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

class SDEMLogger(object):
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

class SwimmingDEMAnalysis(AnalysisStage):
    def __enter__ (self):
        # sys.stdout = SDEMLogger()
        return self

    def __exit__(self, exception_type, exception_value, traceback):
        pass

    def __init__(self, model, parameters = Parameters("{}")):
        sys.stdout = SDEMLogger()
        self.StartTimer()
        self.model = model
        self.main_path = os.getcwd()

        self.SetFluidAlgorithm()
        self.fluid_solution.coupling_analysis = weakref.proxy(self)

        self.pp = self.fluid_solution.pp

        self.SetCouplingParameters(parameters)

        self.SetDispersePhaseAlgorithm()
        self.disperse_phase_solution.coupling_analysis = weakref.proxy(self)

        self.procedures = weakref.proxy(self.disperse_phase_solution.procedures)
        self.report = DEM_procedures.Report()

        # creating a basset_force tool to perform the operations associated
        # with the calculation of this force along the path of each particle
        self.GetBassetForceTools()
        self.disperse_phase_solution.SetAnalyticFaceWatcher()

        # defining member variables for the model_parts (for convenience)
        self.fluid_model_part = self.fluid_solution.fluid_model_part
        self.spheres_model_part = self.disperse_phase_solution.spheres_model_part
        self.cluster_model_part = self.disperse_phase_solution.cluster_model_part
        self.rigid_face_model_part = self.disperse_phase_solution.rigid_face_model_part
        self.DEM_inlet_model_part = self.disperse_phase_solution.DEM_inlet_model_part
        super(SwimmingDEMAnalysis, self).__init__(model, self.pp.fluid_parameters)

    def SetFluidAlgorithm(self):
        import DEM_coupled_fluid_dynamics_analysis
        self.fluid_solution = DEM_coupled_fluid_dynamics_analysis.DEMCoupledFluidDynamicsAnalysis(self.model)
        self.fluid_solution.main_path = self.main_path

    def SetDispersePhaseAlgorithm(self):
        import fluid_coupled_DEM_analysis as DEM_analysis
        self.disperse_phase_solution = DEM_analysis.Solution(self.model, self.pp)

    def ReadDispersePhaseAndCouplingParameters(self):

        with open(self.main_path + '/ProjectParametersDEM.json', 'r') as parameters_file:
            self.pp.CFD_DEM = Parameters(parameters_file.read())

        import dem_default_input_parameters
        dem_defaults = dem_default_input_parameters.GetDefaultInputParameters()

        import swimming_dem_default_input_parameters
        only_swimming_defaults = swimming_dem_default_input_parameters.GetDefaultInputParameters()

        for key in only_swimming_defaults.keys():
            dem_defaults.AddValue(key, only_swimming_defaults[key])

        self.pp.CFD_DEM.ValidateAndAssignDefaults(dem_defaults)

    def SetCouplingParameters(self, parameters):

        # First, read the parameters generated from the interface
        self.ReadDispersePhaseAndCouplingParameters()

        # Second, set the default 'beta' parameters (candidates to be moved to the interface)
        self.SetBetaParameters()

        # Third, set the parameters fed to the particular case that you are running
        self.SetCustomBetaParameters(parameters)

        # Finally adjust some of the parameters for consistency
        #   This function should be reduced to a minimum since,
        #   in principle, there should be no need to change the parameters
        self.SetDerivedParameters()

    def SetDerivedParameters(self):
        self.SetDoSolveDEMVariable()

    def SetAllModelParts(self):
        self.all_model_parts = weakref.proxy(self.disperse_phase_solution.all_model_parts)

        # defining a fluid model
        self.all_model_parts.Add(self.fluid_model_part)

        # defining a model part for the mixed part
        self.all_model_parts.Add(self.model.CreateModelPart("MixedPart"))

        self.mixed_model_part = self.all_model_parts.Get('MixedPart')

    def StartTimer(self):
        self.timer = timer
        self.simulation_start_time = timer.time()

    def SetBetaParameters(self):
        # These are input parameters that have not yet been transferred to the interface
        # import the configuration data as read from the GiD
        Add = self.pp.CFD_DEM.AddEmptyValue
        Add("fluid_already_calculated").SetBool(False)
        Add("do_solve_dem").SetBool(True)
        Add("recovery_echo_level").SetInt(1)
        Add("gradient_calculation_type").SetInt(1)
        Add("pressure_grad_recovery_type").SetInt(0)
        Add("fluid_fraction_grad_type").SetInt(0)
        Add("store_full_gradient_option").SetBool(False)
        Add("store_fluid_pressure_option").SetBool(False)
        Add("laplacian_calculation_type").SetInt(0)
        Add("faxen_terms_type").SetInt(0)
        Add("material_acceleration_calculation_type").SetInt(1)
        Add("faxen_force_type").SetInt(0)
        Add("vorticity_calculation_type").SetInt(5)
        Add("print_FLUID_VEL_PROJECTED_RATE_option").SetBool(False)
        Add("print_MATERIAL_FLUID_ACCEL_PROJECTED_option").SetBool(True)
        Add("basset_force_type").SetInt(0)
        Add("print_BASSET_FORCE_option").SetBool(True)
        Add("basset_force_integration_type").SetInt(2)
        Add("n_init_basset_steps").SetInt(0)
        Add("time_steps_per_quadrature_step").SetInt(1)
        time_steps_per_quadrature_step = self.pp.CFD_DEM["time_steps_per_quadrature_step"].GetInt()
        DEM_dt = self.pp.CFD_DEM["MaxTimeStep"].GetDouble()
        Add("delta_time_quadrature").SetDouble(time_steps_per_quadrature_step * DEM_dt)
        Add("quadrature_order").SetInt(2)
        Add("time_window").SetDouble(0.8)
        Add("number_of_exponentials").SetInt(2)
        time_window = self.pp.CFD_DEM["time_window"].GetDouble()
        quadrature_dt = self.pp.CFD_DEM["delta_time_quadrature"].GetDouble()
        Add("number_of_quadrature_steps_in_window").SetInt(int(time_window / quadrature_dt))
        Add("do_impose_flow_from_field_option").SetBool(False)
        Add("print_MATERIAL_ACCELERATION_option").SetBool(True)
        Add("print_FLUID_ACCEL_FOLLOWING_PARTICLE_PROJECTED_option").SetBool(False)
        Add("print_VORTICITY_option").SetBool(True)
        Add("print_MATERIAL_ACCELERATION_option").SetBool(True)
        Add("print_VISCOSITY_option").SetBool(False)
        Add("print_VELOCITY_GRADIENT_option").SetBool(True)
        Add("print_DISPERSE_FRACTION_option").SetBool(False)
        Add("print_FLUID_FRACTION_GRADIENT_option").SetBool(False)
        Add("print_FLUID_FRACTION_GRADIENT_PROJECTED_option").SetBool(False)
        Add("print_VECTORIAL_ERROR_option").SetBool(False)
        Add("calculate_diffusivity_option").SetBool(False)
        Add("print_CONDUCTIVITY_option").SetBool(False)
        Add("filter_velocity_option").SetBool(False)
        Add("print_PARTICLE_VEL_option").SetBool(False)
        Add("apply_time_filter_to_fluid_fraction_option").SetBool(False)
        Add("full_particle_history_watcher").SetString("Empty")
        Add("prerun_fluid_file_name").SetString("")
        Add("frame_of_reference_type").SetInt(0)
        Add("angular_velocity_of_frame_X").SetDouble(0.0)
        Add("angular_velocity_of_frame_Y").SetDouble(0.0)
        Add("angular_velocity_of_frame_Z").SetDouble(0.0)
        Add("angular_velocity_of_frame_old_X").SetDouble(0.0)
        Add("angular_velocity_of_frame_old_Y").SetDouble(0.0)
        Add("angular_velocity_of_frame_old_Z").SetDouble(0.0)
        Add("acceleration_of_frame_origin_X").SetDouble(0.0)
        Add("acceleration_of_frame_origin_Y").SetDouble(0.0)
        Add("acceleration_of_frame_origin_Z").SetDouble(0.0)
        Add("angular_acceleration_of_frame_X").SetDouble(0.0)
        Add("angular_acceleration_of_frame_Y").SetDouble(0.0)
        Add("angular_acceleration_of_frame_Z").SetDouble(0.0)
        Add("ALE_option").SetBool(False)
        Add("frame_rotation_axis_initial_point").SetVector(Vector([0., 0., 0.]))
        Add("frame_rotation_axis_final_point").SetVector(Vector([0., 0., 1.]))
        Add("angular_velocity_magnitude").SetDouble(1.0)
        Add("print_distance_option").SetBool(False)
        Add("print_SLIP_VELOCITY_option").SetBool(False)

        # Making all time steps be an exact multiple of the smallest time step
        self.Dt_DEM = self.pp.CFD_DEM["MaxTimeStep"].GetDouble()
        self.Dt = self.fluid_solution._GetSolver().settings["time_stepping"]["time_step"].GetDouble()
        self.Dt = int(self.Dt / self.Dt_DEM) * self.Dt_DEM
        self.fluid_solution._GetSolver().settings["time_stepping"]["time_step"].SetDouble(self.Dt)
        self.output_time = self.pp.CFD_DEM["OutputTimeStep"].GetDouble()
        self.output_time = int(self.output_time / self.Dt_DEM) * self.Dt_DEM

        # vestigial variables
        Add("print_steps_per_plot_step").SetInt(1)
        Add("PostCationConcentration").SetBool(False)
        self.pp.viscosity_modification_type = 0.0
        self.domain_size = 3
        self.pp.domain_size = 3
        self.pp.type_of_inlet = 'VelocityImposed' # 'VelocityImposed' or 'ForceImposed'
        self.pp.force = Vector(3)
        self.pp.force[0] = 0
        self.pp.force[1] = 0
        self.pp.force[2] = 1e-10

        self.pp.fluid_fraction_fields = []
        field1 = SDP.FluidFractionFieldUtility.LinearField(0.0,
                                                           [0.0, 0.0, 0.0],
                                                           [-1.0, -1.0, 0.15],
                                                           [1.0, 1.0, 0.3]
                                                          )
        # write down the volume you know it has
        from math import pi
        Add("fluid_domain_volume").SetDouble(0.5 ** 2 * 2 * pi)

        self.pp.fluid_fraction_fields.append(field1)

        # Setting body_force_per_unit_mass_variable_name
        Add("body_force_per_unit_mass_variable_name").SetString('BODY_FORCE')

    def SetDoSolveDEMVariable(self):
        self.do_solve_dem = self.pp.CFD_DEM["do_solve_dem"].GetBool()

        if self.pp.CFD_DEM["flow_in_porous_DEM_medium_option"].GetBool():
            self.do_solve_dem = False

    def SetCustomBetaParameters(self, custom_parameters):
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

    def ReadDispersePhaseModelParts(self,
                                    starting_node_Id=0,
                                    starting_elem_Id=0,
                                    starting_cond_Id=0):
        creator_destructor = self.disperse_phase_solution.creator_destructor
        max_node_Id = creator_destructor.FindMaxNodeIdInModelPart(self.fluid_model_part)
        max_elem_Id = creator_destructor.FindMaxElementIdInModelPart(self.fluid_model_part)
        max_cond_Id = creator_destructor.FindMaxConditionIdInModelPart(self.fluid_model_part)
        self.disperse_phase_solution.BaseReadModelParts(max_node_Id, max_elem_Id, max_cond_Id)

    def Initialize(self):
        Say('Initializing Problem...\n')

        self.run_code = self.GetRunCode()

        # Moving to the recently created folder
        os.chdir(self.main_path)
        [self.post_path, data_and_results, self.graphs_path, MPI_results] = \
        self.procedures.CreateDirectories(str(self.main_path),
                                          str(self.pp.CFD_DEM["problem_name"].GetString()),
                                          self.run_code)
        SDP.CopyInputFilesIntoFolder(self.main_path, self.post_path)
        self.MPI_results = MPI_results
        #self.mixed_model_part = self.all_model_parts.Get('MixedPart')

        vars_man.ConstructListsOfVariables(self.pp)

        self.TransferBodyForceFromDisperseToFluid()

        self.FluidInitialize()
        self.DispersePhaseInitialize()

        self.SetAllModelParts()

        self.SetCutsOutput()
        gid_output_options = self.pp.fluid_parameters["output_processes"]["gid_output"][0]["Parameters"]
        result_file_configuration = gid_output_options["postprocess_parameters"]["result_file_configuration"]
        write_conditions_option = result_file_configuration["gidpost_flags"]["WriteConditionsFlag"].GetString() == "WriteConditionsFlag"
        deformed_mesh_option = result_file_configuration["gidpost_flags"]["WriteDeformedMeshFlag"].GetString() == "WriteDeformed"
        old_gid_output_post_options_dict = {'GiD_PostAscii':'Ascii','GiD_PostBinary':'Binary','GiD_PostAsciiZipped':'AsciiZipped'}
        old_gid_output_multiple_file_option_dict = {'SingleFile':'Single','MultipleFiles':'Multiples'}
        post_mode_key = result_file_configuration["gidpost_flags"]["GiDPostMode"].GetString()
        multiple_files_option_key = result_file_configuration["gidpost_flags"]["MultiFileFlag"].GetString()
        self.pp.GiDMultiFileFlag = old_gid_output_multiple_file_option_dict[multiple_files_option_key]
        self.swimming_DEM_gid_io = \
        swimming_DEM_gid_output.SwimmingDEMGiDOutput(
            file_name = self.pp.CFD_DEM["problem_name"].GetString(),
            vol_output = result_file_configuration["body_output"].GetBool(),
            post_mode = old_gid_output_post_options_dict[post_mode_key],
            multifile = old_gid_output_multiple_file_option_dict[multiple_files_option_key],
            deformed_mesh = deformed_mesh_option,
            write_conditions = write_conditions_option
            )

        self.swimming_DEM_gid_io.initialize_swimming_DEM_results(
            self.spheres_model_part,
            self.cluster_model_part,
            self.rigid_face_model_part,
            self.mixed_model_part
            )

        self.SetDragOutput()

        self.SetPointGraphPrinter()

        self.AssignKinematicViscosityFromDynamicViscosity()

        # coarse-graining: applying changes to the physical properties of the model to adjust for
        # the similarity transformation if required (fluid effects only).
        SDP.ApplySimilarityTransformations(
            self.fluid_model_part,
            self.pp.CFD_DEM["similarity_transformation_type"].GetInt(),
            self.pp.CFD_DEM["model_over_real_diameter_factor"].GetDouble()
            )

        self.SetPostUtils()

        # creating an IOTools object to perform other printing tasks
        self.io_tools = SDP.IOTools(self.pp)

        # creating a projection module for the fluid-DEM coupling
        self.h_min = 0.01
        n_balls = 1
        fluid_volume = 10
        # the variable n_particles_in_depth is only relevant in 2D problems
        self.pp.CFD_DEM.AddEmptyValue("n_particles_in_depth").SetInt(int(math.sqrt(n_balls / fluid_volume)))
        # creating a physical calculations module to analyse the DEM model_part

        dem_physics_calculator = SphericElementGlobalPhysicsCalculator(
            self.spheres_model_part
            )

        if self.pp.CFD_DEM["coupling_level_type"].GetInt():
            default_meso_scale_length_needed = (
                self.pp.CFD_DEM["meso_scale_length"].GetDouble() <= 0.0 and
                self.spheres_model_part.NumberOfElements(0) > 0
                )
            if default_meso_scale_length_needed:
                biggest_size = (
                    2 * dem_physics_calculator.CalculateMaxNodalVariable(
                        self.spheres_model_part,
                        RADIUS
                        )
                    )
                self.pp.CFD_DEM["meso_scale_length"].SetDouble(20 * biggest_size)

            elif self.spheres_model_part.NumberOfElements(0) == 0:
                self.pp.CFD_DEM["meso_scale_length"].SetDouble(1.0)

            self.projection_module = CFD_DEM_coupling.ProjectionModule(
                self.fluid_model_part,
                self.spheres_model_part,
                self.rigid_face_model_part,
                self.pp,
                flow_field=self.GetFieldUtility()
                )

            self.projection_module.UpdateDatabase(self.h_min)

        # creating a custom functions calculator for the implementation of
        # additional custom functions
        self.custom_functions_tool = SDP.FunctionsCalculator(self.pp)

        # creating a stationarity assessment tool
        self.stationarity_tool = SDP.StationarityAssessmentTool(
            self.pp.CFD_DEM["max_pressure_variation_rate_tol"].GetDouble(),
            self.custom_functions_tool
            )

        # creating a debug tool
        self.dem_volume_tool = self.GetVolumeDebugTool()

        #self.SetEmbeddedTools()

        Say('Initialization Complete\n')

        self.report.Prepare(self.timer, self.pp.CFD_DEM["ControlTime"].GetDouble())

        #first_print = True; index_5 = 1; index_10 = 1; index_50 = 1; control = 0.0

        if self.pp.CFD_DEM["ModelDataInfo"].GetBool():
            os.chdir(data_and_results)
            if self.pp.CFD_DEM.ContactMeshOption == "ON":
                coordination_number = self.procedures.ModelData(
                    self.spheres_model_part,
                    self.solver)
                Say('Coordination Number: ' + str(coordination_number) + '\n')
                os.chdir(self.main_path)
            else:
                Say('Activate Contact Mesh for ModelData information\n')

        if self.pp.CFD_DEM["flow_in_porous_medium_option"].GetBool():
            fluid_frac_util = SDP.FluidFractionFieldUtility(
                self.fluid_model_part, self.pp.CFD_DEM.min_fluid_fraction)

            for field in self.pp.fluid_fraction_fields:
                fluid_frac_util.AppendLinearField(field)

            fluid_frac_util.AddFluidFractionField()

        if self.pp.CFD_DEM["flow_in_porous_DEM_medium_option"].GetBool():
            SDP.FixModelPart(self.spheres_model_part)

        # choosing the directory in which we want to work (print to)

        os.chdir(self.post_path)

        ##################################################

        #    I N I T I A L I Z I N G    T I M E    L O O P

        ##################################################
        self.step       = 0
        self.time       = self.pp.fluid_parameters["problem_data"]["start_time"].GetDouble()
        self.Dt         = self.fluid_solution._GetSolver()._ComputeDeltaTime()
        self.final_time = self.pp.CFD_DEM["FinalTime"].GetDouble()
        self.DEM_step   = 0
        self.time_dem   = 0.0
        self.Dt_DEM     = self.spheres_model_part.ProcessInfo.GetValue(DELTA_TIME)
        self.rigid_face_model_part.ProcessInfo[DELTA_TIME] = self.Dt_DEM
        self.cluster_model_part.ProcessInfo[DELTA_TIME] = self.Dt_DEM
        self.stationarity = False

        # setting up loop counters:
        self.fluid_solve_counter          = self.GetFluidSolveCounter()
        self.DEM_to_fluid_counter         = self.GetBackwardCouplingCounter()
        self.derivative_recovery_counter  = self.GetRecoveryCounter()
        self.stationarity_counter         = self.GetStationarityCounter()
        self.print_counter_updated_DEM    = self.GetPrintCounterUpdatedDEM()
        self.print_counter_updated_fluid  = self.GetPrintCounterUpdatedFluid()
        self.debug_info_counter           = self.GetDebugInfo()
        self.particles_results_counter    = self.GetParticlesResultsCounter()
        self.quadrature_counter           = self.GetHistoryForceQuadratureCounter()
        #Phantom
        self.analytic_data_counter        = self.ProcessAnalyticDataCounter()
        self.mat_deriv_averager           = SDP.Averager(1, 3)
        self.laplacian_averager           = SDP.Averager(1, 3)

        self.report.total_steps_expected = int(self.final_time / self.Dt_DEM)

        Say(self.report.BeginReport(self.timer))

        # creating a Post Utils object that executes several post-related tasks
        self.post_utils_DEM = DEM_procedures.PostUtils(self.pp.CFD_DEM, self.spheres_model_part)

        SDP.InitializeVariablesWithNonZeroValues(
            self.fluid_model_part,
            self.spheres_model_part,
            self.pp
            ) # otherwise variables are set to 0 by default

        self.SetUpResultsDatabase()

        # ANALYTICS BEGIN
        self.pp.CFD_DEM.AddEmptyValue("perform_analytics_option").SetBool(False)

        if self.pp.CFD_DEM["perform_analytics_option"].GetBool():
            import analytics
            variables_to_measure = [PRESSURE]
            steps_between_measurements = 100
            gauge = analytics.Gauge(
                self.fluid_model_part,
                self.Dt,
                self.final_time,
                variables_to_measure,
                steps_between_measurements
                )
            point_coors = [0.0, 0.0, 0.01]
            target_node = SDP.FindClosestNode(self.fluid_model_part, point_coors)
            target_id = target_node.Id
            Say(target_node.X, target_node.Y, target_node.Z)
            Say(target_id)
            def condition(node):
                return node.Id == target_id

            gauge.ConstructArrayOfNodes(condition)
            Say(gauge.variables)
        # ANALYTICS END

        import derivative_recovery.derivative_recovery_strategy as derivative_recoverer

        self.recovery = derivative_recoverer.DerivativeRecoveryStrategy(
            self.pp,
            self.fluid_model_part,
            self.custom_functions_tool)

        self.FillHistoryForcePrecalculatedVectors()

        self.PerformZeroStepInitializations()

        self.post_utils.Writeresults(self.time)

    def AddExtraProcessInfoVariablesToFluid(self):
        vars_man.AddExtraProcessInfoVariablesToFluidModelPart(self.pp, self.fluid_model_part)

    def FluidInitialize(self):
        self.fluid_model_part = self.fluid_solution.fluid_model_part
        self.fluid_solution.vars_man = vars_man
        self.fluid_solution.Initialize()

        self.AddExtraProcessInfoVariablesToFluid()

        SDP.AddExtraDofs(
            self.pp,
            self.fluid_model_part,
            self.spheres_model_part,
            self.cluster_model_part,
            self.DEM_inlet_model_part)

    def DispersePhaseInitialize(self):
        vars_man.AddNodalVariables(self.spheres_model_part, self.pp.dem_vars)
        vars_man.AddNodalVariables(self.rigid_face_model_part, self.pp.rigid_faces_vars)
        vars_man.AddNodalVariables(self.DEM_inlet_model_part, self.pp.inlet_vars)
        vars_man.AddExtraProcessInfoVariablesToDispersePhaseModelPart(
            self.pp,
            self.disperse_phase_solution.spheres_model_part)

        self.disperse_phase_solution.Initialize()

    def SetPostUtils(self):
          # creating a Post Utils object that executes several post-related tasks
        self.post_utils = SDP.PostUtils(self.swimming_DEM_gid_io,
                                        self.pp,
                                        self.fluid_model_part,
                                        self.spheres_model_part,
                                        self.cluster_model_part,
                                        self.rigid_face_model_part,
                                        self.mixed_model_part)

    def SetEmbeddedTools(self):
    # creating a distance calculation process for the embedded technology
        # (used to calculate elemental distances defining the structure embedded in the fluid mesh)
        if self.pp.CFD_DEM["embedded_option"].GetBool():
            self.calculate_distance_process = CalculateSignedDistanceTo3DSkinProcess(
                self.rigid_face_model_part,
                self.fluid_model_part
                )
            self.calculate_distance_process.Execute()

    def TheSimulationMustGoOn(self):
        return self.time <= self.final_time

    def GetAnalyticFacesModelParts(self):
        analytic_face_submodelpart_number = 1
        analytic_face_submodelpart_name = self.rigid_face_model_part.GetSubModelPart(str(analytic_face_submodelpart_number))
        return analytic_face_submodelpart_name

    def MakeAnalyticsMeasurements(self):
        self.analytic_face_watcher.MakeMeasurements()
        self.analytic_particle_watcher.MakeMeasurements()

    def RunMainTemporalLoop(self):
        coupling_level_type = self.pp.CFD_DEM["coupling_level_type"].GetInt()
        project_at_every_substep_option = self.pp.CFD_DEM["project_at_every_substep_option"].GetBool()
        coupling_scheme_type = self.pp.CFD_DEM["coupling_scheme_type"].GetString()
        integration_scheme = self.pp.CFD_DEM["TranslationalIntegrationScheme"].GetString()
        dem_inlet_option = self.pp.CFD_DEM["dem_inlet_option"].GetBool()
        interaction_start_time = self.pp.CFD_DEM["interaction_start_time"].GetDouble()

        while self.TheSimulationMustGoOn():

            #self.time = self.time + self.Dt
            self.time = self.fluid_solution._GetSolver().AdvanceInTime(self.time)
            self.step += 1
            self.TellTime(self.time)

            if coupling_scheme_type == "UpdatedDEM":
                time_final_DEM_substepping = self.time + self.Dt

            else:
                time_final_DEM_substepping = self.time

            #self.PerformEmbeddedOperations() TO-DO: it's crashing

            self.UpdateALEMeshMovement(self.time)

            # solving the fluid part
            if self.step >= self.GetFirstStepForFluidComputation():
                self.FluidSolve(
                    self.time,
                    solve_system=self.fluid_solve_counter.Tick() and not self.stationarity
                    )

            # assessing stationarity
                if self.stationarity_counter.Tick():
                    self.AssessStationarity()

            # printing if required

            if self.particles_results_counter.Tick():
                self.io_tools.PrintParticlesResults(
                    self.pp.variables_to_print_in_file,
                    self.time,
                    self.spheres_model_part)

                self.PrintDrag(
                    self.drag_list,
                    self.drag_file_output_list,
                    self.fluid_model_part,
                    self.time)

            if self.print_counter_updated_DEM.Tick():

                if coupling_level_type:
                    self.projection_module.ComputePostProcessResults(self.spheres_model_part.ProcessInfo)

                self.post_utils.Writeresults(self.time_dem)

            # solving the DEM part

            self.derivative_recovery_counter.Activate(self.time > interaction_start_time)

            if self.derivative_recovery_counter.Tick():
                self.RecoverDerivatives()

            Say('Solving DEM... (', self.spheres_model_part.NumberOfElements(0), 'elements )')
            first_dem_iter = True

            for self.time_dem in self.yield_DEM_time(
                    self.time_dem,
                    time_final_DEM_substepping,
                    self.Dt_DEM):
                self.DEM_step += 1
                self.spheres_model_part.ProcessInfo[TIME_STEPS] = self.DEM_step
                self.rigid_face_model_part.ProcessInfo[TIME_STEPS] = self.DEM_step
                self.cluster_model_part.ProcessInfo[TIME_STEPS] = self.DEM_step

                self.PerformInitialDEMStepOperations(self.time_dem)

                it_is_time_to_forward_couple = (
                    self.time >= interaction_start_time and
                    coupling_level_type and
                    (project_at_every_substep_option or first_dem_iter)
                )

                if it_is_time_to_forward_couple:

                    if coupling_scheme_type == "UpdatedDEM":
                        self.ApplyForwardCoupling()

                    else:
                        alpha = 1.0 - (time_final_DEM_substepping - self.time_dem) / self.Dt
                        self.ApplyForwardCoupling(alpha)

                        if self.quadrature_counter.Tick():
                            self.AppendValuesForTheHistoryForce()

                        if integration_scheme in {'Hybrid_Bashforth', 'TerminalVelocityScheme'}:
                            # Advance in space only
                            self.DEMSolve(self.time_dem)
                            self.ApplyForwardCouplingOfVelocityToSlipVelocityOnly(self.time_dem)

                # performing the time integration of the DEM part

                self.spheres_model_part.ProcessInfo[TIME] = self.time_dem
                self.rigid_face_model_part.ProcessInfo[TIME] = self.time_dem
                self.cluster_model_part.ProcessInfo[TIME] = self.time_dem

                if self.do_solve_dem:
                    self.DEMSolve(self.time_dem)

                self.disperse_phase_solution.DEMFEMProcedures.MoveAllMeshes(
                    self.all_model_parts,
                    self.time_dem,
                    self.Dt_DEM)

                #### TIME CONTROL ##################################

                # adding DEM elements by the inlet:
                if dem_inlet_option:
                    # After solving, to make sure that neighbours are already set.
                    self.disperse_phase_solution.DEM_inlet.CreateElementsFromInletMesh(
                        self.spheres_model_part,
                        self.cluster_model_part,
                        self.disperse_phase_solution.creator_destructor)

                if self.print_counter_updated_fluid.Tick():

                    if coupling_level_type:
                        self.projection_module.ComputePostProcessResults(
                            self.spheres_model_part.ProcessInfo)

                    self.post_utils.Writeresults(self.time_dem)

                first_dem_iter = False

                # applying DEM-to-fluid coupling

                if self.DEM_to_fluid_counter.Tick() and self.time >= interaction_start_time:
                    self.projection_module.ProjectFromParticles()

                #Phantom
                self.disperse_phase_solution.RunAnalytics(self.time, is_time_to_print=self.analytic_data_counter.Tick())

            #### PRINTING GRAPHS ####
            os.chdir(self.graphs_path)
            # measuring mean velocities in a certain control volume (the 'velocity trap')
            if self.pp.CFD_DEM["VelocityTrapOption"].GetBool():
                self.post_utils_DEM.ComputeMeanVelocitiesinTrap("Average_Velocity.txt", self.time)

            os.chdir(self.post_path)

            # coupling checks (debugging)
            if self.debug_info_counter.Tick():
                self.dem_volume_tool.UpdateDataAndPrint(
                    self.pp.CFD_DEM["fluid_domain_volume"].GetDouble())

            # printing if required

            if self.particles_results_counter.Tick():
                self.io_tools.PrintParticlesResults(
                    self.pp.variables_to_print_in_file,
                    self.time,
                    self.spheres_model_part)

                self.PrintDrag(
                    self.drag_list,
                    self.drag_file_output_list,
                    self.fluid_model_part,
                    self.time)

    def GetFirstStepForFluidComputation(self):
        return 3

    def CloneTimeStep(self):
        self.fluid_model_part.CloneTimeStep(self.time)

    def DEMSolve(self, time='None'): # time is passed in case it is needed
        self.disperse_phase_solution.solver.Solve()

    def UpdateALEMeshMovement(self, time):
        pass

    def RecoverDerivatives(self):
        self.recovery.Recover()

    def FluidSolve(self, time='None', solve_system=True):
        Say('Solving Fluid... (', self.fluid_model_part.NumberOfElements(0), 'elements )\n')

        if solve_system:
            self.fluid_solution.RunSingleTimeStep()
        else:
            Say("Skipping solving system...\n")

    def PerformZeroStepInitializations(self):
        pass

    def PerformInitialDEMStepOperations(self, time=None):
        pass

    def PerformEmbeddedOperations(self):
        # calculating elemental distances defining the structure embedded in the fluid mesh
        if self.pp.CFD_DEM["embedded_option"].GetBool():
            self.calculate_distance_process.Execute()

        if self.embedded_counter.Tick():
            embedded.ApplyEmbeddedBCsToFluid(self.fluid_model_part)
            embedded.ApplyEmbeddedBCsToBalls(self.spheres_model_part, self.pp.CFD_DEM)

    def AssessStationarity(self):
        Say("Assessing Stationarity...\n")
        self.stationarity = self.stationarity_tool.Assess(self.fluid_model_part)
        self.stationarity_counter.Deactivate(self.stationarity)

    def SetInlet(self):
        if self.pp.CFD_DEM["dem_inlet_option"].GetBool():
            # Constructing the inlet and initializing it
            # (must be done AFTER the self.spheres_model_part Initialize)
            # Note that right now only inlets of a single type are possible.
            # This should be generalized.
            if self.pp.type_of_inlet == 'VelocityImposed':
                self.DEM_inlet = DEM_Inlet(self.DEM_inlet_model_part)
            elif self.pp.type_of_inlet == 'ForceImposed':
                self.DEM_inlet = DEM_Force_Based_Inlet(self.DEM_inlet_model_part, self.pp.force)

            self.DEM_inlet.InitializeDEM_Inlet(self.spheres_model_part, self.creator_destructor)

    def SetAnalyticParticleWatcher(self):
        from analytic_tools import analytic_data_procedures
        self.particle_watcher = AnalyticParticleWatcher()
        self.particle_watcher_analyser = analytic_data_procedures.ParticleWatcherAnalyzer(
            analytic_particle_watcher=self.particle_watcher,
            path=self.main_path)

    def ProcessAnalyticData(self):
        self.disperse_phase_solution.WriteAnalyticDataToFileAndClear()

    def SetInletWatcher(self):
        self.watcher_analyser.SetInlet(self.DEM_inlet)

    def TellTime(self, time):
        Say('TIME = ', time)
        Say('ELAPSED TIME = ', self.timer.time() - self.simulation_start_time, '\n')

    def TellFinalSummary(self, time, step, DEM_step):
        simulation_elapsed_time = self.timer.time() - self.simulation_start_time
        if simulation_elapsed_time and step and DEM_step:
            elapsed_time_per_unit_fluid_step = simulation_elapsed_time / step
            elapsed_time_per_unit_DEM_step = simulation_elapsed_time / DEM_step
        else:
            elapsed_time_per_unit_fluid_step = 0.0
            elapsed_time_per_unit_DEM_step = 0.0
        Say('*************************************************************')
        Say('CALCULATIONS FINISHED. THE SIMULATION ENDED SUCCESSFULLY.')

        Say('Elapsed time: ' + '%.5f'%(simulation_elapsed_time) + ' s ')
        Say('per fluid time step: ' + '%.5f'%(elapsed_time_per_unit_fluid_step) + ' s ')
        Say('per DEM time step: ' + '%.5f'%(elapsed_time_per_unit_DEM_step) + ' s ')
        Say('*************************************************************\n')

    def GetFluidSolveCounter(self):
        return SDP.Counter(is_dead=(self.pp.CFD_DEM["drag_force_type"].GetInt() == 9))

    def GetEmbeddedCounter(self):
        # MA: because I think DISTANCE,1 (from previous time step)
        # is not calculated correctly for step=1
        return SDP.Counter(1, 3, self.pp.CFD_DEM["embedded_option"].GetBool())

    def GetBackwardCouplingCounter(self):
        return SDP.Counter(1, 1, self.pp.CFD_DEM["coupling_level_type"].GetInt() > 1)

    def GetRecoveryCounter(self):
        there_is_something_to_recover = (
            self.pp.CFD_DEM["coupling_level_type"].GetInt() or
            self.pp.CFD_DEM["print_PRESSURE_GRADIENT_option"].GetBool())
        return SDP.Counter(1, 1, there_is_something_to_recover)

    def GetStationarityCounter(self):
        return SDP.Counter(
            steps_in_cycle=self.pp.CFD_DEM["time_steps_per_stationarity_step"].GetInt(),
            beginning_step=1,
            is_active=self.pp.CFD_DEM["stationary_problem_option"].GetBool())

    def GetPrintCounterUpdatedDEM(self):
        counter = SDP.Counter(
            steps_in_cycle=int(self.output_time / self.Dt_DEM + 0.5),
            beginning_step=int(self.Dt / self.Dt_DEM))

        if 'UpdatedDEM' != self.pp.CFD_DEM["coupling_scheme_type"].GetString():
            counter.Kill()
        return counter

    def GetPrintCounterUpdatedFluid(self):
        counter = SDP.Counter(
            steps_in_cycle=int(self.output_time / self.Dt_DEM + 0.5),
            beginning_step=int(self.Dt / self.Dt_DEM))

        if 'UpdatedFluid' != self.pp.CFD_DEM["coupling_scheme_type"].GetString():
            counter.Kill()
        return counter

    def GetDebugInfo(self):
        return SDP.Counter(
            self.pp.CFD_DEM["debug_tool_cycle"].GetInt(),
            1,
            self.pp.CFD_DEM["print_debug_info_option"].GetBool())

    def GetParticlesResultsCounter(self):
        return SDP.Counter(
            self.pp.CFD_DEM["print_particles_results_cycle"].GetInt(),
            1,
            self.pp.CFD_DEM["print_particles_results_option"].GetBool())

    def GetHistoryForceQuadratureCounter(self):
        return SDP.Counter(
            self.pp.CFD_DEM["time_steps_per_quadrature_step"].GetInt(),
            1,
            self.pp.CFD_DEM["basset_force_type"].GetInt())

    def ProcessAnalyticDataCounter(self):
        return SDP.Counter(
            steps_in_cycle=self.pp.CFD_DEM["time_steps_per_analytic_processing_step"].GetInt(),
            beginning_step=1,
            is_active=self.pp.CFD_DEM["do_process_analytic_data"].GetBool())

    def GetVolumeDebugTool(self):
        return SDP.ProjectionDebugUtils(
            self.pp.CFD_DEM["fluid_domain_volume"].GetDouble(),
            self.fluid_model_part,
            self.spheres_model_part,
            self.custom_functions_tool)

    def GetRunCode(self):
        return ""
        #return SDP.CreateRunCode(self.pp)

    def FillHistoryForcePrecalculatedVectors(self):
        # Warning: this estimation is based on a constant time step for DEM.
        # This is usually the case, but could not be so.
        # A more robust implementation is needed!
        N_steps = int(self.final_time / self.pp.CFD_DEM["MaxTimeStep"].GetDouble()) + 20
        not_neglecting_history_force = self.pp.CFD_DEM["basset_force_type"].GetInt() > 0
        using_hinsberg_method = (
            self.pp.CFD_DEM["basset_force_type"].GetInt() >= 3 or
            self.pp.CFD_DEM["basset_force_type"].GetInt() == 1)
        if not_neglecting_history_force:
            self.basset_force_tool.FillDaitcheVectors(
                N_steps,
                self.pp.CFD_DEM["quadrature_order"].GetInt(),
                self.pp.CFD_DEM["time_steps_per_quadrature_step"].GetInt())
        if using_hinsberg_method:
            self.basset_force_tool.FillHinsbergVectors(
                self.spheres_model_part,
                self.pp.CFD_DEM["number_of_exponentials"].GetInt(),
                self.pp.CFD_DEM["number_of_quadrature_steps_in_window"].GetInt())

    def AppendValuesForTheHistoryForce(self):
        using_hinsberg_method = (
            self.pp.CFD_DEM["basset_force_type"].GetInt() >= 3 or
            self.pp.CFD_DEM["basset_force_type"].GetInt() == 1)
        if using_hinsberg_method:
            self.basset_force_tool.AppendIntegrandsWindow(self.spheres_model_part)
        elif self.pp.CFD_DEM["basset_force_type"].GetInt() == 2:
            self.basset_force_tool.AppendIntegrands(self.spheres_model_part)

    def GetBassetForceTools(self):
        self.basset_force_tool = SDP.BassetForceTools()

    def GetFieldUtility(self):
        return None

    def GetResultsCreator(self):
        return None

    def ApplyForwardCoupling(self, alpha='None'):
        self.projection_module.ApplyForwardCoupling(alpha)

    def ApplyForwardCouplingOfVelocityToSlipVelocityOnly(self, time=None):
        self.projection_module.ApplyForwardCouplingOfVelocityToSlipVelocityOnly()

    def PerformFinalOperations(self, time=None):
        os.chdir(self.main_path)
        del self.post_utils
        self.ModifyResultsFolderName(time)

    def ModifyResultsFolderName(self, time):
        pass

    def Finalize(self):

        self.swimming_DEM_gid_io.finalize_results()

        self.PerformFinalOperations(self.time_dem)

        self.FinalizeDragOutput()

        self.fluid_solution.Finalize()

        self.TellFinalSummary(self.step, self.time, self.DEM_step)

    def FinalizeDragOutput(self):
        for i in self.drag_file_output_list:
            i.close()

    def SetCutsOutput(self):
        gid_output_options = self.pp.fluid_parameters["output_processes"]["gid_output"][0]["Parameters"]
        result_file_configuration = gid_output_options["postprocess_parameters"]["result_file_configuration"]

        if not result_file_configuration["body_output"].GetBool():
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

        if self.drag_file_output_list:
            Say('Drag output list:', self.drag_file_output_list)

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

    def TransferBodyForceFromDisperseToFluid(self):
        # setting fluid's body force to the same as DEM's
        if self.pp.CFD_DEM["body_force_on_fluid_option"].GetBool():
            body_force = [self.pp.CFD_DEM["GravityX"].GetDouble(),
                          self.pp.CFD_DEM["GravityY"].GetDouble(),
                          self.pp.CFD_DEM["GravityZ"].GetDouble()]
            modulus_of_body_force = math.sqrt(sum([b**2 for b in body_force]))

            gravity_parameters = self.fluid_solution.parameters['processes']['gravity'][0]['Parameters']
            gravity_parameters['modulus'].SetDouble(modulus_of_body_force)
            [gravity_parameters['direction'][i].SetDouble(b) for i, b in enumerate(body_force)]

    def AssignKinematicViscosityFromDynamicViscosity(self):
        # Eulerian fluid already works with kinematic viscosity
        pass

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

    # To-do: for the moment, provided for compatibility
    def _CreateSolver(self):
        return self.fluid_solution._CreateSolver()

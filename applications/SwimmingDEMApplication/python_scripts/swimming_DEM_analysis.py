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
import variables_management

def Say(*args):
    Logger.PrintInfo("SwimmingDEM", *args)
    Logger.Flush()

# Import MPI modules if needed. This way to do this is only valid when using OpenMPI.
# For other implementations of MPI it will not work.
if "OMPI_COMM_WORLD_SIZE" in os.environ:
    # Kratos MPI
    from KratosMultiphysics.MetisApplication import *
    from KratosMultiphysics.MPISearchApplication import *
    from KratosMultiphysics.mpi import *

    # DEM Application MPI
    import DEM_procedures_mpi as DEM_procedures
    # import DEM_material_test_script_mpi as DEM_material_test_script
    Say('Running under MPI...........\n')
else:
    # DEM Application
    import DEM_procedures
    # import DEM_material_test_script
    Say('Running under OpenMP........\n')

class SDEMLogger(object):
    def __init__(self, do_print_file=False):
        self.terminal = sys.stdout
        self.console_output_file_name = 'console_output.txt'
        self.path_to_console_out_file = os.getcwd()
        self.path_to_console_out_file += '/' + self.console_output_file_name
        self.do_print_file = do_print_file
        if self.do_print_file:
            self.log = open(self.path_to_console_out_file, "a")

    def write(self, message):
        self.terminal.write(message)
        if self.do_print_file:
            self.log.write(message)

    def flush(self):
        #this flush method is needed for python 3 compatibility.
        #this handles the flush command by doing nothing.
        #you might want to specify some extra behavior here.
        pass

    def getvalue(self):
        return self.terminal.getvalue()

class python_parameters:
    def __init__(self):
        pass

class SwimmingDEMAnalysis(AnalysisStage):
    def __enter__ (self):
        return self

    def __exit__(self, exception_type, exception_value, traceback):
        pass

    def __init__(self, model, parameters = Parameters("{}")):
        sys.stdout = SDEMLogger()
        self.StartTimer()
        self.model = model
        self.main_path = os.getcwd()
        self.project_parameters = parameters
        self.vars_man = variables_management.VariablesManager(self.project_parameters)

        # storing some frequently used variables
        self.time_step = self.project_parameters["MaxTimeStep"].GetDouble()
        self.end_time   = self.project_parameters["FinalTime"].GetDouble()
        self.do_print_results = self.project_parameters["do_print_results_option"].GetBool()

        self.SetCouplingParameters(parameters)

        self.SetFluidParameters()

        self.ModifyInputParametersForCoherence()

        self.SetDispersePhaseAlgorithm()

        self.disperse_phase_solution.coupling_analysis = weakref.proxy(self)

        self.SetFluidAlgorithm()
        self.fluid_solution.coupling_analysis = weakref.proxy(self)

        self.procedures = weakref.proxy(self.disperse_phase_solution.procedures)
        self.report = DEM_procedures.Report()

        self.disperse_phase_solution.SetAnalyticFaceWatcher()

        # defining member variables for the model_parts (for convenience)
        self.fluid_model_part = self.fluid_solution.fluid_model_part
        self.spheres_model_part = self.disperse_phase_solution.spheres_model_part
        self.cluster_model_part = self.disperse_phase_solution.cluster_model_part
        self.rigid_face_model_part = self.disperse_phase_solution.rigid_face_model_part
        self.DEM_inlet_model_part = self.disperse_phase_solution.DEM_inlet_model_part
        self.vars_man.ConstructListsOfVariables(self.project_parameters)
        super(SwimmingDEMAnalysis, self).__init__(model, self.project_parameters) # TODO: The DEM jason is now interpreted as the coupling json. This must be changed

    def SetFluidParameters(self):
        self.fluid_parameters = self.project_parameters['fluid_parameters']

        self.vars_man.nodal_results, self.vars_man.gauss_points_results = [], []
        if self.fluid_parameters.Has('sdem_output_processes'):
            gid_output_options = self.fluid_parameters["sdem_output_processes"]["gid_output"][0]["Parameters"]
            result_file_configuration = gid_output_options["postprocess_parameters"]["result_file_configuration"]
            gauss_point_results = result_file_configuration["gauss_point_results"]
            nodal_variables = self.fluid_parameters["sdem_output_processes"]["gid_output"][0]["Parameters"]["postprocess_parameters"]["result_file_configuration"]["nodal_results"]
            self.vars_man.nodal_results = [nodal_variables[i].GetString() for i in range(nodal_variables.size())]
            self.vars_man.gauss_points_results = [gauss_point_results[i].GetString() for i in range(gauss_point_results.size())]

    def SetFluidAlgorithm(self):
        import DEM_coupled_fluid_dynamics_analysis
        self.fluid_solution = DEM_coupled_fluid_dynamics_analysis.DEMCoupledFluidDynamicsAnalysis(self.model, self.project_parameters, self.vars_man)
        self.fluid_solution.main_path = self.main_path

    def SetDispersePhaseAlgorithm(self):
        import fluid_coupled_DEM_analysis as DEM_analysis
        self.disperse_phase_solution = DEM_analysis.FluidCoupledDEMAnalysisStage(self.model, self.project_parameters)

    def ReadDispersePhaseAndCouplingParameters(self):
        self.project_parameters = self.project_parameters # TO-DO: remove this

        import swimming_dem_default_input_parameters
        only_swimming_defaults = swimming_dem_default_input_parameters.GetDefaultInputParameters()

        self.project_parameters.ValidateAndAssignDefaults(only_swimming_defaults)

    def SetCouplingParameters(self, parameters):

        # First, read the parameters generated from the interface
        self.ReadDispersePhaseAndCouplingParameters()

        # Second, set the default 'beta' parameters (candidates to be moved to the interface)
        self.SetBetaParameters()

        # Third, set the parameters fed to the particular case that you are running
        self.SetCustomBetaParameters(parameters)

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

    # Set input parameters that have not yet been transferred to the interface
    # import the configuration data as read from the GiD
    def SetBetaParameters(self):
        Add = self.project_parameters.AddEmptyValue
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
        Add("frame_rotation_axis_initial_point").SetVector(Vector([0., 0., 0.]))
        Add("frame_rotation_axis_final_point").SetVector(Vector([0., 0., 1.]))
        Add("angular_velocity_magnitude").SetDouble(1.0)
        if self.project_parameters["type_of_dem_inlet"].GetString() == 'ForceImposed':
            Add("inlet_force_vector").SetVector(Vector([0., 0., 1.])) # TODO: generalize

        # Setting body_force_per_unit_mass_variable_name
        Add("body_force_per_unit_mass_variable_name").SetString('BODY_FORCE')

    # This step is added to allow modifications to the possibly incompatibilities
    # between the individual parameters coming from each sub-application
    # (i.e., fluid and dem apps)
    def ModifyInputParametersForCoherence(self):
        # Making all time steps exactly commensurable
        output_time = self.project_parameters["OutputTimeStep"].GetDouble()
        self.output_time = int(output_time / self.time_step) * self.time_step
        self.project_parameters["OutputTimeStep"].SetDouble(self.output_time)
        self.fluid_time_step = self.fluid_parameters["solver_settings"]["time_stepping"]["time_step"].GetDouble()
        self.fluid_time_step = int(self.fluid_time_step / self.time_step) * self.time_step
        self.fluid_parameters["solver_settings"]["time_stepping"]["time_step"].SetDouble(self.fluid_time_step)
        self.project_parameters["dem_parameters"]["MaxTimeStep"].SetDouble(self.time_step)

        if self.project_parameters["coupling_level_type"].GetInt() == 0:
            self.project_parameters["project_at_every_substep_option"].SetBool(False)

        # the fluid fraction is not projected from DEM (there may not
        # be a DEM part) but is externally imposed:
        if self.project_parameters["flow_in_porous_medium_option"].GetBool():
            self.project_parameters["coupling_weighing_type"].SetInt(- 1)

        time_steps_per_stationarity_step = self.project_parameters["time_steps_per_stationarity_step"].GetInt()
        self.project_parameters["time_steps_per_stationarity_step"].SetInt(max(1, int(time_steps_per_stationarity_step)))

        if self.project_parameters["coupling_level_type"].GetInt() > 1:
            self.project_parameters["stationary_problem_option"].SetBool(False)

        self.SetDoSolveDEMVariable()

        self.TransferBodyForceFromDisperseToFluid()

    def TransferBodyForceFromDisperseToFluid(self):
        # setting fluid's body force to the same as DEM's
        if self.project_parameters["body_force_on_fluid_option"].GetBool():
            body_force = [self.project_parameters["GravityX"].GetDouble(),
                          self.project_parameters["GravityY"].GetDouble(),
                          self.project_parameters["GravityZ"].GetDouble()]
            modulus_of_body_force = math.sqrt(sum(b**2 for b in body_force))

            gravity_parameters = self.fluid_parameters['processes']['gravity'][0]['Parameters']
            gravity_parameters['modulus'].SetDouble(modulus_of_body_force)
            for i, b in enumerate(body_force):
                gravity_parameters['direction'][i].SetDouble(b)

    def SetDoSolveDEMVariable(self):
        self.do_solve_dem = self.project_parameters["do_solve_dem"].GetBool()

        if self.project_parameters["flow_in_porous_DEM_medium_option"].GetBool():
            self.do_solve_dem = False

    def SetCustomBetaParameters(self, custom_parameters):
        custom_parameters.ValidateAndAssignDefaults(self.project_parameters)
        self.project_parameters = custom_parameters
        # TO DO: remove next lines as info is taken from Parameters object everywhere
        # var_names = [k for k in dictionary.keys()]
        # var_values = [k for k in dictionary.values()]
        # for name, value in zip(var_names, var_values):
        #     self.project_parameters.__setitem__(name, value)

    def Run(self):
        super(SwimmingDEMAnalysis, self).Run()

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
        Say('Initializing simulation...\n')
        self.run_code = self.GetRunCode()

        # Moving to the recently created folder
        os.chdir(self.main_path)
        if self.do_print_results:
            [self.post_path, data_and_results, self.graphs_path, MPI_results] = \
            self.procedures.CreateDirectories(str(self.main_path),
                                            str(self.project_parameters["problem_name"].GetString()),
                                            self.run_code)
            SDP.CopyInputFilesIntoFolder(self.main_path, self.post_path)
            self.MPI_results = MPI_results

        self.FluidInitialize()

        self.DispersePhaseInitialize()

        self.SetAllModelParts()

        if self.project_parameters.Has('sdem_output_processes') and self.do_print_results:
            gid_output_options = self.project_parameters["sdem_output_processes"]["gid_output"][0]["Parameters"]
            result_file_configuration = gid_output_options["postprocess_parameters"]["result_file_configuration"]
            write_conditions_option = result_file_configuration["gidpost_flags"]["WriteConditionsFlag"].GetString() == "WriteConditions"
            deformed_mesh_option = result_file_configuration["gidpost_flags"]["WriteDeformedMeshFlag"].GetString() == "WriteDeformed"
            old_gid_output_post_options_dict = {'GiD_PostAscii':'Ascii','GiD_PostBinary':'Binary','GiD_PostAsciiZipped':'AsciiZipped'}
            old_gid_output_multiple_file_option_dict = {'SingleFile':'Single','MultipleFiles':'Multiples'}
            post_mode_key = result_file_configuration["gidpost_flags"]["GiDPostMode"].GetString()
            multiple_files_option_key = result_file_configuration["gidpost_flags"]["MultiFileFlag"].GetString()

            self.swimming_DEM_gid_io = \
            swimming_DEM_gid_output.SwimmingDEMGiDOutput(
                file_name = self.project_parameters["problem_name"].GetString(),
                vol_output = result_file_configuration["body_output"].GetBool(),
                post_mode = old_gid_output_post_options_dict[post_mode_key],
                multifile = old_gid_output_multiple_file_option_dict[multiple_files_option_key],
                deformed_mesh = deformed_mesh_option,
                write_conditions = write_conditions_option)

            self.swimming_DEM_gid_io.initialize_swimming_DEM_results(
                self.spheres_model_part,
                self.cluster_model_part,
                self.rigid_face_model_part,
                self.mixed_model_part)

        self.SetPointGraphPrinter()

        self.AssignKinematicViscosityFromDynamicViscosity()

        super(SwimmingDEMAnalysis, self).Initialize()

        # coarse-graining: applying changes to the physical properties of the model to adjust for
        # the similarity transformation if required (fluid effects only).
        SDP.ApplySimilarityTransformations(
            self.fluid_model_part,
            self.project_parameters["similarity_transformation_type"].GetInt(),
            self.project_parameters["model_over_real_diameter_factor"].GetDouble()
            )

        if self.do_print_results:
            self.SetPostUtils()

        # creating an IOTools object to perform other printing tasks
        self.io_tools = SDP.IOTools(self.project_parameters)

        dem_physics_calculator = SphericElementGlobalPhysicsCalculator(
            self.spheres_model_part)

        if self.project_parameters["coupling_level_type"].GetInt():
            default_meso_scale_length_needed = (
                self.project_parameters["meso_scale_length"].GetDouble() <= 0.0 and
                self.spheres_model_part.NumberOfElements(0) > 0)

            if default_meso_scale_length_needed:
                biggest_size = (2 * dem_physics_calculator.CalculateMaxNodalVariable(self.spheres_model_part, RADIUS))
                self.project_parameters["meso_scale_length"].SetDouble(20 * biggest_size)

            elif self.spheres_model_part.NumberOfElements(0) == 0:
                self.project_parameters["meso_scale_length"].SetDouble(1.0)

        # creating a custom functions calculator for the implementation of
        # additional custom functions
        fluid_domain_dimension = self.project_parameters["fluid_parameters"]["solver_settings"]["domain_size"].GetInt()
        self.custom_functions_tool = SDP.FunctionsCalculator(fluid_domain_dimension)

        # creating a stationarity assessment tool
        self.stationarity_tool = SDP.StationarityAssessmentTool(
            self.project_parameters["max_pressure_variation_rate_tol"].GetDouble(),
            self.custom_functions_tool
            )

        # creating a debug tool
        self.dem_volume_tool = self.GetVolumeDebugTool()

        #self.SetEmbeddedTools()

        Say('Initialization Complete\n')

        if self.project_parameters["flow_in_porous_DEM_medium_option"].GetBool():
            SDP.FixModelPart(self.spheres_model_part)

        ##################################################

        #    I N I T I A L I Z I N G    T I M E    L O O P

        ##################################################
        self.step = 0
        self.time = self.fluid_parameters["problem_data"]["start_time"].GetDouble()
        self.fluid_time_step = self.fluid_solution._GetSolver()._ComputeDeltaTime()
        self.time_step = self.spheres_model_part.ProcessInfo.GetValue(DELTA_TIME)
        self.rigid_face_model_part.ProcessInfo[DELTA_TIME] = self.time_step
        self.cluster_model_part.ProcessInfo[DELTA_TIME] = self.time_step
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
        # Phantom
        self.disperse_phase_solution.analytic_data_counter = self.ProcessAnalyticDataCounter()
        self.mat_deriv_averager           = SDP.Averager(1, 3)
        self.laplacian_averager           = SDP.Averager(1, 3)

        self.report.total_steps_expected = int(self.end_time / self.time_step)

        Say(self.report.BeginReport(self.timer))

        # creating a Post Utils object that executes several post-related tasks
        self.post_utils_DEM = DEM_procedures.PostUtils(self.project_parameters['dem_parameters'], self.spheres_model_part)

        # otherwise variables are set to 0 by default:
        SDP.InitializeVariablesWithNonZeroValues(self.project_parameters,
                                                 self.fluid_model_part,
                                                 self.spheres_model_part)

        if self.do_print_results:
            self.SetUpResultsDatabase()

        # ANALYTICS BEGIN
        self.project_parameters.AddEmptyValue("perform_analytics_option").SetBool(False)

        if self.project_parameters["perform_analytics_option"].GetBool():
            import analytics
            variables_to_measure = [PRESSURE]
            steps_between_measurements = 100
            gauge = analytics.Gauge(
                self.fluid_model_part,
                self.fluid_time_step,
                self.end_time,
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
            self.project_parameters,
            self.fluid_model_part,
            self.custom_functions_tool)

        self.FillHistoryForcePrecalculatedVectors()

        self.PerformZeroStepInitializations()

        if self.do_print_results:
            self._Print()

    def AddExtraProcessInfoVariablesToFluid(self):
        self.vars_man.AddExtraProcessInfoVariablesToFluidModelPart(self.project_parameters, self.fluid_model_part)

    def FluidInitialize(self):
        self.fluid_model_part = self.fluid_solution.fluid_model_part
        self.fluid_solution.vars_man = self.vars_man
        self.fluid_solution.Initialize()

        self.AddExtraProcessInfoVariablesToFluid()

        SDP.AddExtraDofs(self.fluid_model_part,
                         self.spheres_model_part,
                         self.cluster_model_part,
                         self.DEM_inlet_model_part,
                         self.vars_man)

    def DispersePhaseInitialize(self):
        self.vars_man.__class__.AddNodalVariables(self.spheres_model_part, self.vars_man.dem_vars)
        self.vars_man.__class__.AddNodalVariables(self.rigid_face_model_part, self.vars_man.rigid_faces_vars)
        self.vars_man.__class__.AddNodalVariables(self.DEM_inlet_model_part, self.vars_man.inlet_vars)
        self.vars_man.AddExtraProcessInfoVariablesToDispersePhaseModelPart(self.project_parameters,
                                                                           self.disperse_phase_solution.spheres_model_part)

        self.disperse_phase_solution.Initialize()

    def SetPostUtils(self):
          # creating a Post Utils object that executes several post-related tasks
        self.post_utils = SDP.PostUtils(self.swimming_DEM_gid_io,
                                        self.project_parameters,
                                        self.vars_man,
                                        self.fluid_model_part,
                                        self.spheres_model_part,
                                        self.cluster_model_part,
                                        self.rigid_face_model_part,
                                        self.mixed_model_part)

    def SetEmbeddedTools(self):
    # creating a distance calculation process for the embedded technology
        # (used to calculate elemental distances defining the structure embedded in the fluid mesh)
        if self.project_parameters["embedded_option"].GetBool():
            self.calculate_distance_process = CalculateSignedDistanceTo3DSkinProcess(
                self.rigid_face_model_part,
                self.fluid_model_part
                )
            self.calculate_distance_process.Execute()

    def TheSimulationMustGoOn(self):
        return self.time <= self.end_time

    def GetAnalyticFacesModelParts(self):
        analytic_face_submodelpart_number = 1
        analytic_face_submodelpart_name = self.rigid_face_model_part.GetSubModelPart(str(analytic_face_submodelpart_number))
        return analytic_face_submodelpart_name

    def MakeAnalyticsMeasurements(self):
        self.analytic_face_watcher.MakeMeasurements()
        self.analytic_particle_watcher.MakeMeasurements()

    def RunSolutionLoop(self):
        while self.TheSimulationMustGoOn():
            self.step, self.time = self._GetSolver().AdvanceInTime(self.step, self.time)
            self.InitializeSolutionStep()
            self._GetSolver().Predict()
            self._GetSolver().SolveSolutionStep()
            self.FinalizeSolutionStep()
            self.OutputSolutionStep()

    def InitializeSolutionStep(self):
        self.TellTime()
        self.PerformInitialDEMStepOperations(self.time)
        self.disperse_phase_solution.InitializeSolutionStep()
        if self._GetSolver().CannotIgnoreFluidNow():
            self.fluid_solution.InitializeSolutionStep()
        super(SwimmingDEMAnalysis, self).InitializeSolutionStep()

    def FinalizeSolutionStep(self):
        # printing if required
        if self._GetSolver().CannotIgnoreFluidNow():
            self.fluid_solution.FinalizeSolutionStep()

        self.disperse_phase_solution.FinalizeSolutionStep()

        # applying DEM-to-fluid coupling

        if self.DEM_to_fluid_counter.Tick() and self.time >= self.project_parameters["interaction_start_time"].GetDouble():
            self._GetSolver().projection_module.ProjectFromParticles()

        # coupling checks (debugging)
        if self.debug_info_counter.Tick():
            self.dem_volume_tool.UpdateDataAndPrint(
                self.project_parameters["fluid_domain_volume"].GetDouble())

        super(SwimmingDEMAnalysis, self).FinalizeSolutionStep()

    def OutputSolutionStep(self):
        # printing if required

        if self.print_counter_updated_fluid.Tick():
            self.ComputePostProcessResults()
            self._Print()

        super(SwimmingDEMAnalysis, self).OutputSolutionStep()

    def _Print(self):
        os.chdir(self.post_path)
        import define_output
        self.drag_list = define_output.DefineDragList()
        self.drag_file_output_list = []

        if self.particles_results_counter.Tick():
            self.io_tools.PrintParticlesResults(
                self.vars_man.variables_to_print_in_file,
                self.time,
                self.spheres_model_part)

        self.post_utils.Writeresults(self.time)
        os.chdir(self.main_path)

    def ComputePostProcessResults(self):
        if self.project_parameters["coupling_level_type"].GetInt():
            self._GetSolver().projection_module.ComputePostProcessResults(self.spheres_model_part.ProcessInfo)

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
        if self.project_parameters["embedded_option"].GetBool():
            self.calculate_distance_process.Execute()

        if self.embedded_counter.Tick():
            embedded.ApplyEmbeddedBCsToFluid(self.fluid_model_part)
            embedded.ApplyEmbeddedBCsToBalls(self.spheres_model_part, self.project_parameters)

    def AssessStationarity(self):
        Say("Assessing Stationarity...\n")
        self.stationarity = self.stationarity_tool.Assess(self.fluid_model_part)
        self.stationarity_counter.Deactivate(self.stationarity)

    def SetInlet(self):
        if self.project_parameters["dem_inlet_option"].GetBool():
            # Constructing the inlet and initializing it
            # (must be done AFTER the self.spheres_model_part Initialize)
            # Note that right now only inlets of a single type are possible.
            # This should be generalized.
            if self.project_parameters["type_of_inlet"].GetString() == 'VelocityImposed':
                self.DEM_inlet = DEM_Inlet(self.DEM_inlet_model_part)
            elif self.project_parameters["type_of_inlet"].GetString() == 'ForceImposed':
                self.DEM_inlet = DEM_Force_Based_Inlet(self.DEM_inlet_model_part, self.project_parameters["inlet_force_vector"].GetVector())

            self.disperse_phase_solution.DEM_inlet = self.DEM_inlet
            self.DEM_inlet.InitializeDEM_Inlet(self.spheres_model_part, self.disperse_phase_solution.creator_destructor)

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

    def TellTime(self):
        Say('DEM time: ', str(self.time) + ', step: ', self.step)
        Say('fluid time: ', str(self._GetSolver().next_time_to_solve_fluid) + ', step: ', self._GetSolver().fluid_step)
        Say('ELAPSED TIME = ', self.timer.time() - self.simulation_start_time, '\n')

    def TellFinalSummary(self, time, dem_step, fluid_step, message_n_char_width=60):
        simulation_elapsed_time = self.timer.time() - self.simulation_start_time

        if simulation_elapsed_time and dem_step and fluid_step:
            elapsed_time_per_unit_dem_step = simulation_elapsed_time / dem_step
            elapsed_time_per_unit_fluid_step = simulation_elapsed_time / fluid_step

        else:
            elapsed_time_per_unit_dem_step = 0.0
            elapsed_time_per_unit_fluid_step = 0.0

        final_message = ('\n\n'
                         + '*' * message_n_char_width + '\n'
                         + 'CALCULATIONS FINISHED. THE SIMULATION ENDED SUCCESSFULLY.' + '\n'
                         + 'Total number of DEM steps run: ' + str(dem_step) + '\n'
                         + 'Total number of fluid steps run: ' + str(fluid_step) + '\n'
                         + 'Elapsed time: ' + '%.5f'%(simulation_elapsed_time) + ' s ' + '\n'
                         + ',, per fluid time step: ' + '%.5f'%(elapsed_time_per_unit_fluid_step) + ' s ' + '\n'
                         + ',, per DEM time step: ' + '%.5f'%(elapsed_time_per_unit_dem_step) + ' s ' + '\n'
                         + '*' * message_n_char_width + '\n')

        Say(final_message)

    def GetFluidSolveCounter(self):
        return SDP.Counter(is_dead=(self.project_parameters["drag_force_type"].GetInt() == 9))

    def GetBackwardCouplingCounter(self):
        return SDP.Counter(1, 1, self.project_parameters["coupling_level_type"].GetInt() > 1)

    def GetRecoveryCounter(self):
        there_is_something_to_recover = (
            self.project_parameters["coupling_level_type"].GetInt() or
            self.project_parameters["print_PRESSURE_GRADIENT_option"].GetBool())
        return SDP.Counter(1, 1, there_is_something_to_recover)

    def GetStationarityCounter(self):
        return SDP.Counter(
            steps_in_cycle=self.project_parameters["time_steps_per_stationarity_step"].GetInt(),
            beginning_step=1,
            is_active=self.project_parameters["stationary_problem_option"].GetBool())

    def GetPrintCounterUpdatedDEM(self):
        counter = SDP.Counter(
            steps_in_cycle=int(self.output_time / self.time_step + 0.5),
            beginning_step=int(self.fluid_time_step / self.time_step))

        if 'UpdatedDEM' != self.project_parameters["coupling_scheme_type"].GetString():
            counter.Kill()
        return counter

    def GetPrintCounterUpdatedFluid(self):
        counter = SDP.Counter(
            steps_in_cycle=int(self.output_time / self.time_step + 0.5),
            beginning_step=int(self.output_time / self.time_step),
            is_dead = not self.do_print_results)

        if 'UpdatedFluid' != self.project_parameters["coupling_scheme_type"].GetString():
            counter.Kill()
        return counter

    def GetDebugInfo(self):
        return SDP.Counter(
            self.project_parameters["debug_tool_cycle"].GetInt(),
            1,
            self.project_parameters["print_debug_info_option"].GetBool())

    def GetParticlesResultsCounter(self):
        return SDP.Counter(
            self.project_parameters["print_particles_results_cycle"].GetInt(),
            1,
            self.project_parameters["print_particles_results_option"].GetBool())

    def GetHistoryForceQuadratureCounter(self):
        for prop in self.project_parameters["properties"].values():
            if prop["hydrodynamic_law_parameters"].Has("history_force_parameters"):
                history_force_parameters =  prop["hydrodynamic_law_parameters"]["history_force_parameters"]
                if history_force_parameters.Has("time_steps_per_quadrature_step"):
                    time_steps_per_quadrature_step = history_force_parameters["time_steps_per_quadrature_step"].GetInt()

                    return SDP.Counter(steps_in_cycle=time_steps_per_quadrature_step, beginning_step=1)

        return SDP.Counter(is_dead=True)

    def ProcessAnalyticDataCounter(self):
        return SDP.Counter(
            steps_in_cycle=self.project_parameters["time_steps_per_analytic_processing_step"].GetInt(),
            beginning_step=1,
            is_active=self.project_parameters["do_process_analytic_data"].GetBool())

    def GetVolumeDebugTool(self):
        return SDP.ProjectionDebugUtils(
            self.project_parameters["fluid_domain_volume"].GetDouble(),
            self.fluid_model_part,
            self.spheres_model_part,
            self.custom_functions_tool)

    def GetRunCode(self):
        return ""
        #return SDP.CreateRunCode(self.project_parameters)

    def FillHistoryForcePrecalculatedVectors(self): # TODO: more robust implementation
        # Warning: this estimation is based on a constant time step for DEM.
        # This is usually the case, but could not be so.
        for prop in self.project_parameters["properties"].values():
            if prop["hydrodynamic_law_parameters"].Has("history_force_parameters"):

                if prop["hydrodynamic_law_parameters"]["history_force_parameters"]["name"].GetString() != 'default':
                    total_number_of_steps = int(self.end_time / self.project_parameters["MaxTimeStep"].GetDouble()) + 20
                    history_force_parameters = prop["hydrodynamic_law_parameters"]["history_force_parameters"]
                    time_steps_per_quadrature_step = history_force_parameters["time_steps_per_quadrature_step"].GetInt()
                    self._GetSolver().basset_force_tool.FillDaitcheVectors(
                        total_number_of_steps,
                        history_force_parameters["quadrature_order"].GetInt(),
                        time_steps_per_quadrature_step)

                    if history_force_parameters.Has("mae_parameters"):
                        mae_parameters = history_force_parameters["mae_parameters"]
                        time_window = mae_parameters["window_time_interval"].GetDouble()
                        quadrature_dt = time_steps_per_quadrature_step * self.time_step
                        number_of_quadrature_steps_in_window = int(time_window / quadrature_dt)
                        if mae_parameters["do_use_mae"].GetBool():
                            self._GetSolver().basset_force_tool.FillHinsbergVectors(
                            self.spheres_model_part,
                            mae_parameters["m"].GetInt(),
                            number_of_quadrature_steps_in_window)
                            break

    def GetFieldUtility(self):
        return None

    def ApplyForwardCoupling(self, alpha='None'):
        self._GetSolver().projection_module.ApplyForwardCoupling(alpha)

    def ApplyForwardCouplingOfVelocityToSlipVelocityOnly(self, time=None):
        self._GetSolver().projection_module.ApplyForwardCouplingOfVelocityToSlipVelocityOnly()

    def PerformFinalOperations(self, time=None):
        os.chdir(self.main_path)

        if self.do_print_results:
            del self.post_utils
            self.ModifyResultsFolderName(time)

    def ModifyResultsFolderName(self, time):
        pass

    def Finalize(self):
        Say('Finalizing simulation...\n')
        if self.do_print_results:
            self.swimming_DEM_gid_io.finalize_results()

        self.PerformFinalOperations(self.time)

        self.fluid_solution.Finalize()

        self.TellFinalSummary(self.time, self.step, self._GetSolver().fluid_step)

    def SetPointGraphPrinter(self):
        pass

    def AssignKinematicViscosityFromDynamicViscosity(self):
        # Eulerian fluid already works with kinematic viscosity
        pass

    def GetReturnValue(self):
        return 0.0

    # To-do: for the moment, provided for compatibility
    def _CreateSolver(self):
        import swimming_DEM_solver
        return swimming_DEM_solver.SwimmingDEMSolver(self.model,
                                                     self.project_parameters,
                                                     self.GetFieldUtility(),
                                                     self.fluid_solution._GetSolver(),
                                                     self.disperse_phase_solution._GetSolver(),
                                                     self.vars_man)

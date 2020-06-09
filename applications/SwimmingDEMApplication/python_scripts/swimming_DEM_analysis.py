from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

#TODO: test DEM bounding box

import os
import sys
import math
import time as timer
import weakref

import KratosMultiphysics as Kratos
from KratosMultiphysics import Vector, Logger, Parameters
import KratosMultiphysics.DEMApplication as DEM

from KratosMultiphysics.analysis_stage import AnalysisStage

import KratosMultiphysics.SwimmingDEMApplication.swimming_DEM_procedures as SDP
import KratosMultiphysics.SwimmingDEMApplication.swimming_DEM_gid_output as swimming_DEM_gid_output
import KratosMultiphysics.SwimmingDEMApplication.embedded as embedded
import KratosMultiphysics.SwimmingDEMApplication.variables_management as variables_management

def Say(*args):
    Logger.PrintInfo("SwimmingDEM", *args)
    Logger.Flush()

# Import MPI modules if needed. This way to do this is only valid when using OpenMPI.
# For other implementations of MPI it will not work.

import KratosMultiphysics.DEMApplication.DEM_procedures as DP

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

        self.SetProjectParameters(parameters)

        self.vars_man = variables_management.VariablesManager(self.project_parameters)

        self._GetDEMAnalysis().coupling_analysis = weakref.proxy(self)

        self._GetFluidAnalysis().coupling_analysis = weakref.proxy(self)

        self.procedures = weakref.proxy(self._GetDEMAnalysis().procedures)

        self.report = DP.Report()

        self._GetDEMAnalysis().SetAnalyticFaceWatcher()

        # defining member variables for the model_parts (for convenience)
        self.fluid_model_part = self._GetFluidAnalysis().fluid_model_part
        self.spheres_model_part = self._GetDEMAnalysis().spheres_model_part
        self.cluster_model_part = self._GetDEMAnalysis().cluster_model_part
        self.rigid_face_model_part = self._GetDEMAnalysis().rigid_face_model_part
        self.dem_inlet_model_part = self._GetDEMAnalysis().dem_inlet_model_part
        self.vars_man.ConstructListsOfVariables(self.project_parameters)

        super(SwimmingDEMAnalysis, self).__init__(model, self.project_parameters)

    def SetFluidParameters(self):
        pass

    def SetProjectParameters(self, parameters):
        self.project_parameters = parameters
        self.time_step = self.project_parameters["time_stepping"]["time_step"].GetDouble()
        self.end_time   = self.project_parameters["problem_data"]["end_time"].GetDouble()
        self.do_print_results = self.project_parameters["do_print_results_option"].GetBool()
        self.fluid_parameters = self.project_parameters['fluid_parameters']

        # First, read the parameters generated from the interface
        import KratosMultiphysics.SwimmingDEMApplication.swimming_dem_default_input_parameters as only_swimming_defaults
        import KratosMultiphysics.DEMApplication.dem_default_input_parameters as dem_defaults


        self.project_parameters.ValidateAndAssignDefaults(only_swimming_defaults.GetDefaultInputParameters())
        self.project_parameters["dem_parameters"].ValidateAndAssignDefaults(dem_defaults.GetDefaultInputParameters())

        # Second, set the default 'beta' parameters (candidates to be moved to the interface)
        self.SetBetaParameters()

        # Third, make sure the parameters passed to the different (sub-)analyses is coherent
        # with the general parameters
        self.ModifyInputParametersForCoherence()

    def SetAllModelParts(self):
        self.all_model_parts = weakref.proxy(self._GetDEMAnalysis().all_model_parts)

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
        if self.project_parameters["custom_dem"]["type_of_dem_inlet"].GetString() == 'ForceImposed':
            Add("inlet_force_vector").SetVector(Vector([0., 0., 1.])) # TODO: generalize

        # Setting body_force_per_unit_mass_variable_name
        Add("body_force_per_unit_mass_variable_name").SetString('BODY_FORCE')
        self.project_parameters["dem_parameters"].AddEmptyValue("do_print_results_option").SetBool(self.do_print_results)

    # This step is added to allow modifications to the possibly incompatibilities
    # between the individual parameters coming from each sub-application
    # (i.e., fluid and dem apps)
    def ModifyInputParametersForCoherence(self):
        # Making all time steps exactly commensurable
        output_time = self.project_parameters["output_interval"].GetDouble()
        self.output_time = int(output_time / self.time_step) * self.time_step
        self.project_parameters["output_interval"].SetDouble(self.output_time)

        self.fluid_time_step = self.fluid_parameters["solver_settings"]["time_stepping"]["time_step"].GetDouble()

        if self.fluid_time_step < self.time_step:
            error_message = ('The fluid time step (' + str(self.fluid_time_step)
                             + ') must be larger or equal than the overall time step (' + str(self.time_step)
                             + ')!')
            raise Exception(error_message)

        self.fluid_time_step = int(self.fluid_time_step / self.time_step) * self.time_step
        self.fluid_parameters["solver_settings"]["time_stepping"]["time_step"].SetDouble(self.fluid_time_step)
        self.project_parameters["dem_parameters"]["MaxTimeStep"].SetDouble(self.time_step)
        translational_scheme_name = self.project_parameters["custom_dem"]["translational_integration_scheme"].GetString()
        self.project_parameters["dem_parameters"]["TranslationalIntegrationScheme"].SetString(translational_scheme_name)
        # The fluid fraction is not projected from DEM (there may not
        # be a DEM part) but is externally imposed instead:
        if self.project_parameters["custom_fluid"]["flow_in_porous_medium_option"].GetBool():
            self.project_parameters["coupling"]["coupling_weighing_type"].SetInt(- 1)

        time_steps_per_stationarity_step = self.project_parameters["stationarity"]["time_steps_per_stationarity_step"].GetInt()
        self.project_parameters["stationarity"]["time_steps_per_stationarity_step"].SetInt(max(1, int(time_steps_per_stationarity_step)))

        if self.project_parameters["coupling"]["coupling_level_type"].GetInt() > 1:
            self.project_parameters["stationarity"]["stationary_problem_option"].SetBool(False)

        self.SetDoSolveDEMVariable()

        self.TransferBodyForceFromDisperseToFluid()

    def TransferBodyForceFromDisperseToFluid(self):
        # setting fluid's body force to the same as DEM's
        if self.project_parameters["custom_fluid"]["body_force_on_fluid_option"].GetBool():
            gravity = self.project_parameters["gravity_parameters"]["direction"].GetVector()
            gravity *= self.project_parameters["gravity_parameters"]["modulus"].GetDouble()
            modulus_of_body_force = math.sqrt(sum(b**2 for b in gravity))

            gravity_parameters = self.fluid_parameters['processes']['gravity'][0]['Parameters']
            gravity_parameters['modulus'].SetDouble(modulus_of_body_force)
            gravity_parameters['direction'].SetVector(gravity)

    def SetDoSolveDEMVariable(self):
        self.do_solve_dem = self.project_parameters["custom_dem"]["do_solve_dem"].GetBool()

        if self.project_parameters["custom_fluid"]["flow_in_porous_DEM_medium_option"].GetBool():
            self.do_solve_dem = False

    def Run(self):
        super(SwimmingDEMAnalysis, self).Run()

        return self.GetReturnValue()

    def SetUpResultsDatabase(self):
        pass

    def ReadDispersePhaseModelParts(self,
                                    starting_node_Id=0,
                                    starting_elem_Id=0,
                                    starting_cond_Id=0):
        creator_destructor = self._GetDEMAnalysis().creator_destructor
        max_node_Id = creator_destructor.FindMaxNodeIdInModelPart(self.fluid_model_part)
        max_elem_Id = creator_destructor.FindMaxElementIdInModelPart(self.fluid_model_part)
        max_cond_Id = creator_destructor.FindMaxConditionIdInModelPart(self.fluid_model_part)
        self._GetDEMAnalysis().BaseReadModelParts(max_node_Id, max_elem_Id, max_cond_Id)

    def Initialize(self):
        Say('Initializing simulation...\n')
        self.run_code = self.GetRunCode()

        # Moving to the recently created folder
        os.chdir(self.main_path)
        if self.do_print_results:
            [self.post_path, data_and_results, self.graphs_path, MPI_results] = \
            self.procedures.CreateDirectories(str(self.main_path),
                                            str(self.project_parameters["problem_data"]["problem_name"].GetString()),
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
                file_name = self.project_parameters["problem_data"]["problem_name"].GetString(),
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
            self.project_parameters["similarity"]["similarity_transformation_type"].GetInt(),
            self.project_parameters["similarity"]["model_over_real_diameter_factor"].GetDouble()
            )

        if self.do_print_results:
            self.SetPostUtils()

        # creating an IOTools object to perform other printing tasks
        self.io_tools = SDP.IOTools(self.project_parameters)

        dem_physics_calculator = DEM.SphericElementGlobalPhysicsCalculator(
            self.spheres_model_part)

        if self.project_parameters["coupling"]["coupling_level_type"].GetInt():
            default_meso_scale_length_needed = (
                self.project_parameters["coupling"]["backward_coupling"]["meso_scale_length"].GetDouble() <= 0.0 and
                self.spheres_model_part.NumberOfElements(0) > 0)

            if default_meso_scale_length_needed:
                biggest_size = (2 * dem_physics_calculator.CalculateMaxNodalVariable(self.spheres_model_part, Kratos.RADIUS))
                self.project_parameters["coupling"]["backward_coupling"]["meso_scale_length"].SetDouble(20 * biggest_size)

            elif self.spheres_model_part.NumberOfElements(0) == 0:
                self.project_parameters["coupling"]["backward_coupling"]["meso_scale_length"].SetDouble(1.0)

        # creating a custom functions calculator for the implementation of
        # additional custom functions
        fluid_domain_dimension = self.project_parameters["fluid_parameters"]["solver_settings"]["domain_size"].GetInt()
        self.custom_functions_tool = SDP.FunctionsCalculator(fluid_domain_dimension)

        # creating a debug tool
        self.dem_volume_tool = self.GetVolumeDebugTool()

        #self.SetEmbeddedTools()

        Say('Initialization Complete\n')

        if self.project_parameters["custom_fluid"]["flow_in_porous_DEM_medium_option"].GetBool():
            SDP.FixModelPart(self.spheres_model_part)

        ##################################################

        #    I N I T I A L I Z I N G    T I M E    L O O P

        ##################################################
        self.step = 0
        self.time = self.fluid_parameters["problem_data"]["start_time"].GetDouble()
        self.fluid_time_step = self._GetFluidAnalysis()._GetSolver()._ComputeDeltaTime()
        self.time_step = self.spheres_model_part.ProcessInfo.GetValue(Kratos.DELTA_TIME)
        self.rigid_face_model_part.ProcessInfo[Kratos.DELTA_TIME] = self.time_step
        self.cluster_model_part.ProcessInfo[Kratos.DELTA_TIME] = self.time_step
        self.stationarity = False

        # setting up loop counters:
        self.DEM_to_fluid_counter = self.GetBackwardCouplingCounter()
        self.stationarity_counter = self.GetStationarityCounter()
        self.print_counter = self.GetPrintCounter()
        self.debug_info_counter = self.GetDebugInfo()
        self.particles_results_counter = self.GetParticlesResultsCounter()
        self.quadrature_counter = self.GetHistoryForceQuadratureCounter()
        # Phantom
        self._GetDEMAnalysis().analytic_data_counter = self.ProcessAnalyticDataCounter()
        self.mat_deriv_averager           = SDP.Averager(1, 3)
        self.laplacian_averager           = SDP.Averager(1, 3)

        self.report.total_steps_expected = int(self.end_time / self.time_step)

        Say(self.report.BeginReport(self.timer))

        # creating a Post Utils object that executes several post-related tasks
        self.post_utils_DEM = DP.PostUtils(self.project_parameters['dem_parameters'], self.spheres_model_part)

        # otherwise variables are set to 0 by default:
        SDP.InitializeVariablesWithNonZeroValues(self.project_parameters,
                                                 self.fluid_model_part,
                                                 self.spheres_model_part)

        if self.do_print_results:
            self.SetUpResultsDatabase()

        # ANALYTICS BEGIN
        self.project_parameters.AddEmptyValue("perform_analytics_option").SetBool(False)

        if self.project_parameters["perform_analytics_option"].GetBool():
            import KratosMultiphysics.SwimmingDEMApplication.analytics as analytics
            variables_to_measure = [Kratos.PRESSURE]
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

        self.FillHistoryForcePrecalculatedVectors()

        self.PerformZeroStepInitializations()

        if self.do_print_results:
            self._Print()

    def AddExtraProcessInfoVariablesToFluid(self):
        self.vars_man.AddExtraProcessInfoVariablesToFluidModelPart(self.project_parameters, self.fluid_model_part)

    def FluidInitialize(self):
        self.fluid_model_part = self._GetFluidAnalysis().fluid_model_part
        self._GetFluidAnalysis().vars_man = self.vars_man
        self._GetFluidAnalysis().Initialize()

        self.AddExtraProcessInfoVariablesToFluid()

        SDP.AddExtraDofs(self.fluid_model_part,
                         self.spheres_model_part,
                         self.cluster_model_part,
                         self.dem_inlet_model_part,
                         self.vars_man)

    def DispersePhaseInitialize(self):
        self.vars_man.__class__.AddNodalVariables(self.spheres_model_part, self.vars_man.dem_vars)
        self.vars_man.__class__.AddNodalVariables(self.rigid_face_model_part, self.vars_man.rigid_faces_vars)
        self.vars_man.__class__.AddNodalVariables(self.dem_inlet_model_part, self.vars_man.inlet_vars)
        self.vars_man.AddExtraProcessInfoVariablesToDispersePhaseModelPart(self.project_parameters,
                                                                           self._GetDEMAnalysis().spheres_model_part)

        self._GetDEMAnalysis().Initialize()

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
        if self.project_parameters["custom_fluid"]["embedded_option"].GetBool():
            self.calculate_distance_process = Kratos.CalculateSignedDistanceTo3DSkinProcess(
                self.rigid_face_model_part,
                self.fluid_model_part
                )
            self.calculate_distance_process.Execute()

    def GetAnalyticFacesModelParts(self):
        analytic_face_submodelpart_number = 1
        analytic_face_submodelpart_name = self.rigid_face_model_part.GetSubModelPart(str(analytic_face_submodelpart_number))
        return analytic_face_submodelpart_name

    def MakeAnalyticsMeasurements(self):
        self.analytic_face_watcher.MakeMeasurements()
        self.analytic_particle_watcher.MakeMeasurements()

    def InitializeSolutionStep(self):
        self.TellTime()
        self.PerformInitialDEMStepOperations(self.time)
        self._GetDEMAnalysis().InitializeSolutionStep()
        if self._GetSolver().CannotIgnoreFluidNow():
            self._GetFluidAnalysis().InitializeSolutionStep()
        super(SwimmingDEMAnalysis, self).InitializeSolutionStep()

    def FinalizeSolutionStep(self):
        # printing if required
        if self._GetSolver().CannotIgnoreFluidNow():
            self._GetFluidAnalysis().FinalizeSolutionStep()

        self._GetDEMAnalysis().FinalizeSolutionStep()

        # applying DEM-to-fluid coupling

        if self.DEM_to_fluid_counter.Tick() and self.time >= self.project_parameters["coupling"]["interaction_start_time"].GetDouble():
            self._GetSolver().projection_module.ProjectFromParticles()

        # coupling checks (debugging)
        if self.debug_info_counter.Tick():
            self.dem_volume_tool.UpdateDataAndPrint(
                self.project_parameters["fluid_domain_volume"].GetDouble())

        super(SwimmingDEMAnalysis, self).FinalizeSolutionStep()

    def OutputSolutionStep(self):
        # printing if required

        if self.print_counter.Tick():
            self.ComputePostProcessResults()
            self._Print()

        super(SwimmingDEMAnalysis, self).OutputSolutionStep()

    def _Print(self):
        os.chdir(self.post_path)
        #TODO: review this lines
        # import define_output
        # self.drag_list = define_output.DefineDragList()
        self.drag_file_output_list = []

        if self.particles_results_counter.Tick():
            self.io_tools.PrintParticlesResults(
                self.vars_man.variables_to_print_in_file,
                self.time,
                self.spheres_model_part)

        self.post_utils.Writeresults(self.time)
        os.chdir(self.main_path)

    def ComputePostProcessResults(self):
        if self.project_parameters["coupling"]["coupling_level_type"].GetInt():
            self._GetSolver().projection_module.ComputePostProcessResults(self.spheres_model_part.ProcessInfo)

    def GetFirstStepForFluidComputation(self):
        return 3

    def CloneTimeStep(self):
        self.fluid_model_part.CloneTimeStep(self.time)

    def DEMSolve(self, time='None'): # time is passed in case it is needed
        self._GetDEMAnalysis().solver.SolveSolutionStep()

    def UpdateALEMeshMovement(self, time):
        pass

    def RecoverDerivatives(self):
        self.recovery.Recover()

    def FluidSolve(self, time='None', solve_system=True):
        Say('Solving Fluid... (', self.fluid_model_part.NumberOfElements(0), 'elements )\n')

        if solve_system:
            self._GetFluidAnalysis().RunSingleTimeStep()
        else:
            Say("Skipping solving system...\n")

    def PerformZeroStepInitializations(self):
        pass

    def PerformInitialDEMStepOperations(self, time=None):
        pass

    def PerformEmbeddedOperations(self):
        # calculating elemental distances defining the structure embedded in the fluid mesh
        if self.project_parameters["custom_fluid"]["embedded_option"].GetBool():
            self.calculate_distance_process.Execute()

        if self.embedded_counter.Tick():
            embedded.ApplyEmbeddedBCsToFluid(self.fluid_model_part)
            embedded.ApplyEmbeddedBCsToBalls(self.spheres_model_part, self.project_parameters)

    def AssessStationarity(self):
        self.stationarity = self._GetSolver().AssessStationarity()
        self.stationarity_counter.Deactivate(self.stationarity)

    def SetInlet(self):
        if self.project_parameters["dem_inlet_option"].GetBool():
            # Constructing the inlet and initializing it
            # (must be done AFTER the self.spheres_model_part Initialize)
            # Note that right now only inlets of a single type are possible.
            # This should be generalized.
            if self.project_parameters["type_of_inlet"].GetString() == 'VelocityImposed':
                self.DEM_inlet = DEM_Inlet(self.dem_inlet_model_part)
            elif self.project_parameters["type_of_inlet"].GetString() == 'ForceImposed':
                self.DEM_inlet = DEM_Force_Based_Inlet(self.dem_inlet_model_part, self.project_parameters["inlet_force_vector"].GetVector())

            self._GetDEMAnalysis().DEM_inlet = self.DEM_inlet
            self.DEM_inlet.InitializeDEM_Inlet(self.spheres_model_part, self._GetDEMAnalysis().creator_destructor)

    def SetAnalyticParticleWatcher(self):
        from analytic_tools import analytic_data_procedures
        self.particle_watcher = DEM.AnalyticParticleWatcher()
        self.particle_watcher_analyser = analytic_data_procedures.ParticleWatcherAnalyzer(
            analytic_particle_watcher=self.particle_watcher,
            path=self.main_path)

    def ProcessAnalyticData(self):
        self._GetDEMAnalysis().WriteAnalyticDataToFileAndClear()

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

    def GetBackwardCouplingCounter(self):
        return SDP.Counter(1, 1, self.project_parameters["coupling"]["coupling_level_type"].GetInt() > 1)

    def GetRecoveryCounter(self):
        there_is_something_to_recover = (
            self.project_parameters["coupling"]["coupling_level_type"].GetInt() or
            self.project_parameters["print_PRESSURE_GRADIENT_option"].GetBool())
        return SDP.Counter(1, 1, there_is_something_to_recover)

    def GetStationarityCounter(self):
        return SDP.Counter(
            steps_in_cycle=self.project_parameters["stationarity"]["time_steps_per_stationarity_step"].GetInt(),
            beginning_step=1,
            is_active=self.project_parameters["stationarity"]["stationary_problem_option"].GetBool())

    def GetPrintCounter(self):
        counter = SDP.Counter(steps_in_cycle=int(self.output_time / self.time_step + 0.5),
                              beginning_step=int(self.output_time / self.time_step),
                              is_dead = not self.do_print_results)
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
            steps_in_cycle=self.project_parameters["stationarity"]["time_steps_per_analytic_processing_step"].GetInt(),
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

    def FillHistoryForcePrecalculatedVectors(self): # TODO: more robust implementation
        # Warning: this estimation is based on a constant time step for DEM.
        # This is usually the case, but could not be so.
        for prop in self.project_parameters["properties"].values():
            if prop["hydrodynamic_law_parameters"].Has("history_force_parameters"):

                if prop["hydrodynamic_law_parameters"]["history_force_parameters"]["name"].GetString() != 'default':
                    total_number_of_steps = int(self.end_time / self.project_parameters["time_stepping"]["time_step"].GetDouble()) + 20
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

        self._GetFluidAnalysis().Finalize()

        self.TellFinalSummary(self.time, self.step, self._GetSolver().fluid_step)

    def SetPointGraphPrinter(self):
        pass

    def AssignKinematicViscosityFromDynamicViscosity(self):
        # Eulerian fluid already works with kinematic viscosity
        pass

    def GetReturnValue(self):
        return 0.0

    def _GetDEMAnalysis(self):
        if not hasattr(self, '_disperse_phase_analysis'):
            import KratosMultiphysics.SwimmingDEMApplication.fluid_coupled_DEM_analysis as DEM_analysis
            self._disperse_phase_analysis = DEM_analysis.FluidCoupledDEMAnalysisStage(self.model, self.project_parameters)

        return self._disperse_phase_analysis

    def _GetFluidAnalysis(self):
        if not hasattr(self, '_fluid_phase_analysis'):
            import KratosMultiphysics.SwimmingDEMApplication.DEM_coupled_fluid_dynamics_analysis as fluid_analysis
            self._fluid_phase_analysis = fluid_analysis.DEMCoupledFluidDynamicsAnalysis(self.model, self.project_parameters, self.vars_man)
            self._fluid_phase_analysis.main_path = self.main_path
        return self._fluid_phase_analysis

    # To-do: for the moment, provided for compatibility
    def _CreateSolver(self):
        import KratosMultiphysics.SwimmingDEMApplication.swimming_DEM_solver as swimming_DEM_solver
        return swimming_DEM_solver.SwimmingDEMSolver(self.model,
                                                     self.project_parameters,
                                                     self.GetFieldUtility(),
                                                     self._GetFluidAnalysis()._GetSolver(),
                                                     self._GetDEMAnalysis()._GetSolver(),
                                                     self.vars_man)

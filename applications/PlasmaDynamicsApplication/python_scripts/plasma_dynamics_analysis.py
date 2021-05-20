import os
import sys
import math
import time as timer
import weakref

import KratosMultiphysics as Kratos
from KratosMultiphysics import Logger, Parameters

import KratosMultiphysics.DEMApplication as DEM
import KratosMultiphysics.DEMApplication.DEM_procedures as DP

from KratosMultiphysics.SwimmingDEMApplication.swimming_DEM_analysis import SwimmingDEMAnalysis
from KratosMultiphysics.SwimmingDEMApplication.swimming_DEM_analysis import SDEMLogger

import KratosMultiphysics.PlasmaDynamicsApplication.plasma_dynamics_procedures as PDP
import KratosMultiphysics.PlasmaDynamicsApplication.plasma_dynamics_gid_output as plasma_dynamics_gid_output
import KratosMultiphysics.PlasmaDynamicsApplication.variables_management as variables_management



def Say(*args):
    Logger.PrintInfo("PlasmaDynamics", *args)
    Logger.Flush()



class PlasmaDynamicsLogger(SDEMLogger):
    
    def __init__(self, do_print_file=False):
        super().__init__()



class python_parameters:
    
    def __init__(self):        
        pass



class PlasmaDynamicsAnalysis(SwimmingDEMAnalysis):
    
    def __enter__ (self):      
        return self


    def __exit__(self, exception_type, exception_value, traceback):
        pass


    def __init__(self, model, parameters = Parameters("{}")):
        sys.stdout = PlasmaDynamicsLogger()
        self.StartTimer()
        self.model = model
        self.main_path = os.getcwd()

        self.SetProjectParameters(parameters)

        self.vars_man = variables_management.VariablesManager(self.project_parameters)

        self._GetDEMAnalysis().coupling_analysis = weakref.proxy(self)

        self._GetFluidPhaseAnalysis().coupling_analysis = weakref.proxy(self)

        self.procedures = weakref.proxy(self._GetDEMAnalysis().procedures)

        self.report = DP.Report()

        # self._GetDEMAnalysis().SetAnalyticFaceWatcher()

        self.fluid_model_part = self._GetFluidPhaseAnalysis().fluid_model_part
        self.thermal_model_part = self._GetFluidPhaseAnalysis().thermal_model_part
        self.spheres_model_part = self._GetDEMAnalysis().spheres_model_part
        self.cluster_model_part = self._GetDEMAnalysis().cluster_model_part
        self.rigid_face_model_part = self._GetDEMAnalysis().rigid_face_model_part
        self.dem_inlet_model_part = self._GetDEMAnalysis().dem_inlet_model_part
        self.vars_man.ConstructListsOfVariables(self.project_parameters)

        super(SwimmingDEMAnalysis, self).__init__(model, self.project_parameters)
        
        # super().__init__(model, parameters)
        # self.thermic_model_part = self._GetFluidPhaseAnalysis().thermic_model_part


    def SetProjectParameters(self, parameters):
        self.project_parameters = parameters
        self.time_step = self.project_parameters["time_stepping"]["time_step"].GetDouble()
        self.end_time   = self.project_parameters["problem_data"]["end_time"].GetDouble()
        self.do_print_results = self.project_parameters["do_print_results_option"].GetBool()
        self.fluid_phase_parameters = self.project_parameters['fluid_phase_parameters'] 
        self.fluid_solver_settings = self.fluid_phase_parameters["solver_settings"]["fluid_solver_settings"]
        self.thermal_solver_settings = self.fluid_phase_parameters["solver_settings"]["thermal_solver_settings"]
        
        
        # self.project_parameters["fluid_parameters"] = {
        #     "problem_data"        : self.fluid_phase_parameters["problem_data"],
        #     "output_processes"    : self.fluid_phase_parameters["output_processes"],               
        #     "processes"           : self.fluid_phase_parameters["processes"],                
        #     "solver_settings"     : self.fluid_solver_settings
        #                                       }
        
        # self.fluid_parameters = self.project_parameters["fluid_parameters"]
        
        
        self.fluid_parameters = {
            "problem_data"        : self.fluid_phase_parameters["problem_data"],
            "output_processes"    : self.fluid_phase_parameters["output_processes"],               
            "processes"           : self.fluid_phase_parameters["processes"],                
            "solver_settings"     : self.fluid_solver_settings
                                  }
        # self.project_parameters["fluid_parameters"] = self.fluid_parameters

        import KratosMultiphysics.PlasmaDynamicsApplication.plasma_dynamics_default_input_parameters as only_plasma_defaults
        import KratosMultiphysics.DEMApplication.dem_default_input_parameters as dem_defaults

        self.project_parameters.ValidateAndAssignDefaults(only_plasma_defaults.GetDefaultInputParameters())
        self.project_parameters["dem_parameters"].ValidateAndAssignDefaults(dem_defaults.GetDefaultInputParameters())
   
        self.SetBetaParameters()
        
        self.ModifyInputParametersForCoherence()
        
        
        
        
           
    # TODO
    # Set input parameters that have not yet been transferred to the interface
    # import the configuration data as read from the GiD
    # def SetBetaParameters(self):
    #     pass


    # TODO
    # This step is added to allow modifications to the possibly incompatibilities
    # between the individual parameters coming from each sub-application
    # (i.e., fluid and dem apps)
    def ModifyInputParametersForCoherence(self):
        # Making all time steps exactly commensurable
        output_time = self.project_parameters["output_interval"].GetDouble()
        self.output_time = int(output_time / self.time_step) * self.time_step
        self.project_parameters["output_interval"].SetDouble(self.output_time)
        
        self.fluid_time_step = self.fluid_solver_settings["time_stepping"]["time_step"].GetDouble()

        if self.fluid_time_step < self.time_step:
            error_message = ('The fluid time step (' + str(self.fluid_time_step)
                             + ') must be larger or equal than the overall time step (' + str(self.time_step)
                             + ')!')
            raise Exception(error_message)

        self.fluid_time_step = int(self.fluid_time_step / self.time_step) * self.time_step
        self.fluid_solver_settings["time_stepping"]["time_step"].SetDouble(self.fluid_time_step)
        
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


    def InitializeVariablesWithNonZeroValues(self):       
        PDP.InitializeVariablesWithNonZeroValues(self.project_parameters,
                                                self.fluid_model_part,
                                                self.spheres_model_part)


    # TODO: Create a custom initialization for the application
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
            PDP.CopyInputFilesIntoFolder(self.main_path, self.post_path)
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

            self.plasma_dynamics_gid_io = \
            plasma_dynamics_gid_output.PlasmaDynamicsGiDOutput(
                file_name = self.project_parameters["problem_data"]["problem_name"].GetString(),
                vol_output = result_file_configuration["body_output"].GetBool(),
                post_mode = old_gid_output_post_options_dict[post_mode_key],
                multifile = old_gid_output_multiple_file_option_dict[multiple_files_option_key],
                deformed_mesh = deformed_mesh_option,
                write_conditions = write_conditions_option)

            self.plasma_dynamics_gid_io.initialize_plasma_dynamics_results(
                self.spheres_model_part,
                self.cluster_model_part,
                self.rigid_face_model_part,
                self.mixed_model_part)

        self.SetPointGraphPrinter()

        self.AssignKinematicViscosityFromDynamicViscosity()

        super(SwimmingDEMAnalysis, self).Initialize()

        # coarse-graining: applying changes to the physical properties of the model to adjust for
        # the similarity transformation if required (fluid effects only).
        PDP.ApplySimilarityTransformations(
            self.fluid_model_part,
            self.project_parameters["similarity"]["similarity_transformation_type"].GetInt(),
            self.project_parameters["similarity"]["model_over_real_diameter_factor"].GetDouble()
            )

        if self.do_print_results:
            self.SetPostUtils()

        # creating an IOTools object to perform other printing tasks
        self.io_tools = PDP.IOTools(self.project_parameters)

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
        fluid_domain_dimension = self.project_parameters["fluid_phase_parameters"]["solver_settings"]["domain_size"].GetInt()
        self.custom_functions_tool = PDP.FunctionsCalculator(fluid_domain_dimension)

        # creating a debug tool
        self.dem_volume_tool = self.GetVolumeDebugTool()

        #self.SetEmbeddedTools()

        Say('Initialization Complete\n')

        if self.project_parameters["custom_fluid"]["flow_in_porous_DEM_medium_option"].GetBool():
            PDP.FixModelPart(self.spheres_model_part)


        ##################################################
        #  I N I T I A L I Z I N G    T I M E    L O O P
        ##################################################
        self.step = 0
        self.time = self.fluid_parameters["problem_data"]["start_time"].GetDouble()
        self.fluid_time_step = self._GetFluidPhaseAnalysis()._GetSolver().fluid_solver._ComputeDeltaTime()
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
        self.mat_deriv_averager           = PDP.Averager(1, 3)
        self.laplacian_averager           = PDP.Averager(1, 3)

        self.report.total_steps_expected = int(self.end_time / self.time_step)

        Say(self.report.BeginReport(self.timer))

        # creating a Post Utils object that executes several post-related tasks
        self.post_utils_DEM = DP.PostUtils(self.project_parameters['dem_parameters'], self.spheres_model_part)

        # otherwise variables are set to 0 by default:
        self.InitializeVariablesWithNonZeroValues()

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
            target_node = PDP.FindClosestNode(self.fluid_model_part, point_coors)
            target_id = target_node.Id
            Say(target_node.X, target_node.Y, target_node.Z)
            Say(target_id)
            def condition(node):
                return node.Id == target_id

            gauge.ConstructArrayOfNodes(condition)
            Say(gauge.variables)
        # ANALYTICS END

        self.FillHistoryForcePrecalculatedVectors()

        import KratosMultiphysics.SwimmingDEMApplication.derivative_recovery.derivative_recovery_strategy as derivative_recoverer
        self.recovery = derivative_recoverer.DerivativeRecoveryStrategy(
            self.project_parameters,
            self.fluid_model_part,
            self.custom_functions_tool)

        self.PerformZeroStepInitializations()

        if self.do_print_results:
            self._Print()


    def FluidInitialize(self):
        self.fluid_model_part = self._GetFluidPhaseAnalysis().fluid_model_part
        self._GetFluidPhaseAnalysis().vars_man = self.vars_man
        self._GetFluidPhaseAnalysis().Initialize()

        self.AddExtraProcessInfoVariablesToFluid()

        PDP.AddExtraDofs(self.fluid_model_part,
                         self.spheres_model_part,
                         self.cluster_model_part,
                         self.dem_inlet_model_part,
                         self.vars_man)


    def SetPostUtils(self):
        # creating a Post Utils object that executes several post-related tasks
        self.post_utils = PDP.PostUtils(self.plasma_dynamics_gid_io,
                                        self.project_parameters,
                                        self.vars_man,
                                        self.fluid_model_part,
                                        self.spheres_model_part,
                                        self.cluster_model_part,
                                        self.rigid_face_model_part,
                                        self.mixed_model_part)


    def GetBackwardCouplingCounter(self):
        return PDP.Counter(1, 1, self.project_parameters["coupling"]["coupling_level_type"].GetInt() > 1)


    def GetRecoveryCounter(self):
        there_is_something_to_recover = (
            self.project_parameters["coupling"]["coupling_level_type"].GetInt() or
            self.project_parameters["print_PRESSURE_GRADIENT_option"].GetBool())
        return PDP.Counter(1, 1, there_is_something_to_recover)


    def GetStationarityCounter(self):
        return PDP.Counter(
            steps_in_cycle=self.project_parameters["stationarity"]["time_steps_per_stationarity_step"].GetInt(),
            beginning_step=1,
            is_active=self.project_parameters["stationarity"]["stationary_problem_option"].GetBool())


    def GetPrintCounter(self):
        counter = PDP.Counter(steps_in_cycle=int(self.output_time / self.time_step + 0.5),
                              beginning_step=int(self.output_time / self.time_step),
                              is_dead = not self.do_print_results)
        return counter


    def GetDebugInfo(self):
        return PDP.Counter(
            self.project_parameters["debug_tool_cycle"].GetInt(),
            1,
            self.project_parameters["print_debug_info_option"].GetBool())


    def GetParticlesResultsCounter(self):
        return PDP.Counter(
            self.project_parameters["print_particles_results_cycle"].GetInt(),
            1,
            self.project_parameters["print_particles_results_option"].GetBool())


    def GetHistoryForceQuadratureCounter(self):
        for prop in self.project_parameters["properties"].values():
            if prop["hydrodynamic_law_parameters"].Has("history_force_parameters"):
                history_force_parameters =  prop["hydrodynamic_law_parameters"]["history_force_parameters"]
                if history_force_parameters.Has("time_steps_per_quadrature_step"):
                    time_steps_per_quadrature_step = history_force_parameters["time_steps_per_quadrature_step"].GetInt()
                    return PDP.Counter(steps_in_cycle=time_steps_per_quadrature_step, beginning_step=1)
        return PDP.Counter(is_dead=True)


    def ProcessAnalyticDataCounter(self):
        return PDP.Counter(
            steps_in_cycle=self.project_parameters["stationarity"]["time_steps_per_analytic_processing_step"].GetInt(),
            beginning_step=1,
            is_active=self.project_parameters["do_process_analytic_data"].GetBool())


    def GetVolumeDebugTool(self):
        return PDP.ProjectionDebugUtils(
            self.project_parameters["fluid_domain_volume"].GetDouble(),
            self.fluid_model_part,
            self.spheres_model_part,
            self.custom_functions_tool)


    def Finalize(self):
        Say('Finalizing simulation...\n')
        
        if self.do_print_results:
            self.plasma_dynamics_gid_io.finalize_results()

        self.PerformFinalOperations(self.time)

        self._GetFluidPhaseAnalysis().Finalize()

        self.TellFinalSummary(self.time, self.step, self._GetSolver().fluid_step)


    def _GetDEMAnalysis(self):
        if not hasattr(self, '_disperse_phase_analysis'):
            import KratosMultiphysics.PlasmaDynamicsApplication.fluid_coupled_DEM_analysis as DEM_analysis
            self._disperse_phase_analysis = DEM_analysis.FluidCoupledDEMAnalysis(self.model, self.project_parameters)
        return self._disperse_phase_analysis


    def _GetFluidPhaseAnalysis(self):
        if not hasattr(self, '_fluid_phase_analysis'):
            import KratosMultiphysics.PlasmaDynamicsApplication.DEM_coupled_fluid_analysis as fluid_analysis
            self._fluid_phase_analysis = fluid_analysis.DEMCoupledFluidAnalysis(self.model, self.project_parameters, self.vars_man)
            self._fluid_phase_analysis.main_path = self.main_path
        return self._fluid_phase_analysis


    def _CreateSolver(self):
        import KratosMultiphysics.PlasmaDynamicsApplication.plasma_dynamics_solver as plasma_dynamics_solver
        return plasma_dynamics_solver.PlasmaDynamicsSolver(self.model,
                                                     self.project_parameters,
                                                     self.GetFieldUtility(),
                                                     self._GetFluidPhaseAnalysis()._GetSolver(),
                                                     self._GetDEMAnalysis()._GetSolver(),
                                                     self.vars_man)

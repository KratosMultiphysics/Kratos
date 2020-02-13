from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

from KratosMultiphysics import Logger, Parameters
from KratosMultiphysics.python_solver import PythonSolver
import KratosMultiphysics.SwimmingDEMApplication as SDEM
import math
import KratosMultiphysics.SwimmingDEMApplication.swimming_DEM_procedures as SDP
import KratosMultiphysics.SwimmingDEMApplication.parameters_tools as PT
import KratosMultiphysics.SwimmingDEMApplication.CFD_DEM_coupling as CFD_DEM_coupling
import KratosMultiphysics.SwimmingDEMApplication.derivative_recovery.derivative_recovery_strategy as derivative_recoverer

def Say(*args):
    Logger.PrintInfo("SwimmingDEM", *args)
    Logger.Flush()

class SwimmingDEMSolver(PythonSolver):
    def _ValidateSettings(self, project_parameters):

        default_processes_settings = Parameters("""{
                "python_module" : "calculate_nodal_area_process",
                "kratos_module" : "KratosMultiphysics.SwimmingDEMApplication",
                "process_name"  : "CalculateNodalAreaProcess",
                "Parameters"    : {
                    "model_part_name" : "FluidModelPart",
                    "domain_size" : 3,
                    "fixed_mesh": false
                }
            }

        """)

        if not project_parameters["processes"].Has('non_optional_solver_processes'):
            project_parameters["processes"].AddEmptyArray("non_optional_solver_processes")

        else: # reconstruct non_optional_solver_processes list making sure calculate_nodal_area_process is not added twice
            non_optional_processes_list = list(project_parameters["processes"]["non_optional_solver_processes"])
            project_parameters["processes"].Remove("non_optional_solver_processes")
            project_parameters["processes"].AddEmptyArray("non_optional_solver_processes")

            for process in non_optional_processes_list:
                if process["python_module"].GetString() != 'calculate_nodal_area_process':
                    project_parameters["processes"]["non_optional_solver_processes"].Append(process)

        non_optional_solver_processes = project_parameters["processes"]["non_optional_solver_processes"]
        non_optional_solver_processes.Append(default_processes_settings)
        nodal_area_process_parameters = non_optional_solver_processes[non_optional_solver_processes.size() -1]["Parameters"]
        nodal_area_process_parameters["model_part_name"].SetString(self.fluid_solver.main_model_part.Name)
        nodal_area_process_parameters["domain_size"].SetInt(self.fluid_domain_dimension)
        the_mesh_moves = False
        if self.fluid_solver.settings.Has('move_mesh_flag'):
            the_mesh_moves = self.fluid_solver.settings["move_mesh_flag"].GetBool()
            nodal_area_process_parameters["fixed_mesh"].SetBool(not the_mesh_moves)
        elif self.fluid_solver.settings.Has('time_integration_settings'):
            the_mesh_moves = self.fluid_solver.settings["time_integration_settings"]["move_mesh_flag"].GetBool()
            nodal_area_process_parameters["fixed_mesh"].SetBool(not the_mesh_moves)
        elif self.fluid_solver.settings["solvers"][0]["Parameters"]["time_integration_settings"].Has('move_mesh_flag'):
            the_mesh_moves = self.fluid_solver.settings["solvers"][0]["Parameters"]["time_integration_settings"]["move_mesh_flag"].GetBool()
            nodal_area_process_parameters["fixed_mesh"].SetBool(not the_mesh_moves)
        self.move_mesh_flag = the_mesh_moves
        return project_parameters

    def __init__(self, model, project_parameters, field_utility, fluid_solver, dem_solver, variables_manager):
        # Validate settings
        self.field_utility = field_utility
        self.vars_man = variables_manager
        self.fluid_domain_dimension = project_parameters["fluid_parameters"]["solver_settings"]["domain_size"].GetInt()
        self.fluid_solver = fluid_solver
        self.dem_solver = dem_solver
        self.project_parameters = self._ValidateSettings(project_parameters)
        self.next_time_to_solve_fluid = project_parameters['problem_data']['start_time'].GetDouble()
        self.coupling_level_type = project_parameters["coupling"]["coupling_level_type"].GetInt()
        self.interaction_start_time = project_parameters["coupling"]["interaction_start_time"].GetDouble()
        self.integration_scheme = project_parameters["custom_dem"]["translational_integration_scheme"].GetString()
        self.fluid_dt = fluid_solver.settings["time_stepping"]["time_step"].GetDouble()
        self.do_solve_dem = project_parameters["custom_dem"]["do_solve_dem"].GetBool()
        self.solve_system = not self.project_parameters["custom_fluid"]["fluid_already_calculated"].GetBool()

        self.fluid_step = 0
        self.calculating_fluid_in_current_step = True
        self.first_DEM_iteration = True
        self.SetHistoryForceOptions()
        self.ConstructStationarityTool()
        self.ConstructDerivativeRecoverer()
        self.ConstructHistoryForceUtility()
        # Call the base Python solver constructor
        super(SwimmingDEMSolver, self).__init__(model, project_parameters)

    def ConstructStationarityTool(self):
        self.stationarity = False
        self.stationarity_counter = self.GetStationarityCounter()
        self.stationarity_tool = SDEM.FlowStationarityCheck(self.fluid_solver.main_model_part,
            self.project_parameters["stationarity"]["tolerance"].GetDouble()
            )

    def _ConstructProjectionModule(self):
        # creating a projection module for the fluid-DEM coupling
        self.h_min = 0.01 #TODO: this must be set from interface and the method must be checked for 2D
        n_balls = 1
        fluid_volume = 10
        # the variable n_particles_in_depth is only relevant in 2D problems
        self.project_parameters.AddEmptyValue("n_particles_in_depth").SetInt(int(math.sqrt(n_balls / fluid_volume)))

        projection_module = CFD_DEM_coupling.ProjectionModule(
        self.fluid_solver.main_model_part,
        self.dem_solver.spheres_model_part,
        self.dem_solver.all_model_parts.Get("RigidFacePart"),
        self.project_parameters,
        self.vars_man.coupling_dem_vars,
        self.vars_man.coupling_fluid_vars,
        self.vars_man.time_filtered_vars,
        flow_field=self.field_utility,
        domain_size=self.fluid_domain_dimension
        )

        projection_module.UpdateDatabase(self.h_min)

        return projection_module

    def SetHistoryForceOptions(self):
        self.history_force_on = False
        self.MAE_parameters = Parameters("{}")
        for prop in self.project_parameters["properties"].values(): #TODO: now it only works for one property!
            self.history_force_on = (PT.RecursiveFindParametersWithCondition(
                                     self.project_parameters["properties"], 'history_force_parameters',
                                     condition=lambda value: value['name'].GetString() != 'default'))
            if self.history_force_on:
                self.MAE_parameters = prop["hydrodynamic_law_parameters"]["history_force_parameters"]["mae_parameters"]
            break
        self.do_use_mae = PT.RecursiveFindTrueBoolInParameters(self.MAE_parameters, 'do_use_mae')


    def ConstructDerivativeRecoverer(self):
        self.derivative_recovery_counter = self.GetRecoveryCounter()

        self.recovery = derivative_recoverer.DerivativeRecoveryStrategy(
            self.project_parameters,
            self.fluid_solver.main_model_part,
            SDP.FunctionsCalculator(self.fluid_domain_dimension))

    def ConstructHistoryForceUtility(self):
        self.quadrature_counter = self.GetHistoryForceQuadratureCounter()
        if self.history_force_on:
            self.basset_force_tool = SDEM.BassetForceTools(self.MAE_parameters)

    def GetStationarityCounter(self):
        return SDP.Counter(
            steps_in_cycle=self.project_parameters["stationarity"]["time_steps_per_stationarity_step"].GetInt(),
            beginning_step=self.project_parameters["stationarity"]["time_steps_before_first_assessment"].GetInt(),
            is_active=self.project_parameters["stationarity"]["stationary_problem_option"].GetBool())

    def GetRecoveryCounter(self):
        there_is_something_to_recover = (
            self.project_parameters["coupling"]["coupling_level_type"].GetInt() or
            self.project_parameters["print_PRESSURE_GRADIENT_option"].GetBool())
        return SDP.Counter(1, 1, there_is_something_to_recover)

    def GetHistoryForceQuadratureCounter(self):
        for prop in self.project_parameters["properties"].values():
            if prop["hydrodynamic_law_parameters"].Has("history_force_parameters"):
                history_force_parameters =  prop["hydrodynamic_law_parameters"]["history_force_parameters"]
                if history_force_parameters.Has("time_steps_per_quadrature_step"):
                    time_steps_per_quadrature_step = history_force_parameters["time_steps_per_quadrature_step"].GetInt()

                    return SDP.Counter(steps_in_cycle=time_steps_per_quadrature_step, beginning_step=1)

        return SDP.Counter(is_dead=True)

    def AdvanceInTime(self, time):
        self.time = self.dem_solver.AdvanceInTime(time)
        self.calculating_fluid_in_current_step = bool(time >= self.next_time_to_solve_fluid - 0.5 * self.dem_solver.dt)

        if self.calculating_fluid_in_current_step:
            self.next_time_to_solve_fluid = self.fluid_solver.AdvanceInTime(time)
            self.fluid_step += 1

        return self.time

    def UpdateALEMeshMovement(self, time): # TODO: move to derived solver
        if self.project_parameters["custom_fluid"]["ALE_option"].GetBool():
            self.rotator.RotateMesh(self.fluid_solver.main_model_part, time)
            self._GetProjectionModule().UpdateDatabase(self.CalculateMinElementSize())

    def CalculateMinElementSize(self):
        return self.h_min

    def AssessStationarity(self):
        Say("Assessing Stationarity...\n")
        self.stationarity = self.stationarity_tool.AssessStationarity()
        if not self.stationarity:
            tolerance = self.stationarity_tool.GetTolerance()
            p_dot = self.stationarity_tool.GetCurrentPressureDerivative()
            p_dot_historical = self.stationarity_tool.GetCharacteristicPressureDerivative()
            non_stationarity_measure = self.stationarity_tool.GetTransienceMeasure()

            message = '\nFluid not stationary:\n'
            message += '  * Current average pressure time derivative: ' + str(p_dot) + '\n'
            message += '  * Historic average: ' + str(p_dot_historical) + '\n'
            message += '  * Current transience measure: ' + str(non_stationarity_measure) + ' > ' + str(tolerance) + '\n'
            Say(message)
        self.stationarity_counter.Deactivate(self.stationarity)
        return self.stationarity

    # Compute nodal quantities to be printed that are not generated as part of the
    # solution algorithm. For instance, the pressure gradient, which is not used for
    # the coupling but can be of interest.
    def ComputePostProcessResults(self):
        if self.project_parameters["coupling"]["coupling_level_type"].GetInt():
            self._GetProjectionModule().ComputePostProcessResults(self.dem_solver.spheres_model_part.ProcessInfo)

    def CannotIgnoreFluidNow(self):
        return self.solve_system and self.calculating_fluid_in_current_step

    def Predict(self):
        if self.CannotIgnoreFluidNow():
            self.fluid_solver.Predict()
        self.dem_solver.Predict()

    def ApplyForwardCoupling(self, alpha='None'):
        self._GetProjectionModule().ApplyForwardCoupling(alpha)

    def ApplyForwardCouplingOfVelocityToAuxVelocityOnly(self, alpha=None):
        self._GetProjectionModule().ApplyForwardCouplingOfVelocityToAuxVelocityOnly(alpha)

    def _GetProjectionModule(self):
        if not hasattr(self, 'projection_module'):
            self.projection_module = self._ConstructProjectionModule()
        return self.projection_module

    def SolveSolutionStep(self):
        # update possible movements of the fluid mesh
        self.UpdateALEMeshMovement(self.time)

        # Solving the fluid part
        Say('Solving Fluid... (', self.fluid_solver.main_model_part.NumberOfElements(0), 'elements )\n')
        self.solve_system = not self.project_parameters["custom_fluid"]["fluid_already_calculated"].GetBool() and not self.stationarity

        if self.CannotIgnoreFluidNow():
            self.SolveFluidSolutionStep()
        else:
            Say("Skipping solving system for the fluid phase...\n")

        self.recovery = derivative_recoverer.DerivativeRecoveryStrategy(
            self.project_parameters,
            self.fluid_solver.computing_model_part,
            SDP.FunctionsCalculator(self.fluid_domain_dimension))

        self.derivative_recovery_counter.Activate(self.time > self.interaction_start_time and self.calculating_fluid_in_current_step)

        if self.derivative_recovery_counter.Tick():
            self.recovery.Recover()

        # Solving the disperse-phase component
        Say('Solving DEM... (', self.dem_solver.spheres_model_part.NumberOfElements(0), 'elements )')
        self.SolveDEM()

        return True

    def SolveFluidSolutionStep(self):
        self.fluid_solver.SolveSolutionStep()
        if self.move_mesh_flag:
            self._GetProjectionModule().UpdateDatabase(self.CalculateMinElementSize())
        else: # stationarity can only checked for fixed meshes for the moment
            # Check for stationarity: this is useful for steady-state problems, so that
            # the calculation stops after reaching the solution.
            if self.stationarity_counter.Tick():
                self.AssessStationarity()

    def SolveDEMSolutionStep(self):
        self.dem_solver.SolveSolutionStep()

    def SolveDEM(self):
        #self.PerformEmbeddedOperations() TO-DO: it's crashing

        it_is_time_to_forward_couple = (self.time >= self.interaction_start_time
                                        and self.coupling_level_type)

        alpha = 1.0 - (self.next_time_to_solve_fluid - self.time) / self.fluid_dt

        if (not self.move_mesh_flag
            and (it_is_time_to_forward_couple or self.first_DEM_iteration)):
                self.ApplyForwardCoupling(alpha)

        if self.quadrature_counter.Tick():
            self.AppendValuesForTheHistoryForce()

        if self.integration_scheme in {'Hybrid_Bashforth', 'TerminalVelocityScheme'}:
            # Advance in space only
            if self.do_solve_dem:
                self.SolveDEMSolutionStep()
            if it_is_time_to_forward_couple or self.first_DEM_iteration:
                self.ApplyForwardCouplingOfVelocityToAuxVelocityOnly(alpha)

        # Performing the time integration of the DEM part
        if (self.move_mesh_flag
            and (it_is_time_to_forward_couple or self.first_DEM_iteration)):
            self.ApplyForwardCoupling(alpha)

        if self.do_solve_dem:
            self.SolveDEMSolutionStep()

        self.first_DEM_iteration = False

    def AppendValuesForTheHistoryForce(self):
        if PT.RecursiveFindTrueBoolInParameters(self.MAE_parameters, 'do_use_mae'):
            self.basset_force_tool.AppendIntegrandsWindow(self.dem_solver.spheres_model_part)
        else:
            self.basset_force_tool.AppendIntegrands(self.dem_solver.spheres_model_part)

    def ImportModelPart(self): # TODO: implement this
        pass

    def GetComputingModelPart(self):
        return self.dem_solver.spheres_model_part

from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics
from python_solver import PythonSolver

# Import applications
import KratosMultiphysics.SwimmingDEMApplication as SwimmingDEMApplication
import swimming_DEM_procedures as SDP
import CFD_DEM_coupling
import derivative_recovery.derivative_recovery_strategy as derivative_recoverer
import math

def Say(*args):
    KratosMultiphysics.Logger.PrintInfo("SwimmingDEM", *args)
    KratosMultiphysics.Logger.Flush()

class SwimmingDEMSolver(PythonSolver):
    def _ValidateSettings(self, project_parameters):
        return project_parameters

    def __init__(self, model, project_parameters, fluid_solver, dem_solver, pp):
        # Validate settings
        self.pp = pp
        self.project_parameters = self._ValidateSettings(project_parameters)
        self.fluid_solver = fluid_solver
        self.dem_solver = dem_solver
        self.fluid_step = 0
        self.calculating_fluid_in_current_step = True
        self.next_time_to_solve_fluid = project_parameters['problem_data']['start_time'].GetDouble()
        self.coupling_level_type = project_parameters["coupling_level_type"].GetInt()
        self.coupling_scheme_type = project_parameters["coupling_scheme_type"].GetString()
        self.interaction_start_time = project_parameters["interaction_start_time"].GetDouble()
        self.project_at_every_substep_option = project_parameters["project_at_every_substep_option"].GetBool()
        self.integration_scheme = project_parameters["TranslationalIntegrationScheme"].GetString()
        self.fluid_dt = fluid_solver.settings["time_stepping"]["time_step"].GetDouble()
        self.do_solve_dem = project_parameters["do_solve_dem"].GetBool()
        self.solve_system = not self.project_parameters["fluid_already_calculated"].GetBool()
        self.first_DEM_iteration = True
        self.ConstructStationarityTool()
        self.ConstructDerivativeRecoverer()
        self.ConstructHistoryForceUtility()
        # Call the base Python solver constructor
        super(SwimmingDEMSolver, self).__init__(model, project_parameters)

    def ConstructStationarityTool(self):
        self.stationarity = False
        self.stationarity_counter = self.GetStationarityCounter()
        self.stationarity_tool = SDP.StationarityAssessmentTool(
            self.project_parameters["max_pressure_variation_rate_tol"].GetDouble(),
            SDP.FunctionsCalculator()
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
        self.pp.CFD_DEM,
        self.pp.coupling_dem_vars,
        self.pp.coupling_fluid_vars,
        self.pp.time_filtered_vars,
        flow_field=self.pp.field_utility,
        domain_size=self.pp.domain_size
        )

        projection_module.UpdateDatabase(self.h_min)

        return projection_module

    def ConstructDerivativeRecoverer(self):
        self.derivative_recovery_counter = self.GetRecoveryCounter()
        self.using_hinsberg_method = bool(self.project_parameters["basset_force_type"].GetInt() >= 3 or
                                          self.project_parameters["basset_force_type"].GetInt() == 1)
        self.recovery = derivative_recoverer.DerivativeRecoveryStrategy(
            self.pp,
            self.fluid_solver.main_model_part,
            SDP.FunctionsCalculator(self.pp.domain_size))

    def ConstructHistoryForceUtility(self):
        self.quadrature_counter = self.GetHistoryForceQuadratureCounter()
        self.basset_force_tool = SwimmingDEMApplication.BassetForceTools()

    def GetStationarityCounter(self):
        return SDP.Counter(
            steps_in_cycle=self.project_parameters["time_steps_per_stationarity_step"].GetInt(),
            beginning_step=1,
            is_active=self.project_parameters["stationary_problem_option"].GetBool())

    def GetRecoveryCounter(self):
        there_is_something_to_recover = (
            self.project_parameters["coupling_level_type"].GetInt() or
            self.project_parameters["print_PRESSURE_GRADIENT_option"].GetBool())
        return SDP.Counter(1, 1, there_is_something_to_recover)

    def GetHistoryForceQuadratureCounter(self):
        return SDP.Counter(
            self.pp.CFD_DEM["time_steps_per_quadrature_step"].GetInt(),
            1,
            self.pp.CFD_DEM["basset_force_type"].GetInt())

    def AdvanceInTime(self, step, time):
        self.step, self.time = self.dem_solver.AdvanceInTime(step, time)
        self.calculating_fluid_in_current_step = bool(time >= self.next_time_to_solve_fluid)
        if self.calculating_fluid_in_current_step:
            self.next_time_to_solve_fluid = self.fluid_solver.AdvanceInTime(time)
            self.fluid_step += 1

        return self.step, self.time

    def UpdateALEMeshMovement(self, time): # TODO: move to derived solver
        if self.project_parameters["ALE_option"].GetBool():
            self.rotator.RotateMesh(self.fluid_solver.main_model_part, time)
            self._GetProjectionModule().UpdateDatabase(self.CalculateMinElementSize())

    def CalculateMinElementSize(self):
        return self.h_min

    def AssessStationarity(self):
        Say("Assessing Stationarity...\n")
        self.stationarity = self.stationarity_tool.Assess(self.fluid_solver.main_model_part)
        self.stationarity_counter.Deactivate(self.stationarity)

    # Compute nodal quantities to be printed that are not generated as part of the
    # solution algorithm. For instance, the pressure gradient, which is not used for
    # the coupling but can be of interest.
    def ComputePostProcessResults(self):
        if self.project_parameters["coupling_level_type"].GetInt():
            self._GetProjectionModule().ComputePostProcessResults(self.dem_solver.spheres_model_part.ProcessInfo)

    def CannotIgnoreFluidNow(self):
        return self.solve_system and self.calculating_fluid_in_current_step

    def Predict(self):
        if self.CannotIgnoreFluidNow():
            self.fluid_solver.Predict()

    def ApplyForwardCoupling(self, alpha='None'):
        self._GetProjectionModule().ApplyForwardCoupling(alpha)

    def ApplyForwardCouplingOfVelocityToSlipVelocityOnly(self, time=None):
        self._GetProjectionModule().ApplyForwardCouplingOfVelocityToSlipVelocityOnly()

    def _GetProjectionModule(self):
        if not hasattr(self, 'projection_module'):
            self.projection_module = self._ConstructProjectionModule()
        return self.projection_module

    def SolveSolutionStep(self):
        # update possible movements of the fluid mesh
        self.UpdateALEMeshMovement(self.time)

        # Solving the fluid part
        Say('Solving Fluid... (', self.fluid_solver.main_model_part.NumberOfElements(0), 'elements )\n')
        self.solve_system = not self.project_parameters["fluid_already_calculated"].GetBool() and not self.stationarity

        if self.CannotIgnoreFluidNow():
            self.SolveFluid()
        else:
            Say("Skipping solving system for the fluid phase...\n")

        # Check for stationarity: this is useful for steady-state problems, so that
        # the calculation stops after reaching the solution.
        if self.stationarity_counter.Tick():
            self.AssessStationarity()

        self.derivative_recovery_counter.Activate(self.time > self.interaction_start_time and self.calculating_fluid_in_current_step)

        if self.derivative_recovery_counter.Tick():
            self.recovery.Recover()

        # Solving the disperse-phase component
        Say('Solving DEM... (', self.dem_solver.spheres_model_part.NumberOfElements(0), 'elements )')
        self.SolveDEM()

    def SolveFluid(self):
        self.fluid_solver.SolveSolutionStep()

    def SolveDEM(self):
        #self.PerformEmbeddedOperations() TO-DO: it's crashing

        it_is_time_to_forward_couple = (
            self.time >= self.interaction_start_time and
            self.coupling_level_type and
            (self.project_at_every_substep_option or self.calculating_fluid_in_current_step)
        )

        if it_is_time_to_forward_couple or self.first_DEM_iteration:

            if self.coupling_scheme_type == "UpdatedDEM":
                self.ApplyForwardCoupling()

            else:
                alpha = 1.0 - (self.next_time_to_solve_fluid - self.time) / self.fluid_dt
                self.ApplyForwardCoupling(alpha)

        if self.quadrature_counter.Tick():
            self.AppendValuesForTheHistoryForce()

        if self.integration_scheme in {'Hybrid_Bashforth', 'TerminalVelocityScheme'}:
            # Advance in space only
            self.dem_solver.SolveSolutionStep()
            self.ApplyForwardCouplingOfVelocityToSlipVelocityOnly(self.time)

        # Performing the time integration of the DEM part

        if self.do_solve_dem:
            self.dem_solver.SolveSolutionStep()

        self.first_DEM_iteration = False

    def AppendValuesForTheHistoryForce(self):
        if self.using_hinsberg_method:
            self.basset_force_tool.AppendIntegrandsWindow(self.dem_solver.spheres_model_part)
        elif self.project_parameters["basset_force_type"].GetInt() == 2:
            self.basset_force_tool.AppendIntegrands(self.dem_solver.spheres_model_part)

    def ImportModelPart(self): # TODO: implement this
        pass

    def GetComputingModelPart(self):
        return self.dem_solver.spheres_model_part
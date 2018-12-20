from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics
from python_solver import PythonSolver

# Import applications
import KratosMultiphysics.FluidDynamicsApplication as FluidDynamicsApplication
import KratosMultiphysics.DEMApplication as DEMApplication
import KratosMultiphysics.SwimmingDEMApplication as SwimmingDEMApplication
import swimming_DEM_procedures as SDP
import CFD_DEM_coupling
import derivative_recovery.derivative_recovery_strategy as derivative_recoverer
import math

def Say(*args):
    KratosMultiphysics.Logger.PrintInfo("DEM-FLUID", *args)
    KratosMultiphysics.Logger.Flush()

class SwimmingDEMSolver(PythonSolver):
    def _ValidateSettings(self, project_parameters):
        pass # to-do

        return project_parameters

    def __init__(self, model, project_parameters, fluid_solver, dem_solver, pp):
        # Validate settings

        self.project_parameters = self._ValidateSettings(project_parameters)
        print(self.project_parameters)
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

        # creating a projection module for the fluid-DEM coupling
        self.h_min = 0.01
        n_balls = 1
        fluid_volume = 10
        # the variable n_particles_in_depth is only relevant in 2D problems
        project_parameters.AddEmptyValue("n_particles_in_depth").SetInt(int(math.sqrt(n_balls / fluid_volume)))
        self.pp = pp

        # Call the base Python solver constructor
        super(SwimmingDEMSolver, self).__init__(model, project_parameters)
        self.project_parameters = project_parameters
        self.stationarity = False
        self.stationarity_counter = self.GetStationarityCounter()
        self.stationarity_tool = SDP.StationarityAssessmentTool(
            project_parameters["max_pressure_variation_rate_tol"].GetDouble(),
            SDP.FunctionsCalculator()
            )

        self.projection_module = CFD_DEM_coupling.ProjectionModule(
            self.fluid_solver.main_model_part,
            self.dem_solver.spheres_model_part,
            self.dem_solver.all_model_parts.Get("RigidFacePart"),
            pp.CFD_DEM,
            pp.coupling_dem_vars,
            pp.coupling_fluid_vars,
            pp.time_filtered_vars,
            flow_field=pp.field_utility,
            domain_size=pp.domain_size
            )

        self.derivative_recovery_counter = self.GetRecoveryCounter()

        self.recovery = derivative_recoverer.DerivativeRecoveryStrategy(
            pp,
            self.fluid_solver.main_model_part,
            SDP.FunctionsCalculator(pp.domain_size))

        self.quadrature_counter = self.GetHistoryForceQuadratureCounter()
        self.GetBassetForceTools()

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
        self.calculating_fluid_in_current_step = time == self.next_time_to_solve_fluid
        if self.calculating_fluid_in_current_step:
            self.next_time_to_solve_fluid = self.fluid_solver.AdvanceInTime(time)
            self.fluid_step += 1

        return self.step, self.time

    def UpdateALEMeshMovement(self, time): # TODO: move to derived solver
        if self.project_parameters["ALE_option"].GetBool():
            self.rotator.RotateMesh(self.fluid_solver.main_model_part, time)
            self.projection_module.UpdateDatabase(self.CalculateMinElementSize())

    def CalculateMinElementSize(self):
        return self.h_min

    def AssessStationarity(self):
        Say("Assessing Stationarity...\n")
        self.stationarity = self.stationarity_tool.Assess(self.fluid_solver.main_model_part)
        self.stationarity_counter.Deactivate(self.stationarity)

    def ComputePostProcessResults(self):
        if self.project_parameters["coupling_level_type"].GetInt():
            self.projection_module.ComputePostProcessResults(self.dem_solver.spheres_model_part.ProcessInfo)

    def Predict(self):
        self.fluid_solver.Predict()

    def ApplyForwardCoupling(self, alpha='None'):
        self.projection_module.ApplyForwardCoupling(alpha)

    def ApplyForwardCouplingOfVelocityToSlipVelocityOnly(self, time=None):
        self.projection_module.ApplyForwardCouplingOfVelocityToSlipVelocityOnly()

    def SolveSolutionStep(self):
        self.UpdateALEMeshMovement(self.time)
        # solving the fluid part
        Say('Solving Fluid... (', self.fluid_solver.main_model_part.NumberOfElements(0), 'elements )\n')
        solve_system = not self.project_parameters["fluid_already_calculated"].GetBool() and not self.stationarity

        if solve_system:
            self.fluid_solver.SolveSolutionStep()
            self.FinalizeSolutionStep()
            self.OutputSolutionStep()
        else:
            Say("Skipping solving system for the fluid phase ...\n")

        if self.stationarity_counter.Tick():
            self.AssessStationarity()

        self.ComputePostProcessResults()

        self.derivative_recovery_counter.Activate(self.time > self.interaction_start_time)

        if self.derivative_recovery_counter.Tick():
            self.recovery.Recover()

        Say('Solving DEM... (', self.dem_solver.spheres_model_part.NumberOfElements(0), 'elements )')
        self.SolveDEM(self.time)

    def SolveDEM(self, time):

        it_is_time_to_forward_couple = (
            time >= self.interaction_start_time and
            self.coupling_level_type and
            (self.project_at_every_substep_option or self.calculating_fluid_in_current_step)
        )

        if it_is_time_to_forward_couple:

            if self.coupling_scheme_type == "UpdatedDEM":
                self.ApplyForwardCoupling()

            else:
                alpha = 1.0 - (self.next_time_to_solve_fluid - time) / self.fluid_dt
                self.ApplyForwardCoupling(alpha)

        if self.quadrature_counter.Tick():
            self.AppendValuesForTheHistoryForce()

        if self.integration_scheme in {'Hybrid_Bashforth', 'TerminalVelocityScheme'}:
            # Advance in space only
            self.SolveDEMSolutionStep()
            self.ApplyForwardCouplingOfVelocityToSlipVelocityOnly(time)

        # performing the time integration of the DEM part

        if self.do_solve_dem:
            self.SolveDEMSolutionStep()

        self.dem_solver.FinalizeSolutionStep()

    def AppendValuesForTheHistoryForce(self):
        using_hinsberg_method = (
            self.project_parameters["basset_force_type"].GetInt() >= 3 or
            self.project_parameters["basset_force_type"].GetInt() == 1)
        if using_hinsberg_method:
            self.basset_force_tool.AppendIntegrandsWindow(self.dem_solver.spheres_model_part)
        elif self.project_parameters["basset_force_type"].GetInt() == 2:
            self.basset_force_tool.AppendIntegrands(self.dem_solver.spheres_model_part)

    def SolveDEMSolutionStep(self):
        self.dem_solver.SolveSolutionStep()

    def ImportModelParts(self): # TODO: implement this
        self.fluid_solver.ImportModelPart()
        self.dem_solver.ImportModelPart()

    def GetComputingModelPart(self):
        return self.dem_solver.spheres_model_part

    def GetBassetForceTools(self): # TODO: deprecated
        self.basset_force_tool = SwimmingDEMApplication.BassetForceTools()

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
        self.next_time_to_solve_fluid = project_parameters['problem_data']['start_time'].GetDouble()
        project_parameters = self._ValidateSettings(project_parameters)
        self.fluid_solver = fluid_solver
        self.dem_solver = dem_solver
        self.fluid_step = 0
        self.calculating_fluid_in_current_step = True
        # creating a projection module for the fluid-DEM coupling
        self.h_min = 0.01
        n_balls = 1
        fluid_volume = 10
        # the variable n_particles_in_depth is only relevant in 2D problems
        project_parameters.AddEmptyValue("n_particles_in_depth").SetInt(int(math.sqrt(n_balls / fluid_volume)))
        # creating a physical calculations module to analyse the DEM model_part

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

    def GetStationarityCounter(self):
        return SDP.Counter(
            steps_in_cycle=self.project_parameters["time_steps_per_stationarity_step"].GetInt(),
            beginning_step=1,
            is_active=self.project_parameters["stationary_problem_option"].GetBool())

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

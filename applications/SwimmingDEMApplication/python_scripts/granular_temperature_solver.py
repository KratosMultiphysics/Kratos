from KratosMultiphysics import Logger, Parameters
from KratosMultiphysics.python_solver import PythonSolver
import KratosMultiphysics.SwimmingDEMApplication as SDEM
import math
import KratosMultiphysics.SwimmingDEMApplication.swimming_DEM_procedures as SDP
import KratosMultiphysics.SwimmingDEMApplication.parameters_tools as PT
import KratosMultiphysics.SwimmingDEMApplication.CFD_DEM_coupling as CFD_DEM_coupling
import KratosMultiphysics.SwimmingDEMApplication.derivative_recovery.derivative_recovery_strategy as derivative_recoverer
import KratosMultiphysics.SwimmingDEMApplication.swimming_DEM_solver as swimming_DEM_solver

from swimming_DEM_solver import SwimmingDEMSolver

def Say(*args):
    Logger.PrintInfo("SwimmingDEM", *args)
    Logger.Flush()

class GranularTemperatureSolver(SwimmingDEMSolver):

    def _ValidateSettings(self, project_parameters):
        pass

    def __init__(self, model, project_parameters, field_utility, fluid_solver, dem_solver, variables_manager):
        # Call the base Python solver constructor
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

        self.fluid_step = 0
        self.calculating_fluid_in_current_step = True
        self.first_DEM_iteration = True
        # Call the base Python solver constructor
        super(SwimmingDEMSolver, self).__init__(model, project_parameters)

    def SolveSolutionStep(self):
        if self.CannotIgnoreFluidNow():
            self.SolveFluidSolutionStep()
        else:
            Say("Skipping solving system for the fluid phase...\n")

        # Solving the disperse-phase component
        Say('Solving DEM... (', self.dem_solver.spheres_model_part.NumberOfElements(0), 'elements )')
        self.SolveDEM()

        return True

    def SolveDEM(self):
        #self.PerformEmbeddedOperations() TO-DO: it's crashing
        if self.integration_scheme in {'Hybrid_Bashforth', 'TerminalVelocityScheme'}:
            # Advance in space only
            if self.do_solve_dem:
                self.SolveDEMSolutionStep()

        if self.do_solve_dem:
            self.SolveDEMSolutionStep()

        self.first_DEM_iteration = False

    # Compute nodal quantities to be printed that are not generated as part of the
    # solution algorithm. For instance, the pressure gradient, which is not used for
    # the coupling but can be of interest.
    def ComputePostProcessResults(self):
        if self.project_parameters["coupling"]["coupling_level_type"].GetInt():
            self._GetProjectionModule().ComputePostProcessResults(self.dem_solver.spheres_model_part.ProcessInfo)

    def CannotIgnoreFluidNow(self):
        return False

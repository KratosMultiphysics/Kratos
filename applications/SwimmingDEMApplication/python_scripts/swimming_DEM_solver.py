from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics
from python_solver import PythonSolver

# Import applications
import KratosMultiphysics.FluidDynamicsApplication as FluidDynamicsApplication
import KratosMultiphysics.DEMApplication as DEMApplication
import KratosMultiphysics.SwimmingDEMApplication as SwimmingDEMApplication

class SwimmingDEMSolver(PythonSolver):
    def _ValidateSettings(self, project_parameters):
        pass # to-do

        return project_parameters

    def __init__(self, model, project_parameters, fluid_solver, dem_solver):
        # Validate settings
        self.next_time_to_solve_fluid = project_parameters['problem_data']['start_time'].GetDouble()
        project_parameters = self._ValidateSettings(project_parameters)
        self.fluid_solver = fluid_solver
        self.dem_solver = dem_solver
        self.fluid_step = 0

        # Call the base Python solver constructor
        super(SwimmingDEMSolver, self).__init__(model, project_parameters)

    def AdvanceInTime(self, step, time):
        new_step, new_time = self.dem_solver.AdvanceInTime(step, time)
        if time == self.next_time_to_solve_fluid:
            self.next_time_to_solve_fluid = self.fluid_solver.AdvanceInTime(time)
            self.fluid_step += 1

        return new_step, new_time
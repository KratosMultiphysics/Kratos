import KratosMultiphysics as Kratos
from KratosMultiphysics import Vector
import KratosMultiphysics.SwimmingDEMApplication
from importlib import import_module
import swimming_DEM_solver

import numpy as np
BaseSolver = swimming_DEM_solver.SwimmingDEMSolver
import L2_error_projection_utility as error_projector

class FluidFractionTestSolver(BaseSolver):
    def __init__(self, model, project_parameters, field_utility, fluid_solver, dem_solver, variables_manager):
        self.project_parameters = project_parameters

        super(FluidFractionTestSolver, self).__init__(model,
                                                        project_parameters,
                                                        field_utility,
                                                        fluid_solver,
                                                        dem_solver,
                                                        variables_manager)

    def CannotIgnoreFluidNow(self):
        return self.solve_system and self.calculating_fluid_in_current_step

    def SolveFluidSolutionStep(self):
        self.ImposeVelocity()
        super(FluidFractionTestSolver, self).SolveFluidSolutionStep()

    def ImposeVelocity(self):
        for node in self.fluid_solver.main_model_part.Nodes:
            node.SetSolutionStepValue(Kratos.VELOCITY_Z, 0.0)
            node.Fix(Kratos.VELOCITY_Z)

    def ConstructL2ErrorProjector(self):
        self.L2_error_projector = error_projector.L2ErrorProjectionUtility(self.fluid_solver.main_model_part)

    def ProjectL2Error(self):
        self.velocity_error_projected, self.pressure_error_projected, self.error_model_part = self.L2_error_projector.ProjectL2()
        return self.velocity_error_projected, self.pressure_error_projected, self.error_model_part

    def SolveDEM(self):
        super(FluidFractionTestSolver, self).SolveDEM()

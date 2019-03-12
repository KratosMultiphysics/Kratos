from KratosMultiphysics import *
from KratosMultiphysics.SwimmingDEMApplication import *
import swimming_DEM_solver
from swimming_DEM_solver import Say
BaseSolver = swimming_DEM_solver.SwimmingDEMSolver

class InterpolationTestSolver(BaseSolver):
    def __init__(self, model, project_parameters, field_utility, fluid_solver, dem_solver, variables_manager):
        super(InterpolationTestSolver, self).__init__(model,
                                                      project_parameters,
                                                      field_utility,
                                                      fluid_solver,
                                                      dem_solver,
                                                      variables_manager)
        self.vx = 1.0
        self.vy = 1.0
        self.vz = 1.0

    def CannotIgnoreFluidNow(self):
        return self.calculating_fluid_in_current_step

    def SolveFluidSolutionStep(self):

        for node in self.fluid_solver.main_model_part.Nodes:
            velocity = Vector(3)
            velocity[0] = self.vx * node.X
            velocity[1] = self.vy * node.Y
            velocity[2] = self.vz * node.Z
            velocity *= self.next_time_to_solve_fluid
            node.SetSolutionStepValue(VELOCITY, velocity)

    def SolveDEM(self):
        super(InterpolationTestSolver, self).SolveDEM()

    def SolveDEMSolutionStep(self):
        pass


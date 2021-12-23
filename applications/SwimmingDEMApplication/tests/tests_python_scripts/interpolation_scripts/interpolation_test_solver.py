import KratosMultiphysics as Kratos
from KratosMultiphysics import Vector
import KratosMultiphysics.SwimmingDEMApplication.swimming_DEM_solver as swimming_DEM_solver
BaseSolver = swimming_DEM_solver.SwimmingDEMSolver

class InterpolationTestSolver(BaseSolver):
    def __init__(self, model, project_parameters, field_utility, fluid_solver, dem_solver, variables_manager):
        super().__init__(model,
                        project_parameters,
                        field_utility,
                        fluid_solver,
                        dem_solver,
                        variables_manager)

    def ReturnExactVelocity(self, t, x, y, z):
        interpolate_process_data = self.project_parameters['processes']['check_interpolated_fluid_velocity'][0]
        interpolate_process_parameters = interpolate_process_data['Parameters']
        field_def = [entry.GetString() for entry in interpolate_process_parameters['value']]
        field = [eval(field_def[0]),
                eval(field_def[1]),
                eval(field_def[2])]

        return Vector(field)

    def CannotIgnoreFluidNow(self):
        return self.calculating_fluid_in_current_step

    def SolveFluidSolutionStep(self):

        for node in self.fluid_solver.main_model_part.Nodes:
            velocity = self.ReturnExactVelocity(self.next_time_to_solve_fluid, node.X, node.Y, node.Z)
            node.SetSolutionStepValue(Kratos.VELOCITY, velocity)

    def SolveDEM(self):
        import random

        for node in self.dem_solver.spheres_model_part.Nodes:
            node.X = random.random()
            node.Y = random.random()
            node.Z = random.random()

        super().SolveDEM()

    def SolveDEMSolutionStep(self):
        pass


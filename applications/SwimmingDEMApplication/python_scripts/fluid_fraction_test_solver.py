import KratosMultiphysics as Kratos
import swimming_DEM_solver

BaseSolver = swimming_DEM_solver.SwimmingDEMSolver
import L2_error_calculator_utility as L2_error_calculator
class FluidFractionTestSolver(BaseSolver):
    def __init__(self, model, project_parameters, field_utility, fluid_solver, dem_solver, variables_manager):
        """The default constructor of the class

        Keyword arguments:
        self -- It signifies an instance of a class.
        model -- the container of the fluid model part.
        project_parameters -- The whole project_parameters.
        field_utility -- The utility used (if necessary) to impose the fluid field on nodes
        fluid_solver -- Solver used to calculate the fluid phase
        dem_solver -- Solver used to calculate the particle phase
        variables_manager -- The variable management used (if it is used)
        """
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

    def ConstructL2ErrorCalculator(self):
        self.L2_error_calculator = L2_error_calculator.L2ErrorCalculatorUtility(self.fluid_solver.main_model_part, self.project_parameters)

    def CalculateL2Error(self):
        self.velocity_error_norm, self.pressure_error_norm, self.error_model_part = self.L2_error_calculator.CalculateL2()
        return self.velocity_error_norm, self.pressure_error_norm, self.error_model_part

    def SolveDEM(self):
        super(FluidFractionTestSolver, self).SolveDEM()

import KratosMultiphysics as Kratos
import swimming_DEM_solver

BaseSolver = swimming_DEM_solver.SwimmingDEMSolver
import error_norm_calculator_utility as error_norm_calculator
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

    def ConstructErrorNormCalculator(self):
        self.error_norm_calculator = error_norm_calculator.ErrorNormCalculatorUtility(self.fluid_solver.main_model_part, self.project_parameters)

    def CalculateL2ErrorNorm(self):
        self.velocity_L2_error_norm, self.pressure_L2_error_norm, self.error_model_part = self.error_norm_calculator.CalculateL2Norm()
        return self.velocity_L2_error_norm, self.pressure_L2_error_norm, self.error_model_part

    def CalculateH1ErrorSemiNorm(self):
        self.velocity_H1_error_seminorm, self.pressure_H1_error_seminorm = self.error_norm_calculator.CalculateH1SemiNorm()
        return self.velocity_H1_error_seminorm, self.pressure_H1_error_seminorm

    def _ConstructProjectionModule(self):
        pass

    def ComputePostProcessResults(self):
        pass

    def ApplyForwardCoupling(self, alpha='None'):
        pass

    def ApplyForwardCouplingOfVelocityToAuxVelocityOnly(self, alpha=None):
        pass

    def _GetProjectionModule(self):
        pass

    def SolveSolutionStep(self):
        # update possible movements of the fluid mesh
        self.UpdateALEMeshMovement(self.time)
        # Solving the fluid part
        self.solve_system = not self.project_parameters["custom_fluid"]["fluid_already_calculated"].GetBool() and not self.stationarity
        if self.CannotIgnoreFluidNow():
            self.SolveFluidSolutionStep()

        #self.recovery.Recover()

        # Solving the disperse-phase component
        #self.SolveDEM()

        return True

    def SolveDEM(self):
        super(FluidFractionTestSolver, self).SolveDEM()

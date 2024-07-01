import KratosMultiphysics as Kratos
import swimming_DEM_solver

BaseSolver = swimming_DEM_solver.SwimmingDEMSolver

class FlowPastCylinderSolver(BaseSolver):
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

        super(FlowPastCylinderSolver, self).__init__(model,
                                                        project_parameters,
                                                        field_utility,
                                                        fluid_solver,
                                                        dem_solver,
                                                        variables_manager)

    def CannotIgnoreFluidNow(self):
        return self.solve_system and self.calculating_fluid_in_current_step

    def SolveFluidSolutionStep(self):
        super(FlowPastCylinderSolver, self).SolveFluidSolutionStep()

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

        self.derivative_recovery_counter.SetActivation(self.calculating_fluid_in_current_step)

        if self.derivative_recovery_counter.Tick():
            self.recovery.Recover()

        # Solving the disperse-phase component
        self.SolveDEM()

        return True

    def SolveDEM(self):
        super(FlowPastCylinderSolver, self).SolveDEM()
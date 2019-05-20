import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools
cs_data_structure = cs_tools.cs_data_structure


def Create(parameters):
    return CoupledSolverGaussSeidel(parameters)


class CoupledSolverGaussSeidel(object):
    def __init__(self, parameters):
        self.parameters = parameters
        self.settings = parameters["settings"]

        self.echo_level = self.settings["echo_level"].GetInt()

        self._predictor = cs_tools.CreateInstance(self.parameters["predictor"])
        self._convergence_accelerator = cs_tools.CreateInstance(self.parameters["convergence_accelerator"])
        self._convergence_criterion = cs_tools.CreateInstance(self.parameters["convergence_criterion"])
        self._solver_interfaces[0] = cs_tools.CreateInstance(self.parameters["solver_interfaces"][0])
        self._solver_interfaces[1] = cs_tools.CreateInstance(self.parameters["solver_interfaces"][1])

        self._components = {self._predictor, self._convergence_accelerator, self._convergence_criterion,
                            self._solver_interfaces[0], self._solver_interfaces[1]}

    def Initialize(self):
        for component in self._components:
            component.Initialize()

    def Finalize(self):
        for component in self._components:
            component.Finalize()

    def InitializeSolutionStep(self):
        for solver_name, solver in self.participating_solvers.items():
            solver.InitializeSolutionStep()

    def Predict(self):
        for solver_name, solver in self.participating_solvers.items():
            solver.Predict()

    def SolveSolutionStep(self):
        pass

    def FinalizeSolutionStep(self):
        for solver_name, solver in self.participating_solvers.items():
            solver.FinalizeSolutionStep()

    def OutputSolutionStep(self):
        for solver_name, solver in self.participating_solvers.items():
            solver.OutputSolutionStep()

    def AdvanceInTime(self, current_time):
        new_time = 0.0
        for solver_name, solver in self.participating_solvers.items():
            new_time = max(solver.AdvanceInTime(current_time), new_time)

        if self.start_coupling_time > new_time:
            self.coupling_started = False
        else:
            self.coupling_started = True

        return new_time

    def Check(self):
        for solver_name, solver in self.participating_solvers.items():
            solver.Check()

    def PrintInfo(self):
        cs_tools.PrintInfo("The class", self.__class__.__name__," has the following participants")
        for solver_name, solver in self.participating_solvers.items():
            solver.PrintInfo()
from KratosMultiphysics.CoSimulationApplication.co_simulation_component import CoSimulationComponent
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools
cs_data_structure = cs_tools.cs_data_structure


def Create(parameters):
    return CoupledSolverGaussSeidel(parameters)


class CoupledSolverGaussSeidel(CoSimulationComponent):
    def __init__(self, parameters):
        self.parameters = parameters
        self.settings = parameters["settings"]

        self.echo_level = self.settings["echo_level"].GetInt()

        self._predictor = cs_tools.CreateInstance(self.parameters["predictor"])
        self._convergence_accelerator = cs_tools.CreateInstance(self.parameters["convergence_accelerator"])
        self._convergence_criterion = cs_tools.CreateInstance(self.parameters["convergence_criterion"])
        self._solver_interfaces = []
        self._solver_interfaces.append(cs_tools.CreateInstance(self.parameters["solver_interfaces"][0]))
        self._solver_interfaces.append(cs_tools.CreateInstance(self.parameters["solver_interfaces"][1]))

        self._components = {self._predictor, self._convergence_accelerator, self._convergence_criterion,
                            self._solver_interfaces[0], self._solver_interfaces[1]}

    def Initialize(self):
        for component in self._components:
            component.Initialize()

    def Finalize(self):
        for component in self._components:
            component.Finalize()

    def InitializeSolutionStep(self):
        for component in self._components:
            component.InitializeSolutionStep()

    def SolveSolutionStep(self):
        pass

    def FinalizeSolutionStep(self):
        for component in self._components:
            component.FinalizeSolutionStep()

    def OutputSolutionStep(self):
        for component in self._components:
            component.OutputSolutionStep()

    def Check(self):
        for component in self._components:
            component.Check()

    def PrintInfo(self):
        cs_tools.PrintInfo("The coupled solver ", self.__class__.__name__, " has the following components:")
        for component in self._components:
            component.PrintInfo()

from KratosMultiphysics.CoSimulationApplication.co_simulation_component import CoSimulationComponent
from KratosMultiphysics.CoSimulationApplication.co_simulation_interface import CoSimulationInterface
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools
cs_data_structure = cs_tools.cs_data_structure


def Create(parameters):
    return CoupledSolverGaussSeidel(parameters)


class CoupledSolverGaussSeidel(CoSimulationComponent):
    def __init__(self, parameters):
        super().__init__()

        self.parameters = parameters
        self.settings = parameters["settings"]

        self.echo_level = self.settings["echo_level"].GetInt()
        self.master_solver_index = self.settings["master_solver_index"].GetInt()

        self.predictor = cs_tools.CreateInstance(self.parameters["predictor"])
        self.convergence_accelerator = cs_tools.CreateInstance(self.parameters["convergence_accelerator"])
        self.convergence_criterion = cs_tools.CreateInstance(self.parameters["convergence_criterion"])
        self.solver_wrappers = []
        self.solver_wrappers.append(cs_tools.CreateInstance(self.parameters["solver_wrappers"][0]))
        self.solver_wrappers.append(cs_tools.CreateInstance(self.parameters["solver_wrappers"][1]))
        self.components = [self.predictor, self.convergence_accelerator, self.convergence_criterion,
                            self.solver_wrappers[0], self.solver_wrappers[1]]

        self.x = []

    def Initialize(self):
        super().Initialize()

        for component in self.components[1:]:
            component.Initialize()

        self.x = self.solver_wrappers[self.master_solver_index].GetInterfaceInput()
        self.predictor.Initialize(self.x)

    def Finalize(self):
        super().Finalize()

        for component in self.components:
            component.Finalize()

    def InitializeSolutionStep(self):
        super().InitializeSolutionStep()

        for component in self.components:
            component.InitializeSolutionStep()

    def SolveSolutionStep(self):
        # Coupling iteration loop
        while not self.convergence_criterion.IsSatisfied():
            if not self.convergence_accelerator.IsReady():
                self.x = self.predictor.Predict(self.x)
            else:
                dx = self.convergence_accelerator.Predict(r)
                self.x += dx
            y = self.solver_wrappers[0].SolveSolutionStep(self.x)
            xt = self.solver_wrappers[1].SolveSolutionStep(y)
            print(type(self.x))

            r = xt - self.x
            self.convergence_accelerator.Update(self.x, xt)
            self.convergence_criterion.Update(r)

    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()

        self.predictor.Update(self.x)
        for component in self.components:
            component.FinalizeSolutionStep()

    def OutputSolutionStep(self):
        super().OutputSolutionStep()

        for component in self.components:
            component.OutputSolutionStep()

    def Check(self):
        super().Check()

        for component in self.components:
            component.Check()

    def PrintInfo(self):
        cs_tools.PrintInfo("The coupled solver ", self.__class__.__name__, " has the following components:")
        for component in self.components:
            component.PrintInfo()

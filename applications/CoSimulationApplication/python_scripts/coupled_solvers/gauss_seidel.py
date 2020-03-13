from KratosMultiphysics.CoSimulationApplication.co_simulation_component import CoSimulationComponent
from KratosMultiphysics.CoSimulationApplication.co_simulation_interface import CoSimulationInterface
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools

import numpy as np
import time
import json

cs_data_structure = cs_tools.cs_data_structure
import matplotlib.pyplot as plt

def Create(parameters):
    return CoupledSolverGaussSeidel(parameters)


class CoupledSolverGaussSeidel(CoSimulationComponent):
    def __init__(self, parameters):
        super().__init__()

        self.parameters = parameters
        self.settings = parameters["settings"]

        self.n = self.settings["timestep_start"].GetInt()  # Time step
        self.delta_t = self.settings["delta_t"].GetDouble()  # Time step size

        self.predictor = cs_tools.CreateInstance(self.parameters["predictor"])
        self.convergence_criterion = cs_tools.CreateInstance(self.parameters["convergence_criterion"])
        self.solver_wrappers = []
        for index in range(2):
            # Add timestep_start and delta_t to solver_wrapper settings
            parameters = self.parameters["solver_wrappers"][index]
            if parameters["type"].GetString() == "solver_wrappers.mapped":
                settings = parameters["settings"]["solver_wrapper"]["settings"]
            else:
                settings = parameters["settings"]

            for key in ["timestep_start", "delta_t"]:
                if settings.Has(key):
                    cs_tools.Print(f'WARNING: parameter "{key}" is defined multiple times in JSON file', layout='warning')
                    settings.RemoveValue(key)
                settings.AddValue(key, self.settings[key])

            self.solver_wrappers.append(cs_tools.CreateInstance(parameters))

        self.components = [self.predictor, self.convergence_criterion, self.solver_wrappers[0], self.solver_wrappers[1]]

        self.x = []
        self.iteration = None  # Iteration

        self.save_iterations = False  # Set True in order to save iteration related information
        if self.save_iterations:
            self.complete_solution = None
            self.complete_solution_y = None
            self.iterations = []
            self.start_time = None
            self.stop_time = None
            self.residual = []

    def Initialize(self):
        super().Initialize()

        for component in self.components[1:]:
            component.Initialize()

        # Construct mappers if required
        index_mapped = None
        index_other = None
        for i in range(2):
            type = self.parameters["solver_wrappers"][i]["type"].GetString()
            if type == "solver_wrappers.mapped":
                index_mapped = i
            else:
                index_other = i
        if index_other is None:
            raise ValueError("Not both solvers may be mapped solvers.")
        if index_mapped is not None:
            # Construct input mapper
            interface_input_from = self.solver_wrappers[index_other].GetInterfaceOutput()
            self.solver_wrappers[index_mapped].SetInterfaceInput(interface_input_from)

            # Construct output mapper
            interface_output_to = self.solver_wrappers[index_other].GetInterfaceInput()
            self.solver_wrappers[index_mapped].SetInterfaceOutput(interface_output_to)

        # Initialize variables
        self.x = self.solver_wrappers[1].GetInterfaceOutput()
        self.predictor.Initialize(self.x)

        if self.save_iterations:
            self.start_time = time.time()

    def InitializeSolutionStep(self):
        super().InitializeSolutionStep()

        for component in self.components:
            component.InitializeSolutionStep()

        self.n += 1  # Increment time step
        self.iteration = 0

        # Print timestep
        out = f"=======================================" \
              f"====================\n" \
              f"\tTime step {self.n}\n" \
              f"=======================================" \
              f"====================\n" \
              f"Iteration\tNorm residual"
        cs_tools.PrintInfo(out)

        if self.save_iterations:
            self.residual.append([])

    def SolveSolutionStep(self):
        # Initial value
        self.x = self.predictor.Predict(self.x)
        # First coupling iteration
        y = self.solver_wrappers[0].SolveSolutionStep(self.x)
        xt = self.solver_wrappers[1].SolveSolutionStep(y)
        r = xt - self.x
        self.FinalizeIteration(r)
        # Coupling iteration loop
        while not self.convergence_criterion.IsSatisfied():
            self.x += r
            y = self.solver_wrappers[0].SolveSolutionStep(self.x)
            xt = self.solver_wrappers[1].SolveSolutionStep(y)
            r = xt - self.x
            self.FinalizeIteration(r)

    def FinalizeIteration(self, r):
        self.iteration += 1
        self.convergence_criterion.Update(r)
        # Print iteration information
        norm = np.linalg.norm(r.GetNumpyArray())
        out = f"{self.iteration:<9d}\t{norm:<22.17e}"
        cs_tools.PrintInfo(out)

        if self.save_iterations:
            self.residual[self.n - 1].append(np.linalg.norm(r.GetNumpyArray()))

    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()

        if self.save_iterations:
            timestep_solution = self.x.GetNumpyArray().reshape(-1, 1)
            y = self.solver_wrappers[0].SolveSolutionStep(self.x)
            timestep_solution_y = y.GetNumpyArray().reshape(-1, 1)
            if self.complete_solution is None:
                self.complete_solution = timestep_solution
                self.complete_solution_y = timestep_solution_y
            else:
                self.complete_solution = np.hstack((self.complete_solution, timestep_solution))
                self.complete_solution_y = np.hstack((self.complete_solution_y, timestep_solution_y))
            self.iterations.append(self.iteration)

        self.predictor.Update(self.x)
        for component in self.components:
            component.FinalizeSolutionStep()

    def OutputSolutionStep(self):
        super().OutputSolutionStep()

        for component in self.components:
            component.OutputSolutionStep()

    def Finalize(self):
        super().Finalize()

        for component in self.components:
            component.Finalize()

        if self.save_iterations:
            self.stop_time = time.time()
            type = self.parameters["type"].GetString()
            if self.parameters["settings"].Has("model"):
                model = '_' + self.parameters["settings"]["model"]["type"].GetString()
                if self.parameters["settings"]["model"]["settings"].Has("q") and model == "_coupled_solvers.models.ls":
                    q = '_q' + str(self.parameters["settings"]["model"]["settings"]["q"].GetDouble())
                else:
                    q = ''
            else:
                model = ''
                q = ''
            if self.parameters["settings"].Has("surrogate"):
                sur = '_' + self.parameters["settings"]["surrogate"]["type"].GetString()[23:]
            else:
                sur = ''

            output_name = 'result.' + type + model + q + sur
            output = {"solution": self.complete_solution, "solution_y": self.complete_solution_y,
                      "iterations": self.iterations, "time": self.stop_time - self.start_time,
                      "residual": self.residual}
            np.save(output_name, output)

    def Check(self):
        super().Check()

        for component in self.components:
            component.Check()

    def PrintInfo(self, indent):
        cs_tools.Print('\n', '\t' * indent, "The coupled solver ", self.__class__.__name__, " has the following components:")
        for component in self.components:
            component.PrintInfo(indent + 1)

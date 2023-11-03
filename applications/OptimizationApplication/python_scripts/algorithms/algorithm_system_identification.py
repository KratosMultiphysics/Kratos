from abc import ABC, abstractmethod

import scipy.sparse.linalg as ssl
import numpy as np
import plotly.express as px

from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.utilities.helper_utilities import CallOnAll

import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.controls.control import Control
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.algorithms.standardized_objective import StandardizedObjective
from KratosMultiphysics.OptimizationApplication.controls.master_control import MasterControl
from KratosMultiphysics.OptimizationApplication.algorithms.algorithm import Algorithm
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView
from KratosMultiphysics.OptimizationApplication.utilities.helper_utilities import CallOnAll
from KratosMultiphysics.OptimizationApplication.utilities.opt_convergence import CreateConvergenceCriteria
from KratosMultiphysics.OptimizationApplication.utilities.opt_line_search import CreateLineSearch
from KratosMultiphysics.OptimizationApplication.utilities.logger_utilities import TimeLogger
from KratosMultiphysics.OptimizationApplication.utilities.logger_utilities import OptimizationAlgorithmTimeLogger
from KratosMultiphysics.OptimizationApplication.execution_policies.execution_policy_decorator import ExecutionPolicyDecorator


def Factory(model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
    return AlgorithmSystemIdentification(model, parameters, optimization_problem)


class AlgorithmSystemIdentification(Algorithm):
    #     def __init__(self, optimization_problem: OptimizationProblem) -> None:
    #         self._optimization_problem = optimization_problem

    @classmethod
    def GetDefaultParameters(cls):
        return Kratos.Parameters("""{
            "module"            : "KratosMultiphysics.OptimizationApplication.algorithms",
            "type"              : "PLEASE_PROVIDE_AN_ALGORITHM_CLASS_NAME",
            "objective"         : {},
            "controls"          : [],
            "echo_level"        : 0,
            "settings"          : {
                "echo_level"      : 0,
                "line_search"     : {},
                "conv_settings"   : {}
            }
        }""")

    def __init__(self, model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
        self.model = model
        self.parameters = parameters
        self._optimization_problem = optimization_problem

        parameters.ValidateAndAssignDefaults(self.GetDefaultParameters())

        self.master_control = MasterControl()  # Need to fill it with controls

        for control_name in parameters["controls"].GetStringArray():
            control = optimization_problem.GetControl(control_name)
            self.master_control.AddControl(control)

        settings = parameters["settings"]
        settings.ValidateAndAssignDefaults(self.GetDefaultParameters()["settings"])

        self.echo_level = settings["echo_level"].GetInt()

        ComponentDataView("algorithm", self._optimization_problem).SetDataBuffer(self.GetMinimumBufferSize())

        self.__convergence_criteria = CreateConvergenceCriteria(settings["conv_settings"], self._optimization_problem)

        self.__objective = StandardizedObjective(parameters["objective"], self.master_control, self._optimization_problem)
        self.__control_field: KratosOA.CollectiveExpression = None
        self.__obj_val = None
        self.algorithm_data = None

        self.__line_search_method = CreateLineSearch(settings["line_search"], self._optimization_problem, self.__objective)

    def GetMinimumBufferSize(self) -> int:
        return 2

    def Check(self):
        pass

    def Initialize(self):
        self.converged = False
        self.__obj_val = None
        self.__objective.Initialize()
        self.__objective.Check()
        self.master_control.Initialize()
        self.__control_field = self.master_control.GetControlField()
        self.algorithm_data = ComponentDataView("algorithm", self._optimization_problem)

    def Finalize(self):
        self.__objective.Finalize()
        self.master_control.Finalize()

    def ComputeSearchDirection(self, obj_gradient_expressions_input) -> KratosOA.CollectiveExpression:

        with TimeLogger("AlgorithmSystemIdentification::ComputeSearchDirection", None, "Finished"):

            obj_gradient_expressions_output = obj_gradient_expressions_input.Clone()

            search_direction = None

            for container_expression in obj_gradient_expressions_output.GetContainerExpressions():

                container_expression *= -1

                # # Automated adjustment of the initial step
                # elements = container_expression.GetModelPart().Elements
                # gradient_vector = np.reshape(container_expression.Evaluate(), (len(elements), 1))

                # if search_direction == None:
                #     search_direction = np.zeros(shape=(len(elements)))

                # gauss_newton_likelihood = gradient_vector * self.GetCurrentObjValue()
                # gauss_newton_gradient = gradient_vector@gradient_vector.T

                # search_direction += ssl.lsmr(
                #     gauss_newton_gradient,
                #     gauss_newton_likelihood,
                #     damp=0.0,
                #     # damp=1e-26,
                # )[0]

                # search_direction *= -1.0

                # # Standrad version
                # # search_direction = container_expression.Evaluate() * -1

                # Kratos.Expression.CArrayExpressionIO.Read(container_expression, search_direction)

        self.algorithm_data.GetBufferedData()["search_direction"] = obj_gradient_expressions_output

        return obj_gradient_expressions_output

    def ComputeControlUpdate(self, alpha) -> KratosOA.CollectiveExpression:
        with TimeLogger("AlgorithmSystemIdentification::ComputeControlUpdate", None, "Finished"):
            search_direction = self.algorithm_data.GetBufferedData()["search_direction"]
            if isinstance(alpha, float):
                update = search_direction * alpha
            elif isinstance(alpha, KratosOA.CollectiveExpression):
                update = search_direction.Scale(alpha)
            self.algorithm_data.GetBufferedData()["parameter_update_in_iteration"] = update.Clone()

            # update = self.algorithm_data.GetBufferedData()["search_direction"] * alpha
            # self.algorithm_data.GetBufferedData()["parameter_update_in_iteration"] = update

        return self.UpdateControl

    def UpdateControl(self) -> KratosOA.CollectiveExpression:
        with TimeLogger("AlgorithmSystemIdentification::UpdateControl", None, "Finished"):
            update = self.algorithm_data.GetBufferedData()["parameter_update_in_iteration"]
            self.__control_field += update

            for expression in self.__control_field.GetContainerExpressions():
                Kratos.Expression.CArrayExpressionIO.Read(expression, expression.Evaluate())

            self.algorithm_data.GetBufferedData()["Stiffness_in_iteration"] = self.__control_field

        return self.__control_field

    def Output(self) -> KratosOA.CollectiveExpression:
        with TimeLogger("AlgorithmSystemIdentification::Output", None, "Finished"):
            self.CallOnAllProcesses(["output_processes"], Kratos.OutputProcess.PrintOutput)

    def GetCurrentObjValue(self) -> float:
        return self.__obj_val

    def GetCurrentControlField(self):
        return self.__control_field

    def SolveOptimizationProblem(self) -> bool:
        self.Initialize()
        with TimeLogger("Solve Optimization problem", "Start", "End"):
            self.Solve()
        return self.converged

    def Solve(self) -> bool:
        np.random.seed = 123456789
        iteration_counter = 0

        self.algorithm_data = ComponentDataView("algorithm", self._optimization_problem)
        while not self.converged:
            with OptimizationAlgorithmTimeLogger("AlgorithmSystemIdentification", self._optimization_problem.GetStep()):
                iteration_counter += 1
                self.algorithm_data.GetBufferedData()["iteration"] = iteration_counter

                self.__obj_val = self.__objective.CalculateStandardizedValue(self.__control_field)

                obj_info = self.__objective.GetInfo()
                self.algorithm_data.GetBufferedData()["std_obj_value"] = obj_info["value"]
                self.algorithm_data.GetBufferedData()["rel_obj_change[%]"] = obj_info["rel_change [%]"]
                if "abs_change [%]" in obj_info:
                    self.algorithm_data.GetBufferedData()["abs_obj_change[%]"] = obj_info["abs_change [%]"]

                obj_grad = self.__objective.CalculateStandardizedGradient()
                self.algorithm_data.GetBufferedData()["standardized_gradient"] = obj_grad

                self.ComputeSearchDirection(obj_grad)
                alpha = self.__line_search_method.ComputeStep()

                self.algorithm_data.GetBufferedData()["step_size"] = alpha

                self.ComputeControlUpdate(alpha)

                self.UpdateControl()

                self.Output()

                self.converged = self.__convergence_criteria.IsConverged()

                self._optimization_problem.AdvanceStep()

                self.Finalize()

        return self.converged

    def GetOptimizedObjectiveValue(self) -> float:
        if self.converged:
            return self.__obj_val
        else:
            raise RuntimeError("Optimization problem hasn't been solved.")
import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.algorithms.standardized_objective import StandardizedObjective
from KratosMultiphysics.OptimizationApplication.controls.master_control import MasterControl
from KratosMultiphysics.OptimizationApplication.algorithms.algorithm import Algorithm
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView
from KratosMultiphysics.OptimizationApplication.utilities.opt_convergence import CreateConvergenceCriteria
from KratosMultiphysics.OptimizationApplication.utilities.opt_line_search import CreateLineSearch
from KratosMultiphysics.OptimizationApplication.algorithms.standardized_constraint import StandardizedConstraint
from KratosMultiphysics.LinearSolversApplication.dense_linear_solver_factory import ConstructSolver
from KratosMultiphysics.OptimizationApplication.utilities.helper_utilities import CallOnAll
from KratosMultiphysics.OptimizationApplication.utilities.logger_utilities import time_decorator
from KratosMultiphysics.OptimizationApplication.utilities.logger_utilities import OptimizationAlgorithmTimeLogger


def Factory(model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
    return AlgorithmGradientProjection(model, parameters, optimization_problem)

class AlgorithmGradientProjection(Algorithm):
    """
        A classical steepest descent algorithm to solve unconstrainted optimization problems.
    """

    @classmethod
    def GetDefaultParameters(cls):
        return Kratos.Parameters("""{
            "module"            : "KratosMultiphysics.OptimizationApplication.algorithms",
            "type"              : "PLEASE_PROVIDE_AN_ALGORITHM_CLASS_NAME",
            "objective"         : {},
            "constraints"       : [],
            "controls"          : [],
            "echo_level"        : 0,
            "settings"          : {
                "echo_level"      : 0,
                "line_search"     : {},
                "conv_settings"   : {},
                "linear_solver_settings" : {},
                "correction_size" : 0.0
            }
        }""")

    def __init__(self, model:Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
        self.model = model
        self.parameters = parameters
        self._optimization_problem = optimization_problem
        parameters.ValidateAndAssignDefaults(self.GetDefaultParameters())

        self.master_control = MasterControl() # Need to fill it with controls

        for control_name in parameters["controls"].GetStringArray():
            control = optimization_problem.GetControl(control_name)
            self.master_control.AddControl(control)


        settings = parameters["settings"]
        settings.ValidateAndAssignDefaults(self.GetDefaultParameters()["settings"])

        self.echo_level = settings["echo_level"].GetInt()

        ComponentDataView("algorithm", self._optimization_problem).SetDataBuffer(self.GetMinimumBufferSize())

        self.__convergence_criteria = CreateConvergenceCriteria(settings["conv_settings"], self._optimization_problem)
        self.__line_search_method = CreateLineSearch(settings["line_search"], self._optimization_problem)

        self.__objective = StandardizedObjective(parameters["objective"], self.master_control, self._optimization_problem)
        self.__constraints_list: 'list[StandardizedConstraint]' = []
        for constraint_param in parameters["constraints"].values():
            constraint = StandardizedConstraint(constraint_param, self.master_control, self._optimization_problem)
            self.__constraints_list.append(constraint)
        self.__control_field = None
        self.__obj_val = None

        default_linear_solver_settings = Kratos.Parameters("""{
            "solver_type": "LinearSolversApplication.dense_col_piv_householder_qr"
        }""")
        settings["linear_solver_settings"].ValidateAndAssignDefaults(default_linear_solver_settings)
        self.linear_solver = ConstructSolver(settings["linear_solver_settings"])
        self.correction_size = settings["correction_size"].GetDouble()

    def GetMinimumBufferSize(self) -> int:
        return 2

    def Check(self):
        self.master_control.Check()
        self.__objective.Check()
        CallOnAll(self.__constraints_list, StandardizedConstraint.Check)

    @time_decorator()
    def Initialize(self):
        self.converged = False
        self.__obj_val = None
        self.master_control.Initialize()
        self.__objective.Initialize()
        CallOnAll(self.__constraints_list, StandardizedConstraint.Initialize)
        self.__control_field = self.master_control.GetControlField()
        self.algorithm_data = ComponentDataView("algorithm", self._optimization_problem)

    @time_decorator()
    def Finalize(self):
        self.master_control.Finalize()
        self.__objective.Finalize()
        CallOnAll(self.__constraints_list, StandardizedConstraint.Finalize)

    @time_decorator()
    def ComputeSearchDirection(self, obj_grad: KratosOA.CollectiveExpression, constr_grad: 'list[KratosOA.CollectiveExpression]') -> KratosOA.CollectiveExpression:
        active_constraints_list = [self.__constraints_list[i] for i in range(len(self.__constraints_list)) if self.__constr_value[i] >= 0.0]
        number_of_active_constraints = len(active_constraints_list)
        if not number_of_active_constraints:
            search_direction = obj_grad * -1.0
            correction = obj_grad * 0.0
        else:
            constraint_violations = Kratos.Vector(number_of_active_constraints)
            for i, active_constraint in enumerate(active_constraints_list):
                    constraint_violations[i] = active_constraint.GetScaledViolationValue()

            # compute the projected search direction and correction
            ntn = Kratos.Matrix(number_of_active_constraints, number_of_active_constraints)
            for i in range(number_of_active_constraints):
                for j in range(i, number_of_active_constraints):
                    ntn[i, j] = KratosOA.ExpressionUtils.InnerProduct(constr_grad[i], constr_grad[j])
                    ntn[j, i] = ntn[i, j]

                # get the inverse of ntn
                ntn_inverse = Kratos.Matrix(number_of_active_constraints, number_of_active_constraints)

                # create the identity matrix
                identity_matrix = Kratos.Matrix(number_of_active_constraints, number_of_active_constraints, 0.0)
                for i in range(number_of_active_constraints):
                    identity_matrix[i, i] = 1.0

                # solve for inverse of ntn
                self.linear_solver.Solve(ntn, ntn_inverse, identity_matrix)

                search_direction = - (obj_grad - self.__CollectiveListVectorProduct(constr_grad, ntn_inverse * self.__CollectiveListCollectiveProduct(constr_grad, obj_grad)))
                correction = - self.__CollectiveListVectorProduct(constr_grad, ntn_inverse * constraint_violations)
        correction_norm = KratosOA.ExpressionUtils.NormInf(correction)
        if correction_norm > self.correction_size:
            correction *= self.correction_size / correction_norm
        self.algorithm_data.GetBufferedData()["search_direction"] = search_direction.Clone()
        self.algorithm_data.GetBufferedData()["correction"] = correction.Clone()

    def __CollectiveListCollectiveProduct(self, collective_list: 'list[KratosOA.CollectiveExpression]', other_collective: KratosOA.CollectiveExpression) -> Kratos.Vector:
        result = Kratos.Vector(len(collective_list))

        for i, collective_list_item in enumerate(collective_list):
            result[i] = KratosOA.ExpressionUtils.InnerProduct(collective_list_item, other_collective)
        return result

    def __CollectiveListVectorProduct(self, collective_list: 'list[KratosOA.CollectiveExpression]', vector: Kratos.Vector) -> KratosOA.CollectiveExpression:
        if len(collective_list) != vector.Size():
            raise RuntimeError(f"Collective list size and vector size mismatch. [ Collective list size = {len(collective_list)}, vector size = {vector.Size()}]")
        if len(collective_list) == 0:
            raise RuntimeError("Collective lists cannot be empty.")

        result = collective_list[0] * 0.0
        for i, collective_list_item in enumerate(collective_list):
            result += collective_list_item * vector[i]

        return result

    @time_decorator()
    def ComputeControlUpdate(self, alpha: float) -> KratosOA.CollectiveExpression:
        search_direction = self.algorithm_data.GetBufferedData()["search_direction"]
        update = KratosOA.ExpressionUtils.Scale(search_direction, alpha) + self.algorithm_data.GetBufferedData()["correction"]
        self.algorithm_data.GetBufferedData()["control_field_update"] = update.Clone()

    @time_decorator()
    def UpdateControl(self) -> KratosOA.CollectiveExpression:
        update = self.algorithm_data.GetBufferedData()["control_field_update"]
        self.__control_field += update

    @time_decorator()
    def GetCurrentObjValue(self) -> float:
        return self.__obj_val

    @time_decorator()
    def GetCurrentControlField(self):
        return self.__control_field

    @time_decorator()
    def Output(self) -> KratosOA.CollectiveExpression:
        self.algorithm_data.GetBufferedData()["control_field"] = self.__control_field.Clone()
        for process in self._optimization_problem.GetListOfProcesses("output_processes"):
            if process.IsOutputStep():
                process.PrintOutput()

    @time_decorator()
    def Solve(self):
        while not self.converged:
            with OptimizationAlgorithmTimeLogger("Gradient Projection",self._optimization_problem.GetStep()):
                self._InitializeIteration()

                self.__obj_val = self.__objective.CalculateStandardizedValue(self.__control_field)
                obj_info = self.__objective.GetInfo()
                self.algorithm_data.GetBufferedData()["std_obj_value"] = obj_info["value"]
                self.algorithm_data.GetBufferedData()["rel_change[%]"] = obj_info["rel_change [%]"]
                if "abs_change [%]" in obj_info:
                    self.algorithm_data.GetBufferedData()["abs_change[%]"] = obj_info["abs_change [%]"]

                obj_grad = self.__objective.CalculateStandardizedGradient()

                self.__constr_value = []
                active_constr_grad = []
                for constraint in self.__constraints_list:
                    value = constraint.CalculateStandardizedValue(self.__control_field)
                    self.__constr_value.append(value)
                    constr_name = constraint.GetResponseName()
                    self.algorithm_data.GetBufferedData()[f"std_constr_{constr_name}_value"] = value
                    if value >= 0.0:
                        active_constr_grad.append(constraint.CalculateStandardizedGradient())

                self.ComputeSearchDirection(obj_grad, active_constr_grad)

                alpha = self.__line_search_method.ComputeStep()

                self.ComputeControlUpdate(alpha)

                self._FinalizeIteration()

                self.Output()

                self.UpdateControl()

                self.converged = self.__convergence_criteria.IsConverged()

                self._optimization_problem.AdvanceStep()

        return self.converged

    def GetOptimizedObjectiveValue(self) -> float:
        if self.converged:
            return self.__obj_val
        else:
            raise RuntimeError("Optimization problem hasn't been solved.")

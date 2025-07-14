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
from KratosMultiphysics.OptimizationApplication.utilities.list_collective_expression_utilities import CollectiveListCollectiveProduct
from KratosMultiphysics.OptimizationApplication.utilities.list_collective_expression_utilities import CollectiveListVectorProduct
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem_utilities import OutputGradientFields
from KratosMultiphysics import FindGlobalNodalElementalNeighboursProcess
import numpy as np

def Factory(model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
    return AlgorithmTRIAL(model, parameters, optimization_problem)

class AlgorithmTRIAL(Algorithm):
    """
        A classical steepest descent algorithm to solve unconstrained optimization problems.
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
                "conv_settings"   : {},
                "linear_solver_settings" : {},
                "correction_size" : 0.0
            }
        }""")
                #"line_search"     : {},

    def __init__(self, model:Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
        self.model = model
        self.parameters = parameters
        self._optimization_problem = optimization_problem
        parameters.ValidateAndAssignDefaults(self.GetDefaultParameters())

        self.master_control = MasterControl() # Need to fill it with controls
        self._optimization_problem.AddComponent(self.master_control)

        for control_name in parameters["controls"].GetStringArray():
            control = optimization_problem.GetControl(control_name)
            self.master_control.AddControl(control)


        settings = parameters["settings"]
        settings.ValidateAndAssignDefaults(self.GetDefaultParameters()["settings"])

        self.echo_level = settings["echo_level"].GetInt()

        ComponentDataView("algorithm", self._optimization_problem).SetDataBuffer(self.GetMinimumBufferSize())

        self.__convergence_criteria = CreateConvergenceCriteria(settings["conv_settings"], self._optimization_problem)
        # self.__line_search_method = CreateLineSearch(settings["line_search"], self._optimization_problem)

        self.__objective = StandardizedObjective(parameters["objective"], self.master_control, self._optimization_problem)
        self._optimization_problem.AddComponent(self.__objective)
        self.__constraints_list: 'list[StandardizedConstraint]' = []
        for constraint_param in parameters["constraints"].values():
            constraint = StandardizedConstraint(constraint_param, self.master_control, self._optimization_problem)
            self._optimization_problem.AddComponent(constraint)
            self.__constraints_list.append(constraint)
        self.__control_field = None
        self.__obj_val = None

        default_linear_solver_settings = Kratos.Parameters("""{
            "solver_type": "LinearSolversApplication.dense_col_piv_householder_qr"
        }""")
        settings["linear_solver_settings"].ValidateAndAssignDefaults(default_linear_solver_settings)
        self.linear_solver = ConstructSolver(settings["linear_solver_settings"])
        self.correction_size = settings["correction_size"].GetDouble()

        self.lagrange_parameter = -0.01
        self.lagrange_parameter2 = 1000
        self.alpha = 0.9

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
    def CalculateUpdate(self, shape_derivative: KratosOA.CollectiveExpression):
        constraint_violations = 0.0
        for i, constraint in enumerate(self.__constraints_list):
            if constraint.GetResponseName() == "mass":
                constraint_violations = constraint.GetScaledViolationValue()
                print(f"Violation: {constraint_violations}")


        update = KratosOA.CollectiveExpression()
        for control in self.master_control.GetListOfControls():
            if control.GetName() != "LSM_HJ":
                raise RuntimeError(f"Control {control.GetName()} does not have function \"_ComputePhiGradientNorm\"")

            shape_derivative_el = control._ConvertNodalToElemental(shape_derivative)

            d_volume = constraint_violations / control.GetElements()
            self.lagrange_parameter -= 1/self.lagrange_parameter2 * d_volume
            self.lagrange_parameter2 *= self.alpha
            shape_derivative_el += 1/self.lagrange_parameter2 * d_volume - self.lagrange_parameter

            phi_grad_norm = control._ComputePhiGradientNorm2(shape_derivative_el)
            dt = control._ComputeCFLCondition(shape_derivative_el)
            # get update value
            update.Add(phi_grad_norm * shape_derivative_el * dt)
        
        self.algorithm_data.GetBufferedData()["control_field_update"] = update.Clone()
        
    @time_decorator()
    def ComputeControlUpdate(self) -> KratosOA.CollectiveExpression:
        
        update = self.algorithm_data.GetBufferedData()["search_direction"]
        self.algorithm_data.GetBufferedData()["control_field_update"] = update.Clone()

    @time_decorator()
    def UpdateControl(self) -> KratosOA.CollectiveExpression:
        update = self.algorithm_data.GetBufferedData()["control_field_update"]
        print(f"Update value: {self.__control_field}")
        self.__control_field = KratosOA.ExpressionUtils.Collapse(self.__control_field + update)

    @time_decorator()
    def GetCurrentObjValue(self) -> float:
        return self.__obj_val

    @time_decorator()
    def GetCurrentControlField(self):
        return self.__control_field

    @time_decorator()
    def Output(self) -> KratosOA.CollectiveExpression:
        self.algorithm_data.GetBufferedData()["control_field"] = self.__control_field.Clone()
        OutputGradientFields(self.__objective, self._optimization_problem, True)
        for constraint in self.__constraints_list:
            OutputGradientFields(constraint, self._optimization_problem, constraint.IsActive())
        for process in self._optimization_problem.GetListOfProcesses("output_processes"):
            if process.IsOutputStep():
                process.PrintOutput()

    @time_decorator()
    def Solve(self):
        for control in self.master_control.GetListOfControls():
            self.shape_derivative = control.GetEmptyNodalField()
        shape_collective = KratosOA.CollectiveExpression()
        shape_collective.Add(self.shape_derivative)
        shape_response = self.__objective.GetReponse()
        shape_var = shape_response.GetImplementedPhysicalKratosVariables()[0]
        shape_dict = {shape_var: shape_collective}
        
        while not self.converged:
            with OptimizationAlgorithmTimeLogger("Gradient Projection",self._optimization_problem.GetStep()):
                self._InitializeIteration()

                self.__obj_val = self.__objective.CalculateStandardizedValue(self.__control_field)
                obj_info = self.__objective.GetInfo()
                self.algorithm_data.GetBufferedData()["std_obj_value"] = obj_info["std_value"]
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
                        print(f"CONSTRAINT GRADIENT: {constraint.CalculateStandardizedGradient().Evaluate()}\nCONSTRAINT VALUE: {value}")
                print(f"OBJECTIVE GRADIENT: {obj_grad.Evaluate()}\nOBJECTIVE VALUE: {obj_info["std_value"]}")

                shape_response.CalculateGradient(shape_dict)
                shape_nodal = shape_dict[shape_var].GetContainerExpressions()[0]

                self.CalculateUpdate(shape_nodal)

                self._FinalizeIteration()

                self.Output()

                self.UpdateControl()

                converged = self.__convergence_criteria.IsConverged()

                converged &= all([value <= 0.0 for value in self.__constr_value]) 

                self.converged |=  self.__convergence_criteria.IsMaxIterationsReached()

                self._optimization_problem.AdvanceStep()

        return self.converged

    def GetOptimizedObjectiveValue(self) -> float:
        if self.converged:
            return self.__obj_val
        else:
            raise RuntimeError("Optimization problem hasn't been solved.")

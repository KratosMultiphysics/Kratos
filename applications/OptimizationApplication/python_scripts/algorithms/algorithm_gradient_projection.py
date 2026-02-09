import numpy
import KratosMultiphysics as Kratos
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.algorithms.standardized_objective import StandardizedObjective
from KratosMultiphysics.OptimizationApplication.controls.master_control import MasterControl
from KratosMultiphysics.OptimizationApplication.algorithms.algorithm import Algorithm
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView
from KratosMultiphysics.OptimizationApplication.utilities.opt_line_search import CreateLineSearch
from KratosMultiphysics.OptimizationApplication.algorithms.standardized_constraint import StandardizedConstraint
from KratosMultiphysics.LinearSolversApplication.dense_linear_solver_factory import ConstructSolver
from KratosMultiphysics.OptimizationApplication.utilities.helper_utilities import CallOnAll
from KratosMultiphysics.OptimizationApplication.utilities.logger_utilities import time_decorator
from KratosMultiphysics.OptimizationApplication.utilities.logger_utilities import OptimizationAlgorithmTimeLogger
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem_utilities import OutputGradientFields
from KratosMultiphysics.OptimizationApplication.convergence_criteria.convergence_criterion import ConvergenceCriterion
from KratosMultiphysics.OptimizationApplication.convergence_criteria.constraint_conv_criterion import ConstraintConvCriterion
from KratosMultiphysics.OptimizationApplication.convergence_criteria.combined_conv_criterion import CombinedConvCriterion
from KratosMultiphysics.OptimizationApplication.convergence_criteria.max_iter_conv_criterion import MaxIterConvCriterion
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem_utilities import OptimizationComponentFactory
from KratosMultiphysics.OptimizationApplication.utilities.logger_utilities import ListLogger


def Factory(model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
    return AlgorithmGradientProjection(model, parameters, optimization_problem)

class AlgorithmGradientProjection(Algorithm):
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
        self._optimization_problem.AddComponent(self.master_control)

        for control_name in parameters["controls"].GetStringArray():
            control = optimization_problem.GetControl(control_name)
            self.master_control.AddControl(control)


        settings = parameters["settings"]
        settings.ValidateAndAssignDefaults(self.GetDefaultParameters()["settings"])

        self.echo_level = settings["echo_level"].GetInt()

        ComponentDataView("algorithm", self._optimization_problem).SetDataBuffer(self.GetMinimumBufferSize())

        self.__line_search_method = CreateLineSearch(settings["line_search"], self._optimization_problem)

        self._objective = StandardizedObjective(parameters["objective"], self.master_control, self._optimization_problem)
        self._optimization_problem.AddComponent(self._objective)
        self._constraints_list: 'list[StandardizedConstraint]' = []
        for constraint_param in parameters["constraints"].values():
            constraint = StandardizedConstraint(constraint_param, self.master_control, self._optimization_problem)
            self._optimization_problem.AddComponent(constraint)
            self._constraints_list.append(constraint)
        self._control_field = None
        self._obj_val = None

        self._convergence_criteria = self._CreateConvergenceCriteria(settings["conv_settings"])

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
        self._objective.Check()
        CallOnAll(self._constraints_list, StandardizedConstraint.Check)

    @time_decorator()
    def Initialize(self):
        self.converged = False
        self._obj_val = None
        self.master_control.Initialize()
        self._objective.Initialize()
        CallOnAll(self._constraints_list, StandardizedConstraint.Initialize)
        self._control_field = self.master_control.GetControlField()
        self.algorithm_data = ComponentDataView("algorithm", self._optimization_problem)

        self._convergence_criteria.Initialize()

    @time_decorator()
    def Finalize(self):
        self.master_control.Finalize()
        self._objective.Finalize()
        CallOnAll(self._constraints_list, StandardizedConstraint.Finalize)

        self._convergence_criteria.Finalize()

    @time_decorator()
    def ComputeSearchDirection(self, obj_grad: Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor, constr_grad: 'list[Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor]') -> Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor:
        active_constraints_list = [self._constraints_list[i] for i in range(len(self._constraints_list)) if self.__constr_value[i] >= 0.0]
        number_of_active_constraints = len(active_constraints_list)

        search_direction = Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor(obj_grad, perform_collect_data_recursively=False, perform_store_data_recursively=False)
        correction = Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor(obj_grad, perform_collect_data_recursively=False, perform_store_data_recursively=False)

        if not number_of_active_constraints:
            search_direction.data[:] *= -1.0
            search_direction.StoreData()
            correction.data[:] = 0.0
            correction.StoreData()
        else:
            constraint_violations = Kratos.Vector(number_of_active_constraints)
            for i, active_constraint in enumerate(active_constraints_list):
                    constraint_violations[i] = active_constraint.GetScaledViolationValue()

            # compute the projected search direction and correction
            ntn = Kratos.Matrix(number_of_active_constraints, number_of_active_constraints)
            for i in range(number_of_active_constraints):
                for j in range(i, number_of_active_constraints):
                    ntn[i, j] = numpy.dot(constr_grad[i].data.ravel(), constr_grad[j].data.ravel())
                    ntn[j, i] = ntn[i, j]

            # get the inverse of ntn
            ntn_inverse = Kratos.Matrix(number_of_active_constraints, number_of_active_constraints)

            # create the identity matrix
            identity_matrix = Kratos.Matrix(number_of_active_constraints, number_of_active_constraints, 0.0)
            for i in range(number_of_active_constraints):
                identity_matrix[i, i] = 1.0

            # solve for inverse of ntn
            self.linear_solver.Solve(ntn, ntn_inverse, identity_matrix)

            search_direction.data[:] = - (obj_grad.data[:] - self._CollectiveListVectorProduct(constr_grad, ntn_inverse * self._CollectiveListCollectiveProduct(constr_grad, obj_grad)).data)
            correction.data[:] = - self._CollectiveListVectorProduct(constr_grad, ntn_inverse * constraint_violations).data

        correction_norm = numpy.linalg.norm(correction.data.ravel())
        if correction_norm > self.correction_size:
            correction.data[:] *= self.correction_size / correction_norm

        correction.StoreData()

        self.algorithm_data.GetBufferedData()["search_direction"] = search_direction
        self.algorithm_data.GetBufferedData()["correction"] = correction

    @time_decorator()
    def ComputeControlUpdate(self, alpha: Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor) -> Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor:
        search_direction: Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor = self.algorithm_data.GetBufferedData()["search_direction"]
        update = Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor(search_direction, perform_collect_data_recursively=False, perform_store_data_recursively=False)
        update.data[:] = search_direction.data * alpha.data + self.algorithm_data.GetBufferedData()["correction"].data
        Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor(update, perform_store_data_recursively=False, copy=False).StoreData()
        self.algorithm_data.GetBufferedData()["control_field_update"] = update
        return update

    @time_decorator()
    def UpdateControl(self) -> Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor:
        update: Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor = self.algorithm_data.GetBufferedData()["control_field_update"]
        self._control_field.data[:] += update.data
        Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor(self._control_field, perform_store_data_recursively=False, copy=False).StoreData()

    @time_decorator()
    def GetCurrentObjValue(self) -> float:
        return self._obj_val

    @time_decorator()
    def GetCurrentControlField(self):
        return self._control_field

    @time_decorator()
    def Output(self) -> Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor:
        self.algorithm_data.GetBufferedData()["control_field"] = Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor(self._control_field)
        OutputGradientFields(self._objective, self._optimization_problem, True)
        for constraint in self._constraints_list:
            OutputGradientFields(constraint, self._optimization_problem, constraint.IsActive())
        for process in self._optimization_problem.GetListOfProcesses("output_processes"):
            if process.IsOutputStep():
                process.PrintOutput()

    @time_decorator()
    def Solve(self):
        while not self.converged:
            with OptimizationAlgorithmTimeLogger("Gradient Projection",self._optimization_problem.GetStep()):
                self._InitializeIteration()

                self._obj_val = self._objective.CalculateStandardizedValue(self._control_field)
                obj_info = self._objective.GetInfo()
                self.algorithm_data.GetBufferedData()["std_obj_value"] = obj_info["std_value"]
                self.algorithm_data.GetBufferedData()["rel_change[%]"] = obj_info["rel_change [%]"]
                if "abs_change [%]" in obj_info:
                    self.algorithm_data.GetBufferedData()["abs_change[%]"] = obj_info["abs_change [%]"]

                obj_grad = self._objective.CalculateStandardizedGradient()

                self.__constr_value = []
                active_constr_grad = []
                for constraint in self._constraints_list:
                    value = constraint.CalculateStandardizedValue(self._control_field)
                    self.__constr_value.append(value)
                    constr_name = constraint.GetResponseName()
                    self.algorithm_data.GetBufferedData()[f"std_constr_{constr_name}_value"] = value
                    if value >= 0.0:
                        active_constr_grad.append(constraint.CalculateStandardizedGradient())

                self.ComputeSearchDirection(obj_grad, active_constr_grad)

                alpha = self.__line_search_method.ComputeStep()

                self.ComputeControlUpdate(alpha)

                self._FinalizeIteration()

                self.converged = self._convergence_criteria.IsConverged()

                self.Output()

                self.UpdateControl()

                ListLogger("Convergence info", self._convergence_criteria.GetInfo())

                self._optimization_problem.AdvanceStep()

        return self.converged

    def GetOptimizedObjectiveValue(self) -> float:
        if self.converged:
            return self._obj_val
        else:
            raise RuntimeError("Optimization problem hasn't been solved.")

    def _CreateConvergenceCriteria(self, settings: Kratos.Parameters) -> ConvergenceCriterion:
        """
        Here we create default constraint convergence criteria, if nothing
        is defined in the "constraint_conv_settings".

        """
        default_settings = Kratos.Parameters("""{
            "max_iter"                : 0,
            "objective_conv_settings" : {},
            "constraint_conv_settings": {}
        }""")

        if settings.Has("constraint_conv_settings") and settings["constraint_conv_settings"].IsString():
            default_settings["constraint_conv_settings"].SetString("")

        settings.ValidateAndAssignDefaults(default_settings)

        default_additional_conv_params = Kratos.Parameters("""{
            "type"    : "",
            "module"  : "KratosMultiphysics.OptimizationApplication.convergence_criteria",
            "settings": {}
        }""")

        max_iter_params = Kratos.Parameters("""{
            "max_iter": """ + str(settings["max_iter"].GetInt()) + """
        }""")
        max_iter_conv_criterion = MaxIterConvCriterion(max_iter_params, self._optimization_problem)

        and_conv_criteria: 'list[ConvergenceCriterion]' = []
        if not settings["objective_conv_settings"].IsEquivalentTo(Kratos.Parameters("""{}""")):
            # objective convergence criteria is given.
            current_settings = settings["objective_conv_settings"]
            current_settings.AddMissingParameters(default_additional_conv_params)
            and_conv_criteria.append(OptimizationComponentFactory(self.model, current_settings, self._optimization_problem))

        if settings["constraint_conv_settings"].IsSubParameter() and not settings["constraint_conv_settings"].IsEquivalentTo(Kratos.Parameters("""{}""")):
            # constraint convergence criteria is given
            current_settings = settings["constraint_conv_settings"]
            current_settings.AddMissingParameters(default_additional_conv_params)
            and_conv_criteria.append(OptimizationComponentFactory(self.model, current_settings, self._optimization_problem))
        else:
            if settings["constraint_conv_settings"].IsString():
                if settings["constraint_conv_settings"].GetString() != "none":
                    raise RuntimeError("constraint_conv_settings can only be either \"none\" or sub-paramter.")
            elif settings["constraint_conv_settings"].IsSubParameter():
                # constraint convergence criteria is not given, so creating them here.
                constraint_settings = Kratos.Parameters("""{
                    "component_name": ""
                }""")
                for constraint in self._constraints_list:
                    current_settings = constraint_settings.Clone()
                    current_settings["component_name"].SetString(f"response_function.{constraint.GetName()}")
                    and_conv_criteria.append(ConstraintConvCriterion(current_settings, self._optimization_problem))

        if len(and_conv_criteria) == 0:
            return max_iter_conv_criterion
        elif len(and_conv_criteria) == 1:
            # now create the combined one and return the combined one
            combined_params = Kratos.Parameters("""{
                "operator": "or"
            }""")
            combined_conv_criteria = CombinedConvCriterion(self.model, combined_params, self._optimization_problem)
            combined_conv_criteria.Add(max_iter_conv_criterion)
            combined_conv_criteria.Add(and_conv_criteria[0])
            return combined_conv_criteria
        else:
            # create the "and" combined conv criterion first
            combined_params = Kratos.Parameters("""{
                "operator": "and"
            }""")
            and_combined_conv_criteria = CombinedConvCriterion(self.model, combined_params, self._optimization_problem)
            for and_conv in and_conv_criteria:
                and_combined_conv_criteria.Add(and_conv)

            # now create the "or" combined conv criterion with the max_iter and the previous one
            combined_params = Kratos.Parameters("""{
                "operator": "or"
            }""")
            or_combined_conv_criteria = CombinedConvCriterion(self.model, combined_params, self._optimization_problem)
            or_combined_conv_criteria.Add(max_iter_conv_criterion)
            or_combined_conv_criteria.Add(and_combined_conv_criteria)
            return or_combined_conv_criteria


    @staticmethod
    def _CollectiveListCollectiveProduct(collective_list: 'list[Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor]', other_collective: Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor) -> Kratos.Vector:
        result = Kratos.Vector(len(collective_list))
        for i, collective_list_item in enumerate(collective_list):
            result[i] = numpy.dot(collective_list_item.data.ravel(), other_collective.data.ravel())
        return result

    @staticmethod
    def _CollectiveListVectorProduct(collective_list: 'list[Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor]', vector: Kratos.Vector) -> Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor:
        if len(collective_list) != vector.Size():
            raise RuntimeError(f"Collective list size and vector size mismatch. [ Collective list size = {len(collective_list)}, vector size = {vector.Size()}]")
        if len(collective_list) == 0:
            raise RuntimeError("Collective lists cannot be empty.")

        result = Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor(collective_list[0], perform_collect_data_recursively=False, perform_store_data_recursively=False)
        result.data[:] = 0.0
        for i, collective_list_item in enumerate(collective_list):
            result.data[:] += collective_list_item.data * vector[i]

        result.StoreData()
        return result
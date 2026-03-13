import numpy
import KratosMultiphysics as Kratos
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.algorithms.standardized_objective import StandardizedObjective
from KratosMultiphysics.OptimizationApplication.controls.master_control import MasterControl
from KratosMultiphysics.OptimizationApplication.algorithms.algorithm_gradient_projection import AlgorithmGradientProjection
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView
from KratosMultiphysics.OptimizationApplication.utilities.opt_line_search import CreateLineSearch
from KratosMultiphysics.OptimizationApplication.algorithms.standardized_rgp_constraint import StandardizedRGPConstraint
from KratosMultiphysics.LinearSolversApplication.dense_linear_solver_factory import ConstructSolver
from KratosMultiphysics.OptimizationApplication.utilities.logger_utilities import time_decorator
from KratosMultiphysics.OptimizationApplication.utilities.logger_utilities import OptimizationAlgorithmTimeLogger
from KratosMultiphysics.OptimizationApplication.utilities.logger_utilities import ListLogger
import math


def Factory(model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
    return AlgorithmRelaxedGradientProjection(model, parameters, optimization_problem)

class AlgorithmRelaxedGradientProjection(AlgorithmGradientProjection):
    """
        Relaxed Gradient Projection algorithm to solve constrained optimization problems. Implemnetation is based on:
        Antonau, I., Hojjat, M. & Bletzinger, KU. Relaxed gradient projection algorithm for constrained node-based shape optimization. Struct Multidisc Optim 63, 1633–1651 (2021). https://doi.org/10.1007/s00158-020-02821-y
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
                "history_size"           : 20,
                "max_inner_iter"         : 20,
                "buffer_coeff_update_factor" : 0.1
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

        self.history_size = settings["history_size"].GetInt()

        self.echo_level = settings["echo_level"].GetInt()

        ComponentDataView("algorithm", self._optimization_problem).SetDataBuffer(self.GetMinimumBufferSize())

        # self.__convergence_criteria = CreateConvergenceCriteria(settings["conv_settings"], self._optimization_problem)
        self.__line_search_method = CreateLineSearch(settings["line_search"], self._optimization_problem)

        self._objective = StandardizedObjective(parameters["objective"], self.master_control, self._optimization_problem)
        self._optimization_problem.AddComponent(self._objective)
        self._constraints_list: 'list[StandardizedRGPConstraint]' = []
        for constraint_param in parameters["constraints"].values():
            constraint = StandardizedRGPConstraint(constraint_param, self.master_control, self._optimization_problem, self.GetMinimumBufferSize())
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
        self.max_inner_iter = settings["max_inner_iter"].GetInt()
        self.buffer_coeff_update_factor = settings["buffer_coeff_update_factor"].GetDouble()

    def GetMinimumBufferSize(self) -> int:
        return self.history_size

    @time_decorator()
    def ComputeSearchDirection(self, obj_grad: Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor, constr_grad: 'list[Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor]', w_r: Kratos.Vector, w_c: Kratos.Vector) -> Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor:
        active_constraints_list = self.GetActiveConstraintsList()
        number_of_active_constraints = len(active_constraints_list)

        search_direction = Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor(obj_grad, perform_collect_data_recursively=False, perform_store_data_recursively=False)
        correction = Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor(search_direction)
        projected_direction = Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor(search_direction)

        if not number_of_active_constraints:
            search_direction.data[:] = obj_grad.data * -1.0
            correction.data[:] = 0.0
            projected_direction.data[:] = 0.0
        else:
            # scaling obj gradients
            obj_norm = numpy.linalg.norm(obj_grad.data.ravel(), ord=numpy.inf)
            scaled_obj_grad = Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor(obj_grad, perform_collect_data_recursively=False, perform_store_data_recursively=False)
            if not math.isclose(obj_norm, 0.0, abs_tol=1e-16):
                scaled_obj_grad.data[:] /= obj_norm
                scaled_obj_grad.StoreData()

            scaled_constr_grad: 'list[Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor]' = []
            lagrangian_multiplier = Kratos.Vector(number_of_active_constraints)
            for i in range(number_of_active_constraints):
                # scaling constraints grad
                norm = numpy.linalg.norm(constr_grad[i].data.ravel(), ord=numpy.inf)
                ta = Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor(constr_grad[i], perform_collect_data_recursively=False, perform_store_data_recursively=False)
                scaled_constr_grad.append(ta)
                if not math.isclose(norm, 0.0, abs_tol=1e-16):
                    ta.data[:] /= norm
                    ta.StoreData()

            # compute the projected search direction and correction
            ntn = Kratos.Matrix(number_of_active_constraints, number_of_active_constraints)
            for i in range(number_of_active_constraints):
                for j in range(i, number_of_active_constraints):
                    ntn[i, j] = numpy.dot(scaled_constr_grad[i].data.ravel(), scaled_constr_grad[j].data.ravel())
                    ntn[j, i] = ntn[i, j]

            # get the inverse of ntn
            ntn_inverse = Kratos.Matrix(number_of_active_constraints, number_of_active_constraints)

            # create the identity matrix
            identity_matrix = Kratos.Matrix(number_of_active_constraints, number_of_active_constraints, 0.0)
            for i in range(number_of_active_constraints):
                identity_matrix[i, i] = 1.0

            # solve for inverse of ntn
            self.linear_solver.Solve(ntn, ntn_inverse, identity_matrix)

            lagrangian_multiplier = ntn_inverse * self._CollectiveListCollectiveProduct(scaled_constr_grad, scaled_obj_grad)
            for i in range(number_of_active_constraints):
                if lagrangian_multiplier[i] > 0.0:
                    lagrangian_multiplier[i] = 0.0
                else:
                    lagrangian_multiplier[i] *= w_r[i]
            projected_direction.data[:] = - (scaled_obj_grad.data - self._CollectiveListVectorProduct(scaled_constr_grad, lagrangian_multiplier).data)
            correction.data[:] = - self._CollectiveListVectorProduct(scaled_constr_grad, w_c).data
            search_direction.data[:] = projected_direction.data + correction.data
        self.algorithm_data.GetBufferedData().SetValue("search_direction", search_direction, overwrite=True)
        self.algorithm_data.GetBufferedData().SetValue("correction", correction, overwrite=True)
        self.algorithm_data.GetBufferedData().SetValue("projected_direction", projected_direction, overwrite=True)

    def ComputeBufferCoefficients(self):
        active_constraints_list = self.GetActiveConstraintsList()
        number_of_active_constraints = len(active_constraints_list)
        if number_of_active_constraints == 0:
            return Kratos.Vector(), Kratos.Vector()
        else:
            w_r = Kratos.Vector(number_of_active_constraints)
            w_c = Kratos.Vector(number_of_active_constraints)
            for i, active_constraint in enumerate(active_constraints_list):
                w = active_constraint.ComputeW()
                w_r[i] = active_constraint.Compute_W_relax(w)
                w_c[i] = active_constraint.Compute_W_correction(w)
                print(f"RGP Constraint {active_constraint.GetResponseName()}:: w_r = {w_r[i]}, w_c = {w_c[i]}")
            return w_r, w_c

    @time_decorator()
    def ComputeControlUpdate(self, alpha: Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor) -> Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor:
        search_direction: Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor = self.algorithm_data.GetBufferedData()["search_direction"]
        update = Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor(search_direction, perform_collect_data_recursively=False, perform_store_data_recursively=False)
        update.data[:] = search_direction.data * alpha.data
        update.StoreData()
        self.algorithm_data.GetBufferedData().SetValue("control_field_update", update, overwrite=True)
        return update

    def CheckLinearizedConstraints(self, active_constr_grad: 'list[Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor]', update: Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor, w_r, w_c) -> bool:
        all_feasible = True
        active_constraints_list = self.GetActiveConstraintsList()
        for i, constraint in enumerate(active_constraints_list):
            predicted_value = constraint.GetStandardizedValue() + numpy.dot(active_constr_grad[i].data.ravel(), update.data.ravel())
            print(f"RGP Constraint {constraint.GetResponseName()}:: predicted g_i {predicted_value}")
            if not constraint.IsSatisfied(predicted_value):
                all_feasible = False
                w = constraint.ComputeW()
                w += self.buffer_coeff_update_factor
                w_r[i] = constraint.Compute_W_relax(w)
                w_c[i] = constraint.Compute_W_correction(w)
                print(f"RGP Constraint {constraint.GetResponseName()}:: Corrected w_r = {w_r[i]}, w_c = {w_c[i]}")
            else:
                print(f"RGP Constraint {constraint.GetResponseName()}:: no correction for w_r, w_c ")
        return all_feasible

    def GetActiveConstraintsList(self):
        return [constraint for constraint in self._constraints_list if constraint.IsActive()]

    def GetConstraintsList(self):
        return self._constraints_list

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
                    constraint.UpdateBufferSize()
                    self.__constr_value.append(value)
                    constr_name = constraint.GetResponseName()
                    self.algorithm_data.GetBufferedData()[f"std_constr_{constr_name}_value"] = value
                    msg = "active" if constraint.IsActive() else "non-active"
                    print(f"RGP Constraint {constraint.GetResponseName()} is {msg}")
                    if constraint.IsActive():
                        active_constr_grad.append(constraint.CalculateStandardizedGradient())

                w_r, w_c = self.ComputeBufferCoefficients()

                inner_iter = 0
                while inner_iter < self.max_inner_iter:

                    print(f"\n inner_iter :: {inner_iter}")

                    self.ComputeSearchDirection(obj_grad, active_constr_grad, w_r, w_c)

                    alpha = self.__line_search_method.ComputeStep()

                    update = self.ComputeControlUpdate(alpha)

                    if self.CheckLinearizedConstraints(active_constr_grad, update, w_r, w_c):
                        break

                    inner_iter += 1

                self._FinalizeIteration()

                self.converged = self._convergence_criteria.IsConverged()

                self.Output()

                self.UpdateControl()

                ListLogger("Convergence info", self._convergence_criteria.GetInfo())

                self._optimization_problem.AdvanceStep()

        return self.converged

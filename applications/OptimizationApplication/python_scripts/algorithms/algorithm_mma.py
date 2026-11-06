import numpy as np

import KratosMultiphysics as Kratos
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.algorithms.standardized_objective import StandardizedObjective
from KratosMultiphysics.OptimizationApplication.algorithms.standardized_constraint import StandardizedConstraint
from KratosMultiphysics.OptimizationApplication.controls.master_control import MasterControl
from KratosMultiphysics.OptimizationApplication.responses.response_routine import ResponseRoutine
from KratosMultiphysics.OptimizationApplication.algorithms.algorithm import Algorithm
from KratosMultiphysics.OptimizationApplication.algorithms import mma_math
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView
from KratosMultiphysics.OptimizationApplication.utilities.helper_utilities import CallOnAll
from KratosMultiphysics.OptimizationApplication.utilities.logger_utilities import time_decorator
from KratosMultiphysics.OptimizationApplication.utilities.logger_utilities import OptimizationAlgorithmTimeLogger
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem_utilities import OutputGradientFields
from KratosMultiphysics.OptimizationApplication.convergence_criteria.convergence_criterion import ConvergenceCriterion
from KratosMultiphysics.OptimizationApplication.convergence_criteria.combined_conv_criterion import CombinedConvCriterion
from KratosMultiphysics.OptimizationApplication.convergence_criteria.max_iter_conv_criterion import MaxIterConvCriterion
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem_utilities import OptimizationComponentFactory
from KratosMultiphysics.OptimizationApplication.utilities.logger_utilities import ListLogger


def Factory(model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
    return AlgorithmMMA(model, parameters, optimization_problem)


class AlgorithmMMA(Algorithm):
    """
        Method of Moving Asymptotes (MMA) and its globally convergent variant (GCMMA) for constrained problems:
        minimize f0(x) s.t. f_i(x) <= 0 (inequalities only), xmin <= x <= xmax.

        Based on [svanberg-1987-mma] K. Svanberg, Int. J. Numer. Methods Eng. 24, 359-373 (1987) and
        [svanberg-2007-mma-gcmma-note] K. Svanberg, "MMA and GCMMA - two methods for nonlinear optimization" (2007).
        Each outer iteration builds a convex separable approximation from the values and gradients at x and solves it
        (mma_math.SolveSubproblem) for the next design, so there is no line search. "gcmma" repeats the subproblem
        with a larger conservativeness parameter rho until the approximation lies above the true functions at the new
        point; gradients are evaluated once per outer iteration.

        Requirements and behaviour:
        - Exact gradients: the gradients passed to MMA must be the true dJ/dphi of the control variables. The explicit
          filter's BackwardFilterIntegratedField and the implicit filter's MapGradient return mesh-independent
          sensitivities that are not dJ/dphi, so filtered controls break the approximation (and the GCMMA conservativeness
          check). Use controls without filters, or exact (consistent) filter gradients.
        - "control_bounds" holds one {"control", "lower", "upper"} entry per control (finite bounds are needed because the
          asymptotes and move limits scale with xmax - xmin); the initial design must lie inside them.
        - Artificial variables y_i cost c_i y_i + 0.5 d_i y_i^2. Svanberg advises f0 in [1, 100], constraints scaled so
          that their limits are in [1, 100], and c_i = 1e3...1e4; if y > 0 at the solution (a warning is logged and
          "artificial_variable_max" is stored) increase c by ~100 or rescale the problem.
        - "gcmma": the control and physical fields of step k are restored to the design x_k before output, but the
          primal-solution fields (e.g. displacements) come from the last trial evaluation at x_{k+1}.
        - On the last iteration (converged) the control field is not updated to x_{k+1}: the returned field and
          objective value belong to the last evaluated design x_k.
        - Two constraints on the same response are rejected (use an equivalent expression such as "1.0*mass").
        - Sequential only: raises when run with MPI.
    """

    @classmethod
    def GetDefaultParameters(cls):
        return Kratos.Parameters("""{
            "module"     : "KratosMultiphysics.OptimizationApplication.algorithms",
            "type"       : "PLEASE_PROVIDE_AN_ALGORITHM_CLASS_NAME",
            "objective"  : {},
            "constraints": [],
            "controls"   : [],
            "echo_level" : 0,
            "settings"   : {
                "conv_settings"   : {},
                "variant"         : "mma",
                "control_bounds"  : [],
                "move"            : 0.5,
                "asyinit"         : 0.5,
                "asyincr"         : 1.2,
                "asydecr"         : 0.7,
                "asymin"          : 0.01,
                "asymax"          : 10.0,
                "albefa"          : 0.1,
                "raa0"            : 1e-5,
                "c"               : 1000.0,
                "d"               : 1.0,
                "subsolve_settings": {
                    "tolerance" : 1e-10,
                    "max_iter"  : 200
                },
                "gcmma_settings": {
                    "inner_max_iter"       : 15,
                    "conservativeness_tol" : 1e-6
                }
            }
        }""")

    def __init__(self, model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
        self.model = model
        self.parameters = parameters
        self._optimization_problem = optimization_problem

        parameters.ValidateAndAssignDefaults(self.GetDefaultParameters())

        self.master_control = MasterControl()
        self._optimization_problem.AddComponent(self.master_control)

        for control_name in parameters["controls"].GetStringArray():
            control = optimization_problem.GetControl(control_name)
            self.master_control.AddControl(control)

        settings = parameters["settings"]
        default_settings = self.GetDefaultParameters()["settings"]
        settings.ValidateAndAssignDefaults(default_settings)
        # Nested blocks are completed one by one: "conv_settings" is open, so RecursivelyValidate... would reject its keys.
        settings["subsolve_settings"].ValidateAndAssignDefaults(default_settings["subsolve_settings"])
        settings["gcmma_settings"].ValidateAndAssignDefaults(default_settings["gcmma_settings"])

        ComponentDataView("algorithm", self._optimization_problem).SetDataBuffer(self.GetMinimumBufferSize())

        self._objective = StandardizedObjective(parameters["objective"], self.master_control, self._optimization_problem)
        self._optimization_problem.AddComponent(self._objective)

        self._constraints_list: 'list[StandardizedConstraint]' = []
        for constraint_param in parameters["constraints"].values():
            constraint = StandardizedConstraint(constraint_param, self.master_control, self._optimization_problem)
            if constraint.IsEqualityType():
                raise RuntimeError("AlgorithmMMA supports inequality constraints only.")
            if constraint.GetResponseName() in [other.GetResponseName() for other in self._constraints_list]:
                raise RuntimeError(f"AlgorithmMMA does not support two constraints on the same response \"{constraint.GetResponseName()}\" "
                                   f"(output keys and response data would collide). Use an equivalent expression for the "
                                   f"second one, e.g. \"1.0*{constraint.GetResponseName()}\".")
            self._optimization_problem.AddComponent(constraint)
            self._constraints_list.append(constraint)

        self._convergence_criteria = self.__CreateConvergenceCriteria(settings["conv_settings"])

        self.variant = settings["variant"].GetString()
        if self.variant not in ("mma", "gcmma"):
            raise RuntimeError(f"Unsupported MMA \"variant\" = \"{self.variant}\". Supported options: \"mma\", \"gcmma\".")

        self.move = settings["move"].GetDouble()
        self.asyinit = settings["asyinit"].GetDouble()
        self.asyincr = settings["asyincr"].GetDouble()
        self.asydecr = settings["asydecr"].GetDouble()
        self.asymin = settings["asymin"].GetDouble()
        self.asymax = settings["asymax"].GetDouble()
        self.asymptote_margin = settings["albefa"].GetDouble()
        self.raa0 = settings["raa0"].GetDouble()
        if self.raa0 <= 0.0:
            raise RuntimeError(f"AlgorithmMMA requires \"raa0\" > 0, otherwise variables with zero gradient drift [ raa0 = {self.raa0} ].")

        self.__control_bounds = self.__ReadControlBounds(settings["control_bounds"])

        self.__c_value = settings["c"].GetDouble()
        self.__d_value = settings["d"].GetDouble()
        if self.__c_value < 0.0 or self.__d_value <= 0.0:
            raise RuntimeError(f"AlgorithmMMA requires \"c\" >= 0 and \"d\" > 0 [ c = {self.__c_value}, d = {self.__d_value} ].")

        self.__subsolve_tolerance = settings["subsolve_settings"]["tolerance"].GetDouble()
        self.__subsolve_max_iter = settings["subsolve_settings"]["max_iter"].GetInt()

        self.inner_max_iter = settings["gcmma_settings"]["inner_max_iter"].GetInt()
        self.conservativeness_tol = settings["gcmma_settings"]["conservativeness_tol"].GetDouble()

        self._control_field = None
        self._obj_val = None
        self.__ResetHistory()

    def __ReadControlBounds(self, bounds_settings: Kratos.Parameters) -> 'dict[str, tuple[float, float]]':
        bounds: 'dict[str, tuple[float, float]]' = {}
        for entry in bounds_settings.values():
            if not (entry.Has("control") and entry.Has("lower") and entry.Has("upper")):
                raise RuntimeError(f"Each \"control_bounds\" entry requires \"control\", \"lower\" and \"upper\" [ entry = {entry} ].")
            name, lower, upper = entry["control"].GetString(), entry["lower"].GetDouble(), entry["upper"].GetDouble()
            if name in bounds:
                raise RuntimeError(f"Duplicated \"control_bounds\" entry for control \"{name}\".")
            if not lower < upper:
                raise RuntimeError(f"\"control_bounds\" of control \"{name}\" requires lower < upper [ lower = {lower}, upper = {upper} ].")
            bounds[name] = (lower, upper)

        control_names = [control.GetName() for control in self.master_control.GetListOfControls()]
        for name in control_names:
            if name not in bounds:
                raise RuntimeError(f"AlgorithmMMA requires a \"control_bounds\" entry for control \"{name}\". Asymptotes and move limits scale with (upper - lower), so unbounded controls are not supported.")
        for name in bounds:
            if name not in control_names:
                raise RuntimeError(f"\"control_bounds\" refers to unknown control \"{name}\" [ controls = {control_names} ].")
        return bounds

    def __ResetHistory(self) -> None:
        self.__x_prev: 'np.ndarray | None' = None
        self.__x_prev2: 'np.ndarray | None' = None
        self.__lower: 'np.ndarray | None' = None
        self.__upper: 'np.ndarray | None' = None
        self.__capped_outer_iterations = 0
        self.__inner_cap_warned = False

    def GetMinimumBufferSize(self) -> int:
        return 2

    def Check(self) -> None:
        if Kratos.ParallelEnvironment.GetDefaultDataCommunicator().Size() > 1:
            raise RuntimeError("AlgorithmMMA is sequential and does not support MPI runs.")
        self.master_control.Check()
        self._objective.Check()
        CallOnAll(self._constraints_list, StandardizedConstraint.Check)

    @time_decorator()
    def Initialize(self) -> None:
        self.converged = False
        self._obj_val = None
        self.master_control.Initialize()
        self._objective.Initialize()
        CallOnAll(self._constraints_list, StandardizedConstraint.Initialize)
        self._control_field = self.master_control.GetControlField()
        self.algorithm_data = ComponentDataView("algorithm", self._optimization_problem)
        self._convergence_criteria.Initialize()

        # Per-control bounds expanded to the layout of the combined control field (one block per control).
        xmin_blocks, xmax_blocks = [], []
        for control, ta in zip(self.master_control.GetListOfControls(), self._control_field.GetTensorAdaptors()):
            lower, upper = self.__control_bounds[control.GetName()]
            xmin_blocks.append(np.full(ta.data.size, lower))
            xmax_blocks.append(np.full(ta.data.size, upper))
        self.xmin = np.concatenate(xmin_blocks)
        self.xmax = np.concatenate(xmax_blocks)

        x0 = self._control_field.data.reshape(-1)
        if np.any(x0 < self.xmin) or np.any(x0 > self.xmax):
            raise RuntimeError(f"The initial control field lies outside \"control_bounds\" [ min(x0 - lower) = {np.min(x0 - self.xmin)}, min(upper - x0) = {np.min(self.xmax - x0)} ].")

        m = len(self._constraints_list)
        self.c_vec = np.full(m, self.__c_value)
        self.d_vec = np.full(m, self.__d_value)
        self.__ResetHistory()

    def Finalize(self) -> None:
        self.master_control.Finalize()
        self._objective.Finalize()
        CallOnAll(self._constraints_list, StandardizedConstraint.Finalize)
        self._convergence_criteria.Finalize()

    @time_decorator()
    def ComputeControlUpdate(self, x_new: np.ndarray) -> None:
        update = Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor(self._control_field, perform_collect_data_recursively=False, perform_store_data_recursively=False)
        update.data[:] = x_new.reshape(update.data.shape) - self._control_field.data
        Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor(update, perform_store_data_recursively=False, copy=False).StoreData()
        self.algorithm_data.GetBufferedData()["control_field_update"] = update

    @time_decorator()
    def UpdateControl(self) -> None:
        update: Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor = self.algorithm_data.GetBufferedData()["control_field_update"]
        self._control_field.data[:] += update.data
        Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor(self._control_field, perform_store_data_recursively=False, copy=False).StoreData()

    @time_decorator()
    def Output(self) -> None:
        self.algorithm_data.GetBufferedData()["control_field"] = Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor(self._control_field)
        OutputGradientFields(self._objective, self._optimization_problem, True)
        for constraint in self._constraints_list:
            OutputGradientFields(constraint, self._optimization_problem, constraint.IsActive())
        for process in self._optimization_problem.GetListOfProcesses("output_processes"):
            if process.IsOutputStep():
                process.PrintOutput()

    def GetCurrentObjValue(self) -> float:
        return self._obj_val

    def GetCurrentControlField(self) -> Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor:
        return self._control_field

    def GetOptimizedObjectiveValue(self) -> float:
        if self.converged:
            return self._obj_val
        else:
            raise RuntimeError("Optimization problem hasn't been solved.")

    def __SolveSubproblem(self, p, q, r, alpha, beta):
        return mma_math.SolveSubproblem(p, q, -r[1:], self.__lower, self.__upper, alpha, beta, self.c_vec, self.d_vec,
                                        self.__subsolve_tolerance, self.__subsolve_max_iter)

    def __SolvePlainMMA(self, x, f, df, alpha, beta):
        p, q, r = mma_math.ComputeApproximation(x, f, df, self.__lower, self.__upper, self.xmin, self.xmax, self.raa0)
        x_new, y, lam, residual = self.__SolveSubproblem(p, q, r, alpha, beta)
        return x_new, y, lam, residual, 0

    def __SolveGCMMA(self, x, f, df, alpha, beta):
        rho = mma_math.InitialRho(df, self.xmin, self.xmax)
        for inner_iterations in range(self.inner_max_iter + 1):
            p, q, r = mma_math.ComputeApproximation(x, f, df, self.__lower, self.__upper, self.xmin, self.xmax, rho)
            x_new, y, lam, residual = self.__SolveSubproblem(p, q, r, alpha, beta)
            f_approx = mma_math.EvaluateApproximation(x_new, p, q, r, self.__lower, self.__upper)

            # The trial design must be stored in each control (GetEmptyField + StoreData), not only in a combined wrapper.
            trial_field = self.master_control.GetEmptyField()
            trial_field.data[:] = x_new.reshape(trial_field.data.shape)
            Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor(trial_field, perform_store_data_recursively=False, copy=False).StoreData()
            f_true = np.array([self._objective.CalculateStandardizedValue(trial_field, save_value=False)] +
                              [constraint.CalculateStandardizedValue(trial_field, save_value=False) for constraint in self._constraints_list])

            if mma_math.IsConservative(f_true, f_approx, self.conservativeness_tol):
                return x_new, y, lam, residual, inner_iterations, True
            if inner_iterations < self.inner_max_iter:
                rho_scale = mma_math.RhoScale(x_new, x, self.__lower, self.__upper, self.xmin, self.xmax)
                rho = mma_math.UpdateRho(rho, f_true, f_approx, rho_scale, self.conservativeness_tol)
        return x_new, y, lam, residual, self.inner_max_iter, False

    @time_decorator()
    def Solve(self) -> bool:
        n = self._control_field.data.size
        buffered_data = self.algorithm_data.GetBufferedData()
        while not self.converged:
            with OptimizationAlgorithmTimeLogger("AlgorithmMMA", self._optimization_problem.GetStep()):
                self._InitializeIteration()

                x = self._control_field.data.reshape(-1).copy()

                self._obj_val = self._objective.CalculateStandardizedValue(self._control_field)
                obj_info = self._objective.GetInfo()
                buffered_data["std_obj_value"] = obj_info["std_value"]
                buffered_data["rel_change[%]"] = obj_info["rel_change [%]"]
                if "abs_change [%]" in obj_info:
                    buffered_data["abs_change[%]"] = obj_info["abs_change [%]"]
                f_list = [self._obj_val]
                df_list = [self._objective.CalculateStandardizedGradient().data.reshape(-1).copy()]

                for constraint in self._constraints_list:
                    value = constraint.CalculateStandardizedValue(self._control_field)
                    buffered_data[f"std_constr_{constraint.GetResponseName()}_value"] = value
                    f_list.append(value)
                    df_list.append(constraint.CalculateStandardizedGradient().data.reshape(-1).copy())
                f = np.array(f_list)
                df = np.array(df_list).reshape(len(f_list), n)

                self.__lower, self.__upper = mma_math.UpdateAsymptotes(
                    x, self.__x_prev, self.__x_prev2, self.__lower, self.__upper, self.xmin, self.xmax,
                    self.asyinit, self.asyincr, self.asydecr, self.asymin, self.asymax)
                alpha, beta = mma_math.ComputeMoveLimits(x, self.__lower, self.__upper, self.xmin, self.xmax, self.asymptote_margin, self.move)

                if self.variant == "mma":
                    x_new, y, lam, residual, inner_iterations = self.__SolvePlainMMA(x, f, df, alpha, beta)
                else:
                    x_new, y, lam, residual, inner_iterations, conservative = self.__SolveGCMMA(x, f, df, alpha, beta)
                    self.__CheckInnerCap(conservative)
                    # Trial evaluations left the controls at the trial design; restore x_k for the output of step k.
                    self.master_control.Update(self._control_field)
                    # Cached values are dropped too: responses sharing state would otherwise give gradients from the trial design.
                    CallOnAll([self._objective, *self._constraints_list], ResponseRoutine.ResetCachedValue)
                    buffered_data["gcmma_inner_iterations"] = inner_iterations

                if residual > 1e-6:
                    Kratos.Logger.PrintWarning("AlgorithmMMA", f"The MMA subproblem relative dual residual is {residual:.3e}; the subproblem may not be solved accurately. Usually \"c\" is too large relative to the constraint scale.")
                max_y = float(np.max(y)) if y.size else 0.0
                if max_y > 1e-6:
                    Kratos.Logger.PrintWarning("AlgorithmMMA", f"Artificial variables are active (max y = {max_y:.3e}): the subproblem is infeasible within the move limits or \"c\" is too small. Rescale the problem or increase \"c\".")
                buffered_data["artificial_variable_max"] = max_y

                # True KKT residual of the original problem at x (not the subproblem residual, which is always ~0).
                buffered_data["kkt_norm"] = mma_math.KKTResidual(df[0], df[1:], f[1:], lam, x, self.xmin, self.xmax)

                self.ComputeControlUpdate(x_new)

                self._FinalizeIteration()

                self.converged = self._convergence_criteria.IsConverged()

                self.Output()

                if not self.converged:
                    self.UpdateControl()

                self.__x_prev2 = self.__x_prev
                self.__x_prev = x

                ListLogger("Convergence info", self._convergence_criteria.GetInfo())

                if not self.converged:
                    self._optimization_problem.AdvanceStep()

        return self.converged

    def __CheckInnerCap(self, conservative: bool) -> None:
        """Logs each capped GCMMA inner loop; warns once after 3 consecutive capped outer iterations (inexact gradients?)."""
        self.__capped_outer_iterations = 0 if conservative else self.__capped_outer_iterations + 1
        if not conservative:
            Kratos.Logger.PrintInfo("AlgorithmMMA", f"GCMMA inner loop stopped at inner_max_iter = {self.inner_max_iter} without a conservative approximation.")
        if self.__capped_outer_iterations >= 3 and not self.__inner_cap_warned:
            self.__inner_cap_warned = True
            Kratos.Logger.PrintWarning("AlgorithmMMA", f"GCMMA reached inner_max_iter = {self.inner_max_iter} without a conservative approximation in {self.__capped_outer_iterations} consecutive iterations. The gradients are probably not the exact dJ/dphi (e.g. filtered controls).")

    def __CreateConvergenceCriteria(self, settings: Kratos.Parameters) -> ConvergenceCriterion:
        default_settings = Kratos.Parameters("""{
            "max_iter": 0,
            "type"    : "",
            "module"  : "KratosMultiphysics.OptimizationApplication.convergence_criteria",
            "settings": {}
        }""")
        settings.AddMissingParameters(default_settings)

        max_iter_params = Kratos.Parameters("""{
            "max_iter": """ + str(settings["max_iter"].GetInt()) + """
        }""")
        max_iter_conv_criteria = MaxIterConvCriterion(max_iter_params, self._optimization_problem)
        if settings["type"].GetString() == "":
            return max_iter_conv_criteria
        else:
            additional_conv_params = Kratos.Parameters("""{
                "type"    : "",
                "module"  : "",
                "settings": {}
            }""")
            additional_conv_params["type"].SetString(settings["type"].GetString())
            additional_conv_params["module"].SetString(settings["module"].GetString())
            additional_conv_params["settings"] = settings["settings"]
            additional_conv: ConvergenceCriterion = OptimizationComponentFactory(self.model, additional_conv_params, self._optimization_problem)

            combined_params = Kratos.Parameters("""{
                "operator": "or"
            }""")
            combined_conv_criteria = CombinedConvCriterion(self.model, combined_params, self._optimization_problem)
            combined_conv_criteria.Add(max_iter_conv_criteria)
            combined_conv_criteria.Add(additional_conv)
            return combined_conv_criteria

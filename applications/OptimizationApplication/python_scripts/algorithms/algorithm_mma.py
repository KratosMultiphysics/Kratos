import numpy as np

import KratosMultiphysics as Kratos
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.algorithms.standardized_objective import StandardizedObjective
from KratosMultiphysics.OptimizationApplication.algorithms.standardized_constraint import StandardizedConstraint
from KratosMultiphysics.OptimizationApplication.controls.master_control import MasterControl
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
        A native implementation of the Method of Moving Asymptotes (MMA) and its
        globally convergent variant (GCMMA), for constrained optimization problems.

        Based on: K. Svanberg, "The Method of Moving Asymptotes - A New Method for
        Structural Optimization", Int. J. Numer. Methods Eng., Vol. 24, 359-373 (1987),
        and the general conservativeness mechanism of its GCMMA extension.

        Unlike the other native algorithms in this module, MMA solves a convex
        subproblem at every iteration to obtain the full new design point directly
        (rather than a search direction plus a separately computed step size), so it
        does not use "opt_line_search"; both the plain and globally convergent
        variants are selected via "settings.variant".
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
                "echo_level"           : 0,
                "conv_settings"        : {},
                "variant"              : "mma",
                "move"                 : 0.5,
                "asyinit"              : 0.5,
                "asyincr"              : 1.2,
                "asydecr"              : 0.7,
                "asymin"               : 0.01,
                "asymax"               : 10.0,
                "albefa"               : 0.1,
                "controls_lower_bound" : "",
                "controls_upper_bound" : "",
                "a0"                   : 1.0,
                "a"                    : 0.0,
                "c"                    : 1000.0,
                "d"                    : 1.0,
                "subsolve_settings"    : {
                    "epsilon_init"             : 1.0,
                    "epsilon_min"              : 1e-7,
                    "epsilon_reduction_factor" : 0.1,
                    "residual_tol_factor"      : 0.9,
                    "fraction_to_boundary"     : 0.99,
                    "max_newton_iter"          : 50
                },
                "gcmma_settings": {
                    "raa_init_factor"      : 0.1,
                    "raa_init_floor"       : 1e-6,
                    "raa_increase_factor"  : 2.0,
                    "raa_max"              : 1e5,
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
        settings.ValidateAndAssignDefaults(self.GetDefaultParameters()["settings"])

        self.echo_level = settings["echo_level"].GetInt()

        ComponentDataView("algorithm", self._optimization_problem).SetDataBuffer(self.GetMinimumBufferSize())

        self._objective = StandardizedObjective(parameters["objective"], self.master_control, self._optimization_problem)
        self._optimization_problem.AddComponent(self._objective)

        self._constraints_list: 'list[StandardizedConstraint]' = []
        for constraint_param in parameters["constraints"].values():
            constraint = StandardizedConstraint(constraint_param, self.master_control, self._optimization_problem)
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
        self.albefa = settings["albefa"].GetDouble()

        lower_bound_str = settings["controls_lower_bound"].GetString()
        upper_bound_str = settings["controls_upper_bound"].GetString()
        if lower_bound_str == "" or upper_bound_str == "":
            raise RuntimeError(
                "AlgorithmMMA requires finite \"controls_lower_bound\" and \"controls_upper_bound\" settings to be "
                "specified. Unlike unconstrained algorithms, MMA's moving asymptotes and move limits are scaled "
                "relative to (controls_upper_bound - controls_lower_bound) and are mathematically undefined for an "
                "unbounded control space.")
        self.__lower_bound_value = float(lower_bound_str)
        self.__upper_bound_value = float(upper_bound_str)

        self.a0 = settings["a0"].GetDouble()
        self.__a_value = settings["a"].GetDouble()
        self.__c_value = settings["c"].GetDouble()
        self.__d_value = settings["d"].GetDouble()

        subsolve_settings = settings["subsolve_settings"]
        self.subsolve_settings = {
            "epsilon_init": subsolve_settings["epsilon_init"].GetDouble(),
            "epsilon_min": subsolve_settings["epsilon_min"].GetDouble(),
            "epsilon_reduction_factor": subsolve_settings["epsilon_reduction_factor"].GetDouble(),
            "residual_tol_factor": subsolve_settings["residual_tol_factor"].GetDouble(),
            "fraction_to_boundary": subsolve_settings["fraction_to_boundary"].GetDouble(),
            "max_newton_iter": subsolve_settings["max_newton_iter"].GetInt(),
        }

        gcmma_settings = settings["gcmma_settings"]
        self.raa_init_factor = gcmma_settings["raa_init_factor"].GetDouble()
        self.raa_init_floor = gcmma_settings["raa_init_floor"].GetDouble()
        self.raa_increase_factor = gcmma_settings["raa_increase_factor"].GetDouble()
        self.raa_max = gcmma_settings["raa_max"].GetDouble()
        self.inner_max_iter = gcmma_settings["inner_max_iter"].GetInt()
        self.conservativeness_tol = gcmma_settings["conservativeness_tol"].GetDouble()

        self._control_field = None
        self._obj_val = None

    def GetMinimumBufferSize(self) -> int:
        return 2

    def Check(self) -> None:
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

        n = self._control_field.data.size
        m = len(self._constraints_list)
        self.xmin = np.full(n, self.__lower_bound_value)
        self.xmax = np.full(n, self.__upper_bound_value)
        self.a_vec = np.full(m, self.__a_value)
        self.c_vec = np.full(m, self.__c_value)
        self.d_vec = np.full(m, self.__d_value)

        # MMA's own recursion state (previous design points and asymptotes) is kept as
        # plain instance attributes -- private numerical solver state, not a cross-iteration
        # CSV-logged quantity -- mirroring how algorithm_adam.py keeps its momentum history.
        self.xold1 = None
        self.xold2 = None
        self.low = None
        self.upp = None

    def Finalize(self) -> None:
        self.master_control.Finalize()
        self._objective.Finalize()
        CallOnAll(self._constraints_list, StandardizedConstraint.Finalize)
        self._convergence_criteria.Finalize()

    @time_decorator()
    def ComputeControlUpdate(self, xmma: np.ndarray) -> None:
        update = Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor(self._control_field, perform_collect_data_recursively=False, perform_store_data_recursively=False)
        update.data[:] = xmma.reshape(update.data.shape) - self._control_field.data
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

    def __SolvePlainMMA(self, x, f0, df0, f, df, alfa, beta):
        p0, q0 = mma_math.compute_pq(df0, x, self.low, self.upp)
        p, q = mma_math.compute_pq(df, x, self.low, self.upp)
        r = mma_math.compute_r(f, x, p, q, self.low, self.upp)
        b = -r
        xmma, _, _, _, kkt_norm = mma_math.solve_mma_subproblem(
            x, alfa, beta, self.low, self.upp, p0, q0, p, q, b, self.a0, self.a_vec, self.c_vec, self.d_vec,
            self.subsolve_settings)
        return xmma, kkt_norm, 0

    def __SolveGCMMA(self, x, f0, df0, f, df, alfa, beta):
        n = x.size
        m = self.a_vec.size

        raa0 = mma_math.initial_raa(df0, self.xmin, self.xmax, self.raa_init_factor, self.raa_init_floor)
        raa = mma_math.initial_raa(df, self.xmin, self.xmax, self.raa_init_factor, self.raa_init_floor) \
            if m > 0 else np.empty(0)

        xmma = x
        kkt_norm = 0.0
        inner_iterations = 0
        reached_cap = True
        for inner_iterations in range(1, self.inner_max_iter + 1):
            p0, q0 = mma_math.compute_pq_gcmma(df0, x, self.xmin, self.xmax, self.low, self.upp, raa0)
            r0 = mma_math.compute_r(f0, x, p0, q0, self.low, self.upp)
            if m > 0:
                p, q = mma_math.compute_pq_gcmma(df, x, self.xmin, self.xmax, self.low, self.upp, raa)
                r = mma_math.compute_r(f, x, p, q, self.low, self.upp)
            else:
                p, q, r = np.empty((0, n)), np.empty((0, n)), np.empty(0)
            b = -r

            xmma, _, _, _, kkt_norm = mma_math.solve_mma_subproblem(
                x, alfa, beta, self.low, self.upp, p0, q0, p, q, b, self.a0, self.a_vec, self.c_vec, self.d_vec,
                self.subsolve_settings)

            f0_approx = mma_math.evaluate_approximation(xmma, p0, q0, r0, self.low, self.upp)
            f_approx = mma_math.evaluate_approximation(xmma, p, q, r, self.low, self.upp) if m > 0 else np.empty(0)

            # Building the trial field from a raw combined-tensor-adaptor wrapper around
            # self._control_field would only mutate that wrapper's own top-level buffer, not
            # each control's underlying storage that MasterControl.Update() (called inside
            # CalculateStandardizedValue) actually inspects. master_control.GetEmptyField() +
            # an explicit StoreData() is the pattern this codebase already uses for exactly
            # this situation (see standardized_NLOPT_objective.py's UpdateMasterControlAndLogFields).
            trial_field = self.master_control.GetEmptyField()
            trial_field.data[:] = xmma.reshape(trial_field.data.shape)
            Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor(trial_field, perform_store_data_recursively=False, copy=False).StoreData()

            f0_real = self._objective.CalculateStandardizedValue(trial_field, save_value=False)
            f_real = np.array([constraint.CalculateStandardizedValue(trial_field, save_value=False)
                                for constraint in self._constraints_list])

            if mma_math.check_conservativeness(f0_real, f_real, f0_approx, f_approx, self.conservativeness_tol):
                reached_cap = False
                break

            # Trial evaluations above never advance the optimization problem's step -- only the
            # single call at the bottom of the outer Solve() loop does. Do not add an
            # AdvanceStep() call anywhere in this inner loop.
            raa0, raa = mma_math.update_raa(raa0, raa, f0_real, f_real, f0_approx, f_approx,
                                             self.raa_increase_factor, self.raa_max)

        if reached_cap:
            Kratos.Logger.PrintWarning("AlgorithmMMA", f"GCMMA inner loop reached inner_max_iter="
                                        f"{self.inner_max_iter} without satisfying conservativeness; "
                                        "accepting the last trial point anyway.")

        return xmma, kkt_norm, inner_iterations

    @time_decorator()
    def Solve(self) -> bool:
        n = self._control_field.data.size
        while not self.converged:
            with OptimizationAlgorithmTimeLogger("AlgorithmMMA", self._optimization_problem.GetStep()):
                self._InitializeIteration()

                outer_iter = self._optimization_problem.GetStep() + 1
                x = self._control_field.data.reshape(-1).copy()

                self._obj_val = self._objective.CalculateStandardizedValue(self._control_field)
                obj_info = self._objective.GetInfo()
                self.algorithm_data.GetBufferedData()["std_obj_value"] = obj_info["std_value"]
                self.algorithm_data.GetBufferedData()["rel_obj[%]"] = obj_info["rel_change [%]"]
                if "abs_change [%]" in obj_info:
                    self.algorithm_data.GetBufferedData()["abs_obj[%]"] = obj_info["abs_change [%]"]
                df0 = self._objective.CalculateStandardizedGradient().data.reshape(-1).copy()
                f0 = self._obj_val

                f_list = []
                df_list = []
                for constraint in self._constraints_list:
                    value = constraint.CalculateStandardizedValue(self._control_field)
                    self.algorithm_data.GetBufferedData()[f"std_constr_{constraint.GetResponseName()}_value"] = value
                    f_list.append(value)
                    df_list.append(constraint.CalculateStandardizedGradient().data.reshape(-1).copy())
                f = np.array(f_list)
                df = np.array(df_list).reshape(len(f_list), n) if f_list else np.empty((0, n))

                self.low, self.upp = mma_math.update_asymptotes(
                    x, self.xold1, self.xold2, self.xmin, self.xmax, self.low, self.upp, outer_iter,
                    self.asyinit, self.asyincr, self.asydecr, self.asymin, self.asymax)
                alfa, beta = mma_math.compute_move_limits(x, self.low, self.upp, self.xmin, self.xmax, self.albefa, self.move)

                if self.variant == "mma":
                    xmma, kkt_norm, inner_iterations = self.__SolvePlainMMA(x, f0, df0, f, df, alfa, beta)
                else:
                    xmma, kkt_norm, inner_iterations = self.__SolveGCMMA(x, f0, df0, f, df, alfa, beta)

                self.algorithm_data.GetBufferedData()["kkt_norm"] = kkt_norm
                if self.variant == "gcmma":
                    self.algorithm_data.GetBufferedData()["gcmma_inner_iterations"] = inner_iterations

                self.ComputeControlUpdate(xmma)

                self._FinalizeIteration()

                self.converged = self._convergence_criteria.IsConverged()

                self.Output()

                self.UpdateControl()

                self.xold2 = self.xold1
                self.xold1 = x

                ListLogger("Convergence info", self._convergence_criteria.GetInfo())

                if not self.converged:
                    self._optimization_problem.AdvanceStep()

        return self.converged

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

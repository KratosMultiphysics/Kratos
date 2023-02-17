from typing import Union

import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.LinearSolversApplication.dense_linear_solver_factory import ConstructSolver
from KratosMultiphysics.OptimizationApplication.optimization_info import OptimizationInfo
from KratosMultiphysics.OptimizationApplication.algorithms.algorithm import Algorithm
from KratosMultiphysics.OptimizationApplication.utilities.control_transformation_technique import ControlTransformationTechnique
from KratosMultiphysics.OptimizationApplication.utilities.response_function_implementor import ObjectiveResponseFunctionImplementor
from KratosMultiphysics.OptimizationApplication.utilities.response_function_implementor import ConstraintResponseFunctionImplementor
from KratosMultiphysics.OptimizationApplication.utilities.helper_utilities import WriteCollectiveVariableDataHolderToOptmizationInfo

class GradientProjectionAlgorithm(Algorithm):
    @classmethod
    def GetDefaultParameters(cls):
        return Kratos.Parameters("""{
            "module"            : "KratosMultiphysics.OptimizationApplication.algorithms",
            "type"              : "PLEASE_PROVIDE_AN_ALGORITHM_CLASS_NAME",
            "model_part_name"   : "OptimizationModelPart",
            "control_names"     : [],
            "objectives"        : [],
            "constraints"       : [],
            "echo_level"        : 0,
            "settings"          : {
                "max_correction_share"                  : 0.75,
                "relative_tolerance"                    : 1e-3,
                "equality_constraint_relative_tolerance": 1e-3,
                "projection_step_size"                  : 1e-3,
                "correction_step_size"                  : 1e-3,
                "linear_solver_settings"                : {},
                "echo_level"                            : 0
            }
        }""")

    def __init__(self, model: Kratos.Model, parameters: Kratos.Parameters, optimization_info: OptimizationInfo):
        super().__init__(model, parameters, optimization_info)

        self.optimization_info = optimization_info
        algorithm_settings = parameters["settings"]
        algorithm_settings.ValidateAndAssignDefaults(self.GetDefaultParameters()["settings"])

        self.max_correction_share = algorithm_settings["max_correction_share"].GetDouble()
        self.relative_tolerance = algorithm_settings["relative_tolerance"].GetDouble()
        self.equality_constraint_relative_tolerance = algorithm_settings["equality_constraint_relative_tolerance"].GetDouble()
        self.projection_step_size = algorithm_settings["projection_step_size"].GetDouble()
        self.correction_step_size = algorithm_settings["correction_step_size"].GetDouble()
        self.echo_level = algorithm_settings["echo_level"].GetInt()
        self.algorithm_info_prefix = "problem_data/algorithm_data/info"

        default_linear_solver_settings = Kratos.Parameters("""{
            "solver_type": "LinearSolversApplication.dense_col_piv_householder_qr"
        }""")
        algorithm_settings["linear_solver_settings"].ValidateAndAssignDefaults(default_linear_solver_settings)
        self.linear_solver = ConstructSolver(algorithm_settings["linear_solver_settings"])

    def GetMinimumBufferSize(self) -> int:
        return 2

    def Check(self):
        if len(self.GetObjectives()) > 1:
            raise RuntimeError(f"{self.__class__.__name__} only supports single objective optimizations.")

    def Initialize(self):
        super().Initialize()
        # since we only have one objective
        self.objective = self.GetObjectives()[0]

        # set the objective and constraint info
        self.optimization_info.SetValue(f"{self.algorithm_info_prefix}/objective/weight", 1.0)
        self.optimization_info.SetValue(f"{self.algorithm_info_prefix}/objective/improvement", 1.0)
        self.optimization_info.SetValue(f"{self.algorithm_info_prefix}/objective/feasible_value", 0.0)
        self.optimization_info.SetValue(f"{self.algorithm_info_prefix}/objective/feasible_value", 0.0, 1)
        self.optimization_info.SetValue(f"{self.algorithm_info_prefix}/objective/relative_change", 0.0)
        self.optimization_info.SetValue(f"{self.algorithm_info_prefix}/objective/relative_change", 0.0, 1)
        self.optimization_info.SetValue(f"{self.algorithm_info_prefix}/constraints/feasible", True)
        self.optimization_info.SetValue(f"{self.algorithm_info_prefix}/constraints/feasible", True, 1)

        for constraint in self.GetConstraints():
            self.optimization_info.SetValue(f"{self.algorithm_info_prefix}/constraints/{constraint.GetName()}/weight", 1.0)

    def SolveSolutionStep(self) -> bool:
        self.__ComputeResponseValues()

        self.__ComputeResponseWeights()

        objective_weight = self.optimization_info.GetValue(f"{self.algorithm_info_prefix}/objective/weight")

        # calculate objective sensitivities
        compound_objective_sensitivities = KratosOA.CollectiveVariableDataHolder()
        collective_sizes: 'list[int]' = []
        for control_transformation_technique in self.GetControlTransformationTechniques():
            objective_sensitivities = KratosOA.CollectiveVariableDataHolder()
            self.__ComputeResponseSensitivityForControlTransformationTechnique(control_transformation_technique, self.objective, objective_sensitivities)
            objective_sensitivities *= (objective_weight / KratosOA.ContainerVariableDataHolderUtils.NormL2(objective_sensitivities))
            compound_objective_sensitivities.AddCollectiveVariableDataHolder(objective_sensitivities)
            collective_sizes.append((len(collective_sizes), len(collective_sizes) + len(objective_sensitivities.GetVariableDataHolders())))

        compound_objective_sensitivities /= KratosOA.ContainerVariableDataHolderUtils.NormL2(compound_objective_sensitivities)

        # compute constraint sensitivities
        compound_constraints_sensitivities: 'list[KratosOA.CollectiveVariableDataHolder()]' = []
        for constraint in self.GetConstraints():
            if constraint.IsActive():
                compound_constraint_sensitivities = KratosOA.CollectiveVariableDataHolder()
                for control_transformation_technique in self.GetControlTransformationTechniques():
                    constraint_sensitivities = KratosOA.CollectiveVariableDataHolder()
                    self.__ComputeResponseSensitivityForControlTransformationTechnique(control_transformation_technique, constraint, constraint_sensitivities)
                    constraint_sensitivities /= KratosOA.ContainerVariableDataHolderUtils.NormL2(constraint_sensitivities)
                    compound_constraint_sensitivities.AddCollectiveVariableDataHolder(constraint_sensitivities)

                # add the sensitivities to the list
                compound_constraint_sensitivities /= KratosOA.ContainerVariableDataHolderUtils.NormL2(compound_constraint_sensitivities)
                compound_constraints_sensitivities.append(compound_constraint_sensitivities)

        # computes search direction for all control transformation techniques
        projection, correction, search_direction = self.__ComputeSearchDirectionCorrection(compound_objective_sensitivities, compound_constraints_sensitivities)

        # computes control update
        for (start_index, end_index), control_transformation_technique in zip(collective_sizes, self.GetControlTransformationTechniques()):
            control_projection = KratosOA.CollectiveVariableDataHolder(projection.GetVariableDataHolders()[start_index: end_index])
            control_correction = KratosOA.CollectiveVariableDataHolder(correction.GetVariableDataHolders()[start_index: end_index])
            control_search_direction = KratosOA.CollectiveVariableDataHolder(search_direction.GetVariableDataHolders()[start_index: end_index])

            self.__ComputeControlUpdateForControlTransformationTechniques(control_transformation_technique, control_projection, control_correction, control_search_direction)

    def __ComputeResponseValues(self):
        # calculates and prints the response data
        self.__PrintInfo(1, self.objective.GetResponseInfo(), "")

        # calculate constraint values and prints info
        for constraint in self.GetConstraints():
            self.__PrintInfo(1, constraint.GetResponseInfo(), "")

    def __ComputeResponseSensitivityForControlTransformationTechnique(self, control_transformation_technique: ControlTransformationTechnique, response_function: Union[ObjectiveResponseFunctionImplementor, ConstraintResponseFunctionImplementor], output_sensitivity_collective: KratosOA.CollectiveVariableDataHolder):
        self.__PrintInfo(1, f"Computing \"{response_function.GetName()}\" sensitivities for \"{control_transformation_technique.GetName()}\" control...")

        control = control_transformation_technique.GetControl()

        control.CalculateSensitivity(response_function, output_sensitivity_collective)

        WriteCollectiveVariableDataHolderToOptmizationInfo(
            self.optimization_info,
            output_sensitivity_collective.Clone(),
            f"problem_data/response_data/<model_part_name>/{response_function.GetName()}/sensitivities/{control_transformation_technique.GetName()}/raw")

        control_transformation_technique.TransformSensitivity(output_sensitivity_collective)

        WriteCollectiveVariableDataHolderToOptmizationInfo(
            self.optimization_info,
            output_sensitivity_collective.Clone(),
            f"problem_data/response_data/<model_part_name>/{response_function.GetName()}/sensitivities/{control_transformation_technique.GetName()}/transformed")

        self.__PrintInfo(1, f"Computed \"{response_function.GetName()}\" sensitivities for \"{control_transformation_technique.GetName()}\" control.")

    def __ComputeSearchDirectionCorrection(self, objective_sensitivities: KratosOA.CollectiveVariableDataHolder, constraints_sensitivities: 'list[KratosOA.CollectiveVariableDataHolder]'):
        self.__PrintInfo(1, f"Computing search direction...")
        # create active constraints list
        active_constraints_list = [constraint for constraint in self.GetConstraints() if constraint.IsActive()]
        number_of_active_constraints = len(active_constraints_list)

        # if there are active constraints
        if number_of_active_constraints > 0:
            # create violations vector
            constraint_violations = Kratos.Vector(number_of_active_constraints)
            for i, active_constraint in enumerate(active_constraints_list):
                weight = self.optimization_info.GetValue(f"{self.algorithm_info_prefix}/constraints/{active_constraint.GetName()}/weight")
                constraint_violations[i] = weight * self.__CaomputeDivision(active_constraint.GetStandardizedValue(), active_constraint.GetReferenceValue())

            # compute the projected search direction and correction
            ntn = Kratos.Matrix(number_of_active_constraints, number_of_active_constraints)
            for i in range(number_of_active_constraints):
                for j in range(i, number_of_active_constraints):
                    ntn[i, j] = KratosOA.ContainerVariableDataHolderUtils.InnerProduct(constraints_sensitivities[i], constraints_sensitivities[j])
                    ntn[j, i] = ntn[i, j]

            # get the inverse of ntn
            ntn_inverse = Kratos.Matrix(number_of_active_constraints, number_of_active_constraints)

            # create the identity matrix
            identity_matrix = Kratos.Matrix(number_of_active_constraints, number_of_active_constraints, 0.0)
            for i in range(number_of_active_constraints):
                identity_matrix[i, i] = 1.0

            # solve for inverse of ntn
            self.linear_solver.Solve(ntn, ntn_inverse, identity_matrix)

            projection = - (objective_sensitivities - self.__CollectiveListVectorProduct(constraints_sensitivities, ntn_inverse * self.__CollectiveListCollectiveProduct(constraints_sensitivities, objective_sensitivities)))
            correction = - self.__CollectiveListVectorProduct(constraints_sensitivities, ntn_inverse * constraint_violations)
            projection_norm: float = KratosOA.ContainerVariableDataHolderUtils.NormL2(projection)
            correction_norm: float = KratosOA.ContainerVariableDataHolderUtils.NormL2(correction)

            search_direction = (projection * self.projection_step_size / projection_norm) + (correction * self.correction_step_size)
        else:
            # when no constraints are active
            projection_norm: float = KratosOA.ContainerVariableDataHolderUtils.NormL2(objective_sensitivities)
            correction_norm = 0.0

            projection = - objective_sensitivities * self.projection_step_size / projection_norm
            correction = projection.CloneWithDataInitializedToZero()
            search_direction = projection.Clone()

        # calculate search direction norm
        search_direction_norm: float = KratosOA.ContainerVariableDataHolderUtils.NormL2(search_direction)

        self.optimization_info.SetValue(f"{self.algorithm_info_prefix}/projection_norm", projection_norm)
        self.optimization_info.SetValue(f"{self.algorithm_info_prefix}/correction_norm", correction_norm)
        self.optimization_info.SetValue(f"{self.algorithm_info_prefix}/search_direction_norm", search_direction_norm)

        Kratos.Logger.PrintInfo("   ==== Algorithm info", f"projection_step_size: {self.projection_step_size}")
        Kratos.Logger.PrintInfo("                      ", f"correction_step_size: {self.correction_step_size}")
        Kratos.Logger.PrintInfo("                      ", f"projection_norm     : {projection_norm}")
        Kratos.Logger.PrintInfo("                      ", f"correction_norm     : {correction_norm}")
        Kratos.Logger.PrintInfo("                      ", f"search_norm         : {search_direction_norm}")

        self.__PrintInfo(1, f"Computed search direction.")

        return projection, correction, search_direction

    def __ComputeControlUpdateForControlTransformationTechniques(self, control_transformation_technique: ControlTransformationTechnique, projection: KratosOA.CollectiveVariableDataHolder, correction: KratosOA.CollectiveVariableDataHolder, search_direction: KratosOA.CollectiveVariableDataHolder):
        self.__PrintInfo(1, f"Computing update for \"{control_transformation_technique.GetName()}\" control...")

        WriteCollectiveVariableDataHolderToOptmizationInfo(
            self.optimization_info,
            projection.Clone(),
            f"problem_data/control_data/<model_part_name>/{control_transformation_technique.GetName()}/projection")

        WriteCollectiveVariableDataHolderToOptmizationInfo(
            self.optimization_info,
            correction.Clone(),
            f"problem_data/control_data/<model_part_name>/{control_transformation_technique.GetName()}/correction")

        WriteCollectiveVariableDataHolderToOptmizationInfo(
            self.optimization_info,
            search_direction.Clone(),
            f"problem_data/control_data/<model_part_name>/{control_transformation_technique.GetName()}/search_direction")

        control_transformation_technique.TransformUpdate(search_direction)
        control_transformation_technique.SetControlUpdate(search_direction)

        self.__PrintInfo(1, f"Computed update for \"{control_transformation_technique.GetName()}\" control.")

    def __CollectiveListCollectiveProduct(self, collective_list: 'list[KratosOA.CollectiveVariableDataHolder]', other_collective: KratosOA.CollectiveVariableDataHolder) -> Kratos.Vector:
        result = Kratos.Vector(len(collective_list))

        for i, collective_list_item in enumerate(collective_list):
            result[i] = KratosOA.ContainerVariableDataHolderUtils.InnerProduct(collective_list_item, other_collective)
        return result

    def __CollectiveListVectorProduct(self, collective_list: 'list[KratosOA.CollectiveVariableDataHolder]', vector: Kratos.Vector) -> KratosOA.CollectiveVariableDataHolder:
        if len(collective_list) != vector.Size():
            raise RuntimeError(f"Collective list size and vector size mismatch. [ Collective list size = {len(collective_list)}, vector size = {vector.Size()}]")
        if len(collective_list) == 0:
            raise RuntimeError(f"Collective lists cannot be empty.")

        result = collective_list[0].CloneWithDataInitializedToZero()
        for i, collective_list_item in enumerate(collective_list):
            result += collective_list_item * vector[i]

        return result

    def IsConverged(self) -> bool:
        if self.optimization_info["step"] > 1:
            # check for objective convergence
            objective_relative_change = abs(self.objective.GetStandardizedValue() / self.objective.GetStandardizedValue(1) - 1.0)
            is_converged = objective_relative_change < self.relative_tolerance

            # check for constraint convergence
            for constraint in self.GetConstraints():
                if constraint.GetResponseType() == "=":
                    constraint_relative_change = abs(constraint.GetStandardizedValue() / constraint.GetReferenceValue())
                    is_converged = is_converged and constraint_relative_change < self.equality_constraint_relative_tolerance
                else:
                    is_converged = is_converged and not constraint.IsActive()

            return is_converged
        else:
            return False

    def PrintConvergence(self):
        if self.optimization_info["step"] > 1:
            # check for objective convergence
            objective_relative_change = abs(self.objective.GetStandardizedValue() / self.objective.GetStandardizedValue(1) - 1.0)
            Kratos.Logger.PrintInfo("   ==== Convergence check", f"Objective \"{self.objective.GetName()}\" relative change: {objective_relative_change}")

            # check for constraint convergence
            for constraint in self.GetConstraints():
                if constraint.GetResponseType() == "=":
                    constraint_relative_change = abs(constraint.GetStandardizedValue() / constraint.GetReferenceValue())
                    Kratos.Logger.PrintInfo("                         ", f"Equality constraint \"{constraint.GetName()}\" relative change: {constraint_relative_change}")
                else:
                    Kratos.Logger.PrintInfo("                         ", f"Inequality constraint \"{constraint.GetName()}\" active: {constraint.IsActive()}")

    def __PrintInfo(self, required_echo_level: int, message: str, title="GradientProjectionAlgorithm"):
        if self.echo_level >= required_echo_level:
            Kratos.Logger.PrintInfo(title, message)

    def __CaomputeDivision(self, numerator: float, check_denominator: float, other_denominators: float = 1.0):
        if check_denominator != 0.0:
            return numerator / (check_denominator * other_denominators)
        else:
            return numerator / other_denominators

    def __ComputeResponseWeights(self):
        step = self.optimization_info["step"]

        if (step > 1):
            # first compute the objective weight
            relative_change_percentage = (self.objective.GetStandardizedValue() / self.objective.GetStandardizedValue(1) - 1.0) * 100.0
            self.optimization_info.SetValue(f"{self.algorithm_info_prefix}/objective/relative_change", relative_change_percentage)
            if relative_change_percentage < 0.0:
                ratio = 1.0 / abs(relative_change_percentage)
                ratio = max(0.8, min(1.2, ratio))
            else:
                ratio = 1.2

            objective_weight = self.optimization_info.GetValue(f"{self.algorithm_info_prefix}/objective/weight", 1)
            objective_weight *= ratio
            self.optimization_info.SetValue(f"{self.algorithm_info_prefix}/objective/weight", objective_weight)
        else:
            self.optimization_info.SetValue(f"{self.algorithm_info_prefix}/objective/weight", 1.0)
            self.optimization_info.SetValue(f"{self.algorithm_info_prefix}/objective/relative_change", 0.0)

        # now compute the weights for constraints
        current_is_feasible = True
        for constraint in self.GetConstraints():
            if constraint.IsActive():
                constraint_type = constraint.GetResponseType()
                current_violation = self.__CaomputeDivision(constraint.GetStandardizedValue(), constraint.GetReferenceValue())
                current_is_feasible = current_is_feasible and not 100.0 * abs(current_violation) > 0.5

                if (step > 1):
                    previous_violation = self.__CaomputeDivision(constraint.GetStandardizedValue(1), constraint.GetReferenceValue())
                    relative_change = self.__CaomputeDivision(constraint.GetValue(), constraint.GetValue(1)) - 1.0
                    weight = self.optimization_info.GetValue(f"{self.algorithm_info_prefix}/constraints/{constraint.GetName()}/weight", 1)

                    if constraint_type == "=":
                        if ((100 * abs(previous_violation) > 1) and (100 * abs(current_violation) > 1) and (((current_violation > 0.0) and (previous_violation < 0)) or ((current_violation < 0.0) and (previous_violation > 0)))):
                            weight *= 0.95
                        elif ((abs(current_violation) > abs(previous_violation)) and (100 * abs(previous_violation) > 0.5)):
                            weight *= 1.25
                        elif ((abs(current_violation) < abs(previous_violation)) and (100 * abs(current_violation) > 0.5) and (100 * abs(relative_change) < 1.0)):
                            weight *= 1.25
                        elif ((100 * abs(previous_violation) < 0.1) and (100 * abs(current_violation) < 0.1)):
                            weight *= 0.95
                    else:
                        if constraint.IsActive(1):
                            if ((abs(current_violation) > abs(previous_violation)) and (100 * abs(previous_violation) > 1)):
                                weight *= 1.25
                            if ((abs(current_violation) < abs(previous_violation)) and (100 * abs(relative_change) < 1.0)):
                                weight *= 1.25

                    weight = max(1.0, weight)
                    self.optimization_info.SetValue(f"{self.algorithm_info_prefix}/constraints/{constraint.GetName()}/weight", weight)
                else:
                    self.optimization_info.SetValue(f"{self.algorithm_info_prefix}/constraints/{constraint.GetName()}/weight", 1.0)

        self.optimization_info.SetValue(f"{self.algorithm_info_prefix}/constraints/feasible", current_is_feasible)
        previous_is_feasible = self.optimization_info.GetValue(f"{self.algorithm_info_prefix}/constraints/feasible", 1)
        last_feasible_objective_value = self.optimization_info.GetValue(f"{self.algorithm_info_prefix}/objective/feasible_value", 1)

        improvement = self.optimization_info.GetValue(f"{self.algorithm_info_prefix}/objective/improvement", 1)
        # apply contraction if necessary
        if (step > 2 and current_is_feasible and previous_is_feasible and (self.objective.GetStandardizedValue() > last_feasible_objective_value)):
            improvement *= 0.9

        if (step > 2 and current_is_feasible and previous_is_feasible and (self.objective.GetStandardizedValue() < last_feasible_objective_value)):
            improvement *= 1.01
        self.optimization_info.SetValue(f"{self.algorithm_info_prefix}/objective/improvement", improvement)

        if step > 1:
            previous_relative_change_percentage = self.optimization_info.GetValue(f"{self.algorithm_info_prefix}/objective/relative_change", 1)
            previous_is_feasible = self.optimization_info.GetValue(f"{self.algorithm_info_prefix}/constraints/feasible", 1)
            scale = 1.0
            ratio = improvement / abs(relative_change_percentage)

            if (ratio > 1.0 and relative_change_percentage < 0.0 and previous_relative_change_percentage < 0.0 and current_is_feasible and previous_is_feasible):
                scale = ratio
            if (ratio < 1.0 and relative_change_percentage < 0.0 and previous_relative_change_percentage < 0.0 and current_is_feasible and previous_is_feasible):
                scale = ratio

            scale = min(1.2, max(0.8, scale))

            if (relative_change_percentage > 0.0 and current_is_feasible and previous_is_feasible):
                scale = 0.8

            self.projection_step_size *= scale

        if current_is_feasible:
            self.optimization_info.SetValue(f"{self.algorithm_info_prefix}/objective/feasible_value", self.objective.GetStandardizedValue())
        else:
            if step == 1:
                self.optimization_info.SetValue(f"{self.algorithm_info_prefix}/objective/feasible_value", 0.0)
            else:
                self.optimization_info.SetValue(f"{self.algorithm_info_prefix}/objective/feasible_value", self.optimization_info.GetValue(f"{self.algorithm_info_prefix}/objective/feasible_value", 1))
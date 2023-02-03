from typing import Union

import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.LinearSolversApplication.dense_linear_solver_factory import ConstructSolver
from KratosMultiphysics.OptimizationApplication.optimization_info import OptimizationInfo
from KratosMultiphysics.OptimizationApplication.algorithms.algorithm import Algorithm
from KratosMultiphysics.OptimizationApplication.utilities.control_transformation_technique import ControlTransformationTechnique
from KratosMultiphysics.OptimizationApplication.utilities.response_function_implementor import ObjectiveResponseFunctionImplementor
from KratosMultiphysics.OptimizationApplication.utilities.response_function_implementor import ConstraintResponseFunctionImplementor
from KratosMultiphysics.OptimizationApplication.utilities.helper_utils import ContainerVariableDataHolderUnion


class GradientProjectionAlgorithm(Algorithm):
    def __init__(self, model: Kratos.Model, parameters: Kratos.Parameters, optimization_info: OptimizationInfo):
        default_settings = Kratos.Parameters("""{
            "max_correction_share"  : 0.75,
            "relative_tolerance"    : 1e-3,
            "step_size"             : 1e-3,
            "linear_solver_settings": {},
            "echo_level"            : 0
        }""")
        parameters.ValidateAndAssignDefaults(default_settings)

        self.optimization_info = optimization_info

        self.max_correction_share = parameters["max_correction_share"].GetDouble()
        self.relative_tolerance = parameters["relative_tolerance"].GetDouble()
        self.step_size = parameters["step_size"].GetDouble()
        self.echo_level = parameters["echo_level"].GetInt()

        default_linear_solver_settings = Kratos.Parameters("""{
            "solver_type": "LinearSolversApplication.dense_col_piv_householder_qr"
        }""")
        parameters["linear_solver_settings"].ValidateAndAssignDefaults(default_linear_solver_settings)
        self.linear_solver = ConstructSolver(parameters["linear_solver_settings"])

    def GetMinimumBufferSize(self) -> int:
        return 2

    def Check(self):
        if len(self.GetObjectives()) > 1:
            raise RuntimeError(f"{self.GetName()} algorithm of type {self.__class__.__name__} only supports single objective optimizations.")

    def Initialize(self):
        # since we only have one objective
        self.objective = self.GetObjectives()[0]

    def InitializeSolutionStep(self):
        self.optimization_info[self.GetName()] = {
            "objectives": {},
            "constraints": {},
            "controls": {}
        }

    def SolveSolutionStep(self) -> bool:
        # get algorithm specific data from optimization info
        algorithm_data = self.optimization_info[self.GetName()]

        self.__ComputeResponseValues()

        # calculate response gradients
        for controller in self.GetControllers():
            self.__PrintInfo(1, f"Computing sensitivities for {controller.GetName()} control...")

            # create controller optimization info
            algorithm_data["controls"][controller] = {
                "objectives": {},
                "constraints": {}
            }

            control_data = algorithm_data["controls"][controller]

            # compute objective sensitivities
            self.__ComputeResponseSensitivityForControlWrapper(controller, self.objective, control_data["objectives"])

            # compute constraint sensitivities
            for constraint, constraint_data in algorithm_data["constraints"].items():
                if constraint_data["is_active"]:
                    self.__ComputeResponseSensitivityForControlWrapper(controller, constraint, control_data["constraints"])

            # compute control update
            self.__ComputeControlUpdatesForControlWrapper(controller)

            self.__PrintInfo(1, f"Computed sensitivities for {controller.GetName()} control.")

    def IsConverged(self) -> bool:
        if self.optimization_info["step"] > 1:
            # check for objective convergence
            is_converged = abs(self.optimization_info[self.GetName()]["objectives"][self.objective]["value"] / self.optimization_info.GetSolutionStepData(1)
                               [self.GetName()]["objectives"][self.objective]["value"] - 1.0) < self.relative_tolerance

            # check for constraint convergence
            for constraint_data in self.optimization_info[self.GetName()]["constraints"].values():
                is_converged = is_converged and not constraint_data["is_active"]

            return is_converged
        else:
            return False

    def __PrintInfo(self, required_echo_level: int, message: str, title="GradientProjectionAlgorithm"):
        if self.echo_level >= required_echo_level:
            Kratos.Logger.PrintInfo(title, message)

    def __ComputeResponseValues(self):
        algorithm_data = self.optimization_info[self.GetName()]

        # calculate objective value
        algorithm_data["objectives"][self.objective] = {
            "value": self.objective.CalculateValue(),
            "standardized_value": self.objective.CalculateStandardizedValue()
        }

        relative_change = 0.0
        if self.optimization_info["step"] > 1:
            relative_change = (self.objective.CalculateValue() / self.optimization_info.GetSolutionStepData(1)[self.GetName()]["objectives"][self.objective]["value"] - 1.0) * 100.0

        self.__PrintInfo(1, self.objective.GetResponseInfo() + f"\n\t\t rel_change [%]: {relative_change}", "")

        # calculate constraint values
        for constraint in self.GetConstraints():
            algorithm_data["constraints"][constraint] = {
                "value": constraint.CalculateValue(),
                "standardized_value": constraint.CalculateStandardizedValue(),
                "is_active": constraint.IsActive()
            }
            self.__PrintInfo(1, constraint.GetResponseInfo(), "")

    def __ComputeResponseSensitivityForControlWrapper(self, controller: ControlTransformationTechnique, response_function: Union[ObjectiveResponseFunctionImplementor, ConstraintResponseFunctionImplementor], optimization_data: dict):
        control = controller.GetControl()
        for control_model_part in control.GetModelParts():
            optimization_data[control_model_part] = {}

            # calculate raw sensitivities
            raw_sensitivity_container = control.CreateContainerVariableDataHolder(control_model_part)
            response_function.CalculateStandardizedSensitivity(control.GetControlSensitivityVariable(), raw_sensitivity_container)

            # calculate transformed sensitivities
            transformed_sensitivities_container = raw_sensitivity_container.Clone()
            controller.TransformSensitivity(transformed_sensitivities_container)

            optimization_data[control_model_part][response_function] = {
                "raw_sensitivities": raw_sensitivity_container,
                "transformed_sensitivities": transformed_sensitivities_container
            }

    def __ComputeControlUpdatesForControlWrapper(self, controller: ControlTransformationTechnique):
        algorithm_data = self.optimization_info[self.GetName()]
        algorithm_data["controls"]["update"] = {}

        # get active constraints list
        active_constraints_data = [(constraint, constraint_data["standardized_value"]) for constraint, constraint_data in algorithm_data["constraints"].items() if constraint_data["is_active"]]

        # get active constraints values
        active_constraint_values = Kratos.Vector([data[1] for data in active_constraints_data])

        control = controller.GetControl()
        for control_model_part in control.GetModelParts():
            # get the objective sensitivities
            transformed_objective_sensitivities: ContainerVariableDataHolderUnion = algorithm_data["controls"][
                controller]["objectives"][control_model_part][self.objective]["transformed_sensitivities"]

            number_of_active_constraints = len(active_constraints_data)

            # get the active constraint sensitivities
            if number_of_active_constraints == 0:
                self.__PrintInfo(1, "No constraints active, use negative objective gradient as search direction.")
                # calculate control_update, search_direction, search_correction
                search_direction = transformed_objective_sensitivities * (-1.0)
                search_direction *= float(self.step_size / KratosOA.ContainerVariableDataHolderUtils.NormInf(search_direction))

                control_update = search_direction.Clone()

                search_correction = control.CreateContainerVariableDataHolder(control_model_part)
                search_correction.SetDataForContainerVariableToZero(control.GetControlSensitivityVariable())
            else:
                control_model_part_constraint_data = algorithm_data["controls"][controller]["constraints"][control_model_part]
                transformed_constraint_sensitivities = [control_model_part_constraint_data[data[0]]["transformed_sensitivities"] for data in active_constraints_data]

                # compute the projected search direction and correction
                ntn = Kratos.Matrix(number_of_active_constraints, number_of_active_constraints)
                for i in range(number_of_active_constraints):
                    for j in range(i, number_of_active_constraints):
                        ntn[i, j] = KratosOA.ContainerVariableDataHolderUtils.InnerProduct(transformed_constraint_sensitivities[i], transformed_constraint_sensitivities[j])
                        ntn[j, i] = ntn[i, j]

                # get the inverse of ntn
                ntn_inverse = Kratos.Matrix(number_of_active_constraints, number_of_active_constraints)

                if control_model_part.GetCommunicator().MyPID() == 0:
                    # since this is a small matrix, this is done serially
                    identity_matrix = Kratos.Matrix(number_of_active_constraints, number_of_active_constraints, 0.0)
                    for i in range(number_of_active_constraints):
                        identity_matrix[i, i] = 1.0

                    self.linear_solver.Solve(ntn, ntn_inverse, identity_matrix)

                # now broadcast the ntn_inverse to all ranks
                ntn_inverse = control_model_part.GetCommunicator().GetDataCommunicator().Broadcast(ntn_inverse, 0)

                # computing N^T\cdot\nabla{F}
                nt_nabla_f = Kratos.Vector(number_of_active_constraints)
                for i in range(number_of_active_constraints):
                    nt_nabla_f[i] = KratosOA.ContainerVariableDataHolderUtils.InnerProduct(transformed_objective_sensitivities, transformed_constraint_sensitivities[i])

                ntn_inv_nt_nabla_f = ntn_inverse * nt_nabla_f
                ntn_inv_constraint_values = ntn_inverse * active_constraint_values

                search_direction = -transformed_objective_sensitivities
                search_correction = control.CreateContainerVariableDataHolder(control_model_part)
                search_correction.SetDataForContainerVariableToZero(control.GetControlUpdateVariable())

                for i in range(number_of_active_constraints):
                    search_direction += transformed_constraint_sensitivities[i] * ntn_inv_nt_nabla_f[i]
                    search_correction -= transformed_constraint_sensitivities[i] * ntn_inv_constraint_values[i]

                search_direction_norm_inf = KratosOA.ContainerVariableDataHolderUtils.NormInf(search_direction)
                search_correction_norm_inf = KratosOA.ContainerVariableDataHolderUtils.NormInf(search_correction)
                if abs(search_correction_norm_inf) > 1e-12:
                    if search_direction_norm_inf < self.max_correction_share * self.step_size:
                        search_direction *= (self.step_size - search_correction_norm_inf) / search_direction_norm_inf
                    else:
                        search_correction *= self.max_correction_share * self.step_size / search_correction_norm_inf
                        search_direction *= (1 - self.max_correction_share) * self.step_size / search_direction_norm_inf
                else:
                    search_direction *= self.step_size / search_direction_norm_inf

                control_update = search_direction + search_correction

            # calculate transformed updates
            transformed_updates_container = control_update.Clone()
            controller.TransformUpdate(transformed_updates_container)

            # now set the control update
            controller.SetControlUpdate(transformed_updates_container.Clone())

            # now record computed values in optimization info for later visualization
            algorithm_data["controls"]["update"][controller] = {
                "search_direction": search_direction,
                "search_correction": search_correction,
                "raw_control_update": control_update,
                "transformed_control_update": transformed_updates_container
            }

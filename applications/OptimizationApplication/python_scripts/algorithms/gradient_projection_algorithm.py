from typing import Union

import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.LinearSolversApplication.dense_linear_solver_factory import ConstructSolver
from KratosMultiphysics.OptimizationApplication.optimization_info import OptimizationInfo
from KratosMultiphysics.OptimizationApplication.controls.control_wrapper import ControlWrapper
from KratosMultiphysics.OptimizationApplication.algorithms.algorithm import Algorithm
from KratosMultiphysics.OptimizationApplication.algorithms.algorithm_wrapper import ObjectiveResponseFunctionWrapper
from KratosMultiphysics.OptimizationApplication.algorithms.algorithm_wrapper import ConstraintResponseFunctionWrapper
from KratosMultiphysics.OptimizationApplication.utilities.container_data import ContainerData

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
        self.objective: ObjectiveResponseFunctionWrapper = self.GetObjectives()[0]

    def SolveSolutionStep(self) -> bool:
        # iniitalize optimization info for the algorithm
        self.optimization_info[self.GetName()] = {
            "objectives": {},
            "constraints": {},
            "controls": {}
        }

        # get algorithm specific data from optimization info
        algorithm_data = self.optimization_info[self.GetName()]

        # calculate objective value
        algorithm_data["objectives"][self.objective] = {
            "value": self.objective.CalculateValue(),
            "standardized_value": self.objective.CalculateStandardizedValue()
        }
        self.__PrintInfo(1, self.objective.GetResponseInfo(), "")

        # calculate constraint values
        constraint: ConstraintResponseFunctionWrapper
        for constraint in self.GetConstraints():
            algorithm_data["constraints"][constraint] = {
                "value": constraint.CalculateValue(),
                "standardized_value": constraint.CalculateStandardizedValue(),
                "is_active": constraint.IsActive()
            }
            self.__PrintInfo(1, constraint.GetResponseInfo(), "")

        # calculate response gradients
        control_wrapper: ControlWrapper
        for control_wrapper in self.GetControllers():
            self.__PrintInfo(1, f"Computing sensitivities for {control_wrapper.GetName()} control...")

            # create control_wrapper optimization info
            algorithm_data["controls"][control_wrapper] = {
                "objectives" : {},
                "constraints": {}
            }

            control_data = algorithm_data["controls"][control_wrapper]

            # compute objective sensitivities
            self.__ComputeResponseSensitivityForControlWrapper(control_wrapper, self.objective, control_data["objectives"])

            # compute constraint sensitivities
            for constraint, constraint_data in algorithm_data["constraints"].items():
                if constraint_data["is_active"]:
                    self.__ComputeResponseSensitivityForControlWrapper(control_wrapper, constraint, control_data["constraints"])

            # compute control update
            self.__ComputeControlUpdatesForControlWrapper(control_wrapper)

            self.__PrintInfo(1, f"Computed sensitivities for {control_wrapper.GetName()} control.")

    def IsConverged(self) -> bool:
        if self.optimization_info["step"] > 1:
            # check for objective convergence
            is_converged = abs(self.optimization_info[self.GetName()]["objectives"][self.objective]["value"] / self.optimization_info.GetSolutionStepData(1)[self.GetName()]["objectives"][self.objective]["value"] - 1.0) < self.relative_tolerance

            # check for constraint convergence
            for constraint_data in self.optimization_info[self.GetName()]["constraints"].values():
                is_converged = is_converged and not constraint_data["is_active"]

            return is_converged
        else:
            return False

    def __PrintInfo(self, required_echo_level: int, message: str, title = "GradientProjectionAlgorithm"):
        if self.echo_level >= required_echo_level:
            Kratos.Logger.PrintInfo(title, message)

    def __ComputeResponseSensitivityForControlWrapper(self, control_wrapper: ControlWrapper, response_function: Union[ObjectiveResponseFunctionWrapper, ConstraintResponseFunctionWrapper], optimization_data: dict):
        control = control_wrapper.GetControl()
        control_model_part: Kratos.ModelPart
        for control_model_part in control.GetModelParts():
            optimization_data[control_model_part] = {}

            # calculate raw sensitivities
            raw_sensitivity_container = ContainerData(control_model_part, control.GetContainerType())
            response_function.CalculateStandardizedSensitivity(control.GetControlSensitivityVariable(), raw_sensitivity_container)

            # calculate modified sensitivities
            modified_sensitivities_container = raw_sensitivity_container.Clone()
            control_wrapper.ModifySensitivities(modified_sensitivities_container)

            optimization_data[control_model_part][response_function] = {
                "raw_sensitivities": raw_sensitivity_container,
                "modified_sensitivities": modified_sensitivities_container
            }

    def __ComputeControlUpdatesForControlWrapper(self, control_wrapper: ControlWrapper):
        algorithm_data = self.optimization_info[self.GetName()]
        algorithm_data["controls"]["update"] = {}

        # get active constraints list
        active_constraints_data = [(constraint, constraint_data["standardized_value"]) for constraint, constraint_data in algorithm_data["constraints"].items() if constraint_data["is_active"]]

        # get active constraints values
        active_constraint_values_vector = Kratos.Vector([data[1] for data in active_constraints_data])

        control = control_wrapper.GetControl()
        control_model_part: Kratos.ModelPart
        for control_model_part in control.GetModelParts():
            # get the objective sensitivities
            modified_objective_sensitivities: ContainerData = algorithm_data["controls"][control_wrapper]["objectives"][control_model_part][self.objective]["modified_sensitivities"]

            # get the active constraint sensitivities
            if len(active_constraints_data) == 0:
                self.__PrintInfo(1, "No constraints active, use negative objective gradient as search direction.")
                search_direction = modified_objective_sensitivities * (-1.0)
                search_direction *= self.step_size / search_direction.NormInf()
                control_update = search_direction.Clone()
                search_correction = ContainerData(control_model_part, control.GetContainerType())
                search_correction.SetData(Kratos.Vector(control_update.GetData().Size(), 0.0))
            else:
                control_model_part_constraint_data = algorithm_data["controls"][control_wrapper]["constraints"][control_model_part]
                modified_constraint_sensitivities = [control_model_part_constraint_data[data[0]]["modified_sensitivities"].GetData() for data in active_constraints_data]

                # compute the projected search direction and correction
                search_direction = ContainerData(control_model_part, control.GetContainerType())
                search_correction = ContainerData(control_model_part, control.GetContainerType())
                KratosOA.GradientProjectionSolverUtils.CalculateProjectedSearchDirectionAndCorrection(
                    search_direction.GetData(),
                    search_correction.GetData(),
                    self.linear_solver,
                    active_constraint_values_vector,
                    modified_objective_sensitivities.GetData(),
                    modified_constraint_sensitivities)

                # calculate control update
                control_update = ContainerData(control_model_part, control.GetContainerType())
                KratosOA.GradientProjectionSolverUtils.CalculateControlUpdate(
                    control_update.GetData(),
                    search_direction.GetData(),
                    search_correction.GetData(),
                    self.step_size,
                    self.max_correction_share)

            # now set the control update
            control_wrapper.SetControlUpdate(control_update.Clone())

            # now record computed values in optimization info for later visualization
            algorithm_data["controls"]["update"][control_wrapper] = {
                "search_direction": search_direction,
                "search_correction": search_correction,
                "control_update": control_update
            }



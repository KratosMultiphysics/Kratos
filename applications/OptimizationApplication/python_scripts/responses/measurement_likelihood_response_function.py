from __future__ import annotations
import json

import numpy as np

import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.execution_policies.execution_policy_decorator import ExecutionPolicyDecorator
from KratosMultiphysics.OptimizationApplication.responses.response_function import ResponseFunction
from KratosMultiphysics.OptimizationApplication.responses.response_function import SupportedSensitivityFieldVariableTypes
from KratosMultiphysics.OptimizationApplication.utilities.model_part_utilities import ModelPartUtilities

from KratosMultiphysics.StructuralMechanicsApplication import structural_mechanics_analysis
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication


def Factory(model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem) -> ResponseFunction:
    if not parameters.Has("name"):
        raise RuntimeError(f"MeasurementLikelihoodResponseFunction instantiation requires a \"name\" in parameters [ parameters = {parameters}].")
    if not parameters.Has("settings"):
        raise RuntimeError(f"MeasurementLikelihoodResponseFunction instantiation requires a \"settings\" in parameters [ parameters = {parameters}].")

    return MeasurementLikelihoodResponseFunction(parameters["name"].GetString(), model, parameters["settings"], optimization_problem)


class MeasurementLikelihoodResponseFunction(ResponseFunction):
    def __init__(self, name: str, model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
        super().__init__(name)

        default_settings = Kratos.Parameters("""{
            "primal_analysis_name"           : "",
            "measurement_data_file"          : "MeasurementData.json",
            "measurement_standard_deviation"  : 1,
            "perturbation_size"              : 1e-8,
            "adjoint_analysis_settings_file_name": "adjoint_parameters.json",
            "evaluated_model_part_names"     : [
                "PLEASE_PROVIDE_A_MODEL_PART_NAME"
            ]
        }""")
        parameters.ValidateAndAssignDefaults(default_settings)

        self.model_part_names = parameters["evaluated_model_part_names"].GetStringArray()
        self.perturbation_size = parameters["perturbation_size"].GetDouble()

        self.measurement_data_file = parameters["measurement_data_file"].GetString()
        self.measurement_std = parameters["measurement_standard_deviation"].GetDouble()

        self.adjoint_analysis_file_name = parameters["adjoint_analysis_settings_file_name"].GetString()

        self.model = model
        self.model_part: Kratos.ModelPart = None
        self.adjoint_model = Kratos.Model()

        self.primal_analysis_execution_policy_decorator: ExecutionPolicyDecorator = optimization_problem.GetExecutionPolicy(parameters["primal_analysis_name"].GetString())

        if len(self.model_part_names) == 0:
            raise RuntimeError("No model parts were provided for MeasurementLikelihoodResponseFunction.")

    def GetImplementedPhysicalKratosVariables(self) -> list[SupportedSensitivityFieldVariableTypes]:
        return [Kratos.YOUNG_MODULUS]

    def Initialize(self) -> None:
        model_parts_list = [self.model[model_part_name] for model_part_name in self.model_part_names]
        root_model_part = model_parts_list[0].GetRootModelPart()
        _, self.model_part = ModelPartUtilities.UnionModelParts(root_model_part, model_parts_list, False)

        file = open(self.measurement_data_file)
        self.measurement_data = json.load(file)

        self.sensor_positions_model_part = self.model.CreateModelPart("sensor_positions")

        for sensor in self.measurement_data["load_cases"][0]["sensors_infos"]:

            point = Kratos.Point(sensor["position_of_mesh_node"])
            search_configuration = Kratos.Configuration.Initial
            search_tolerance = 1e-6

            point_locator = Kratos.BruteForcePointLocator(self.model_part)

            found_node_id = point_locator.FindNode(
                point, search_configuration, search_tolerance
            )

            sensor["mesh_node_id"] = found_node_id
            # Add node to sensor_positions_model_part
            self.sensor_positions_model_part.AddNode(self.model_part.GetNode(found_node_id))

        # Update measurement data file with newly added infos
        with open(self.measurement_data_file, "w") as outfile:
            json.dump(self.measurement_data, outfile)

    def Check(self) -> None:
        pass

    def Finalize(self) -> None:
        pass

    def GetEvaluatedModelPart(self) -> Kratos.ModelPart:
        if self.model_part is None:
            raise RuntimeError("Please call MeasurementLikelihoodResponseFunction::Initialize first.")
        return self.model_part

    def GetAnalysisModelPart(self) -> Kratos.ModelPart:
        return self.primal_analysis_execution_policy_decorator.GetExecutionPolicy().GetAnalysisModelPart()

    def CalculateValue(self) -> float:
        # LL = 0.5*(1/covariance) * error.T @ error  # Without constant terms, and for diag. covariance matrix
        # LL = 1/2 * 1/cov * error**2

        objective_value = 0
        for sensor in self.measurement_data["load_cases"][0]["sensors_infos"]:
            measured_displacement = sensor["measured_value"]
            measured_direction = Kratos.Array3(sensor["measurement_direction_normal"])

            mesh_node = self.model_part.GetNode(sensor["mesh_node_id"])
            node_displacement = mesh_node.GetSolutionStepValue(Kratos.DISPLACEMENT)
            in_measurement_direction_projected_vector = (measured_direction[0]*node_displacement[0])+(measured_direction[1]*node_displacement[1])+(measured_direction[2]*node_displacement[2])

            objective_value += (measured_displacement-in_measurement_direction_projected_vector)**2

        objective_value *= 0.5 * 1/self.measurement_std

        return objective_value

    def CalculateGradient(self, physical_variable_collective_expressions: dict[SupportedSensitivityFieldVariableTypes, KratosOA.ContainerExpression.CollectiveExpressions]) -> None:
        # first calculate the gradients
        merged_model_part_map = ModelPartUtilities.GetMergedMap(self.model_part, physical_variable_collective_expressions, False)
        intersected_model_part_map = ModelPartUtilities.GetIntersectedMap(self.GetAnalysisModelPart(), merged_model_part_map, True)

        ###################################################################

        for variable, collective_expression in physical_variable_collective_expressions.items():

            with open(self.adjoint_analysis_file_name, 'r') as parameter_file:
                adjoint_parameters = Kratos.Parameters(parameter_file.read())

            adjoint_parameters["solver_settings"]["sensitivity_settings"]["element_data_value_sensitivity_variables"].SetStringArray([variable.Name()])

            adjoint_analysis = structural_mechanics_analysis.StructuralMechanicsAnalysis(self.adjoint_model, adjoint_parameters)
            adjoint_analysis.Run()

            Kratos.VariableUtils().CopyModelPartElementalVar(KratosOA.YOUNG_MODULUS_SENSITIVITY, self.adjoint_model.GetModelPart("Structure"), self.model["Structure"])

            # now fill the collective expressions

            for expression in collective_expression.GetContainerExpressions():
                if isinstance(expression, KratosOA.ContainerExpression.ElementPropertiesExpression):
                    element_data_expression = Kratos.ContainerExpression.ElementNonHistoricalExpression(expression.GetModelPart())
                    element_data_expression.Read(StructuralMechanicsApplication.YOUNG_MODULUS_SENSITIVITY)
                    expression.CopyFrom(element_data_expression)
                else:
                    expression.Read(Kratos.KratosGlobals.GetVariable(variable.Name() + "_SENSITIVITY"))

        #####################################################################

        # Likelihood Gradient
        # covariance = 0.001
        # S = -(-1/covariance) * b@A
        # Gradient = error @ d_disp/d_Youngs_modulus

        # TODO remove after testing
        # KratosOA.ResponseUtils.LinearStrainEnergyResponseUtils.CalculateGradient(
        #     list(merged_model_part_map.keys()),
        #     list(merged_model_part_map.values()),
        #     list(intersected_model_part_map.values()),
        #     self.perturbation_size)

        # # now fill the collective expressions
        # for variable, collective_expression in physical_variable_collective_expressions.items():
        #     collective_expression.Read(Kratos.KratosGlobals.GetVariable(variable.Name() + "_SENSITIVITY"))

    def __str__(self) -> str:
        return f"Response [type = {self.__class__.__name__}, name = {self.GetName()}, model part name = {self.model_part.FullName()}]"

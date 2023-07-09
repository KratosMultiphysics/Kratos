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
        raise RuntimeError(f"MeasurementResidualResponseFunction instantiation requires a \"name\" in parameters [ parameters = {parameters}].")
    if not parameters.Has("settings"):
        raise RuntimeError(f"MeasurementResidualResponseFunction instantiation requires a \"settings\" in parameters [ parameters = {parameters}].")

    return MeasurementResidualResponseFunction(parameters["name"].GetString(), model, parameters["settings"], optimization_problem)


class MeasurementResidualResponseFunction(ResponseFunction):
    def __init__(self, name: str, model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
        super().__init__(name)

        default_settings = Kratos.Parameters("""{
            "primal_analysis_name"           : "",
            "measurement_data_file"          : "MeasurementData.json",
            "measurement_standard_deviation"  : 1.0,
            "perturbation_size"              : 1e-8,
            "adjoint_analysis_settings_file_name": "adjoint_parameters.json",
            "evaluated_model_part_names"     : [
                "PLEASE_PROVIDE_A_MODEL_PART_NAME"
            ]
        }""")
        parameters.ValidateAndAssignDefaults(default_settings)

        self.model = model
        self.primal_model_part: Kratos.ModelPart = None

        evaluated_model_parts = [model[model_part_name] for model_part_name in parameters["evaluated_model_part_names"].GetStringArray()]
        if len(evaluated_model_parts) == 0:
            raise RuntimeError(f"No model parts were provided for LinearStrainEnergyResponseFunction. [ response name = \"{self.GetName()}\"]")
        self.primal_model_part = ModelPartUtilities.GetOperatingModelPart(ModelPartUtilities.OperationType.UNION, f"response_{self.GetName()}", evaluated_model_parts, False)

        self.perturbation_size = parameters["perturbation_size"].GetDouble()

        self.measurement_data_file = parameters["measurement_data_file"].GetString()
        self.measurement_std = parameters["measurement_standard_deviation"].GetDouble()

        self.adjoint_analysis_file_name = parameters["adjoint_analysis_settings_file_name"].GetString()

        self.primal_analysis_execution_policy_decorator: ExecutionPolicyDecorator = optimization_problem.GetExecutionPolicy(parameters["primal_analysis_name"].GetString())

    def GetImplementedPhysicalKratosVariables(self) -> list[SupportedSensitivityFieldVariableTypes]:
        return [Kratos.YOUNG_MODULUS]

    def Initialize(self) -> None:
        if not KratosOA.ModelPartUtils.CheckModelPartStatus(self.primal_model_part, "element_specific_properties_created"):
            KratosOA.OptimizationUtils.CreateEntitySpecificPropertiesForContainer(self.primal_model_part, self.primal_model_part.Elements)
            KratosOA.ModelPartUtils.LogModelPartStatus(self.primal_model_part, "element_specific_properties_created")

        with open(self.adjoint_analysis_file_name, 'r') as parameter_file:
            self.adjoint_parameters = Kratos.Parameters(parameter_file.read())

        file = open(self.measurement_data_file)
        self.measurement_data = json.load(file)

        ModelPartUtilities.ExecuteOperationOnModelPart(self.primal_model_part)

        # self.sensor_positions_model_part = self.model.CreateModelPart("sensor_positions")

        # for sensor in self.measurement_data["load_cases"][0]["sensors_infos"]:

        #     point = Kratos.Point(sensor["position_of_mesh_node"])
        #     search_configuration = Kratos.Configuration.Initial
        #     search_tolerance = 1e-6

        #     point_locator = Kratos.BruteForcePointLocator(self.model_part)

        #     found_node_id = point_locator.FindNode(
        #         point, search_configuration, search_tolerance
        #     )

        #     sensor["mesh_node_id"] = found_node_id
        #     # Add node to sensor_positions_model_part
        #     self.sensor_positions_model_part.AddNode(self.model_part.GetNode(found_node_id))

        # Update measurement data file with newly added infos
        # with open(self.measurement_data_file, "w") as outfile:
        #     json.dump(self.measurement_data, outfile)

    def Check(self) -> None:
        KratosOA.PropertiesVariableExpressionIO.Check(Kratos.Expression.ElementExpression(self.primal_model_part), Kratos.YOUNG_MODULUS)

    def Finalize(self) -> None:
        pass

    def GetEvaluatedModelPart(self) -> Kratos.ModelPart:
        if self.primal_model_part is None:
            raise RuntimeError("Please call MeasurementResidualResponseFunction::Initialize first.")
        return self.primal_model_part

    def GetAnalysisModelPart(self) -> Kratos.ModelPart:
        return self.primal_analysis_execution_policy_decorator.GetAnalysisModelPart()

    def CalculateValue(self) -> float:
        self.primal_analysis_execution_policy_decorator.Execute()

        objective_value = 0
        for sensor in self.measurement_data["load_cases"][0]["sensors_infos"]:
            measured_displacement = sensor["measured_value"]
            measured_direction = Kratos.Array3(sensor["measurement_direction_normal"])

            mesh_node = self.primal_model_part.GetNode(sensor["mesh_node_id"])
            node_displacement = mesh_node.GetSolutionStepValue(Kratos.DISPLACEMENT)
            in_measurement_direction_projected_vector = (measured_direction[0]*node_displacement[0])+(measured_direction[1]*node_displacement[1])+(measured_direction[2]*node_displacement[2])

            objective_value += (measured_displacement-in_measurement_direction_projected_vector)**2

        objective_value *= 0.5 * 1/self.measurement_std

        return objective_value

    def CalculateGradient(self, physical_variable_collective_expressions: "dict[SupportedSensitivityFieldVariableTypes, KratosOA.CollectiveExpression]") -> None:
        # first merge all the model parts
        merged_model_part_map = ModelPartUtilities.GetMergedMap(physical_variable_collective_expressions, False)

        # now get the intersected model parts
        intersected_model_part_map = ModelPartUtilities.GetIntersectedMap(self.primal_model_part, merged_model_part_map, True)

        for variable, collective_expression in physical_variable_collective_expressions.items():

            self.adjoint_parameters["solver_settings"]["sensitivity_settings"]["element_data_value_sensitivity_variables"].SetStringArray([variable.Name()])

            self.adjoint_model = Kratos.Model()
            adjoint_analysis = structural_mechanics_analysis.StructuralMechanicsAnalysis(self.adjoint_model, self.adjoint_parameters)
            adjoint_model_part: Kratos.ModelPart = adjoint_analysis._GetSolver().GetComputingModelPart()

            adjoint_analysis.Initialize()

            # create element specific properties in adjoint model part
            KratosOA.OptimizationUtils.CreateEntitySpecificPropertiesForContainer(adjoint_model_part, adjoint_model_part.Elements)

            # read primal E
            primal_youngs_modulus = Kratos.Expression.ElementExpression(self.primal_model_part)
            KratosOA.PropertiesVariableExpressionIO.Read(primal_youngs_modulus, Kratos.YOUNG_MODULUS)

            # assign primal E to adjoint
            adjoint_young_modulus = Kratos.Expression.ElementExpression(adjoint_model_part)
            adjoint_young_modulus.SetExpression(primal_youngs_modulus.GetExpression())
            KratosOA.PropertiesVariableExpressionIO.Write(adjoint_young_modulus, Kratos.YOUNG_MODULUS)

            adjoint_analysis.RunSolutionLoop()
            adjoint_analysis.Finalize()

            # now fill the collective expressions
            for expression in collective_expression.GetContainerExpressions():
                if isinstance(expression, Kratos.Expression.ElementExpression):
                    adjoint_young_modulus_sensitivity = Kratos.Expression.ElementExpression(adjoint_model_part)
                    Kratos.Expression.VariableExpressionIO.Read(adjoint_young_modulus_sensitivity, KratosOA.YOUNG_MODULUS_SENSITIVITY)
                    expression.SetExpression(adjoint_young_modulus_sensitivity.GetExpression())
                else:
                    KratosOA.PropertiesVariableExpressionIO.Read(expression, Kratos.KratosGlobals.GetVariable(variable.Name() + "_SENSITIVITY"))

            # Calculate via finite differencing
            # for expression in collective_expression.GetContainerExpressions():
            #     gradient = self.CalculateGradientWithFiniteDifferencing(physical_variable_collective_expressions)
            #     np_gradient = np.array(gradient)

            #     adjoint_young_modulus_sensitivity = Kratos.Expression.ElementExpression(expression.GetModelPart())
            #     Kratos.Expression.CArrayExpressionIO.Read(adjoint_young_modulus_sensitivity, np_gradient)
            #     expression.SetExpression(adjoint_young_modulus_sensitivity.GetExpression())

    def CalculateGradientWithFiniteDifferencing(self, physical_variable_collective_expressions: "dict[SupportedSensitivityFieldVariableTypes, KratosOA.CollectiveExpression]") -> "list(float)":
        # for element in self.model_part.Elements:
        #     E = element.Properties[Kratos.YOUNG_MODULUS]
        #     element.Properties[Kratos.YOUNG_MODULUS] = E

        # # checking
        # check_expression = Kratos.Expression.ElementExpression(expression.GetModelPart())
        # KratosOA.PropertiesVariableExpressionIO.Read(check_expression, StructuralMechanicsApplication.YOUNG_MODULUS_SENSITIVITY)
        # numpy_array = check_expression.Evaluate()
        # Kratos.Expression.CArrayExpressionIO.Read(check_expression, numpy_array)
        # Kratos.Expression.CArrayExpressionIO.Move(check_expression, numpy_array)
        # Kratos.Expression.CArrayExpressionIO.Write(check_expression, numpy_array)

        # first merge all the model parts
        merged_model_part_map = ModelPartUtilities.GetMergedMap(physical_variable_collective_expressions, False)

        # now get the intersected model parts
        intersected_model_part_map = ModelPartUtilities.GetIntersectedMap(self.primal_model_part, merged_model_part_map, True)

        model_part = self.primal_analysis_execution_policy_decorator.GetAnalysisModelPart()

        if not KratosOA.ModelPartUtils.CheckModelPartStatus(model_part, "element_specific_properties_created"):
            KratosOA.OptimizationUtils.CreateEntitySpecificPropertiesForContainer(model_part, model_part.Elements)
            KratosOA.ModelPartUtils.LogModelPartStatus(model_part, "element_specific_properties_created")

        perturbation = model_part.GetElement(1).Properties[Kratos.YOUNG_MODULUS] * 1e-6

        for variable, collective_expression in physical_variable_collective_expressions.items():

            reference_value = self.CalculateValue()

            gradient = []

            for i, element in enumerate(model_part.Elements):

                old_stiffness = element.Properties[Kratos.YOUNG_MODULUS]
                element.Properties[Kratos.YOUNG_MODULUS] += perturbation

                objective_value = self.CalculateValue()

                gradient.append((reference_value-objective_value)/perturbation)
                element.Properties[Kratos.YOUNG_MODULUS] = old_stiffness

            return gradient

    def __str__(self) -> str:
        return f"Response [type = {self.__class__.__name__}, name = {self.GetName()}, model part name = {self.primal_model_part.FullName()}]"

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
            "measurement_data_file"          : "MeasurementData.json",
            "measurement_standard_deviation"  : 1.0,
            "perturbation_size"              : 1e-8,
            "analysis_sets": [{
                "primal_analysis_name": "",
                "adjoint_analysis_settings_file_name": "",
                "load_case_in_measurement_file": -1
            }],
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

        self.analysis_sets = parameters["analysis_sets"]

        if self.analysis_sets.size() == 0:
            raise RuntimeError("MeasurementResidualResponseFunction:: Please provide analysis sets for the response function.")

        self.primal_analysis_execution_policy_decorators: "List(ExecutionPolicyDecorator)" = []
        for analysis_set in self.analysis_sets:
            self.primal_analysis_execution_policy_decorators.append(optimization_problem.GetExecutionPolicy(analysis_set["primal_analysis_name"].GetString()))

    def GetImplementedPhysicalKratosVariables(self) -> list[SupportedSensitivityFieldVariableTypes]:
        return [Kratos.YOUNG_MODULUS]

    def Initialize(self) -> None:
        if not KratosOA.ModelPartUtils.CheckModelPartStatus(self.primal_model_part, "element_specific_properties_created"):
            KratosOA.OptimizationUtils.CreateEntitySpecificPropertiesForContainer(self.primal_model_part, self.primal_model_part.Elements)
            KratosOA.ModelPartUtils.LogModelPartStatus(self.primal_model_part, "element_specific_properties_created")

        file = open(self.measurement_data_file)
        self.measurement_data = json.load(file)

        ModelPartUtilities.ExecuteOperationOnModelPart(self.primal_model_part)

    def Check(self) -> None:
        KratosOA.PropertiesVariableExpressionIO.Check(Kratos.Expression.ElementExpression(self.primal_model_part), Kratos.YOUNG_MODULUS)

    def Finalize(self) -> None:
        pass

    def GetEvaluatedModelPart(self) -> Kratos.ModelPart:
        if self.primal_model_part is None:
            raise RuntimeError("Please call MeasurementResidualResponseFunction::Initialize first.")
        return self.primal_model_part

    def GetAnalysisModelPart(self) -> Kratos.ModelPart:
        return self.primal_analysis_execution_policy_decorators[0].GetAnalysisModelPart()

    def CalculateValue(self) -> float:

        objective_value = 0

        for i, execution_policy in enumerate(self.primal_analysis_execution_policy_decorators):
            execution_policy.Execute()

            load_case = self.analysis_sets[i]["load_case_in_measurement_file"].GetInt()

            for sensor in self.measurement_data["load_cases"][load_case]["sensors_infos"]:
                measured_displacement = sensor["measured_value"]
                measured_direction = Kratos.Array3(sensor["measurement_direction_normal"])

                mesh_node = self.primal_model_part.GetNode(sensor["mesh_node_id"])

                node_displacement = mesh_node.GetSolutionStepValue(Kratos.KratosGlobals.GetVariable(sensor["type_of_sensor"]))
                in_measurement_direction_projected_vector = (measured_direction[0]*node_displacement[0])+(measured_direction[1]*node_displacement[1])+(measured_direction[2]*node_displacement[2])

                objective_value += (measured_displacement-in_measurement_direction_projected_vector)**2

        objective_value *= 0.5 * 1/self.measurement_std

        return objective_value

    def transform_sensitivities_to_continuous_space(self, expression: Kratos.Expression) -> None:
        sensitivities = expression.Evaluate()

        for i, element in enumerate(expression.GetModelPart().GetElements()):
            occupied_space = -1
            geometry = element.GetGeometry()
            try:
                occupied_space = geometry.Area()
            except:
                pass

            try:
                occupied_space = geometry.Volume()
            except:
                pass

            if occupied_space == -1:
                raise RuntimeError(
                    f"MeasurementResidualResponseFunction:: Error in normalize_sensitivities_to_element_area. The Geometry object {geometry} of the provided element {element} can not calculate an area or volume.")

            sensitivities[i] /= occupied_space

        Kratos.Expression.CArrayExpressionIO.Read(expression, sensitivities)

    def CalculateGradient(self, physical_variable_collective_expressions: "dict[SupportedSensitivityFieldVariableTypes, KratosOA.CollectiveExpression]") -> None:
        # first merge all the model parts
        merged_model_part_map = ModelPartUtilities.GetMergedMap(physical_variable_collective_expressions, False)

        # now get the intersected model parts
        intersected_model_part_map = ModelPartUtilities.GetIntersectedMap(self.primal_model_part, merged_model_part_map, True)

        print("MeasurementResidualResponseFunction:: Start gradient calculation")
        for variable, collective_expression in physical_variable_collective_expressions.items():

            if variable.Name() != "YOUNG_MODULUS":
                raise RuntimeError("MeasurementResidualResponseFunction:: CalculateGradient currently only supports 'YOUNG_MODULUS' as property")

            summed_gradients: np.ndarray = np.ndarray([])
            current_sensitivities: np.ndarray = np.ndarray([])

            for analysis_set in self.analysis_sets:

                with open(analysis_set["adjoint_analysis_settings_file_name"].GetString(), 'r') as parameter_file:
                    self.adjoint_parameters = Kratos.Parameters(parameter_file.read())

                self.adjoint_parameters["solver_settings"]["sensitivity_settings"]["element_data_value_sensitivity_variables"].SetStringArray([variable.Name()])
                self.adjoint_parameters["solver_settings"]["response_function_settings"]["measurement_load_case_to_use"].SetInt(analysis_set["load_case_in_measurement_file"].GetInt())

                adjoint_analysis = Kratos.Model()
                adjoint_analysis = structural_mechanics_analysis.StructuralMechanicsAnalysis(adjoint_analysis, self.adjoint_parameters)
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

                adjoint_young_modulus_sensitivity = Kratos.Expression.ElementExpression(adjoint_model_part)
                Kratos.Expression.VariableExpressionIO.Read(adjoint_young_modulus_sensitivity, KratosOA.YOUNG_MODULUS_SENSITIVITY)

                current_sensitivities = adjoint_young_modulus_sensitivity.Evaluate()
                if summed_gradients.size == 1:
                    summed_gradients = current_sensitivities
                else:
                    summed_gradients += current_sensitivities

            # now fill the collective expressions
            summed_gradients /= self.analysis_sets.size()
            for expression in collective_expression.GetContainerExpressions():
                Kratos.Expression.CArrayExpressionIO.Read(expression, summed_gradients)
                self.transform_sensitivities_to_continuous_space(expression)

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

        model_part = self.GetAnalysisModelPart()

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

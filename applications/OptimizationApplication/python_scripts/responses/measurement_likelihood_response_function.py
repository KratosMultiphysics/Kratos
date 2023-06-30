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

        self.model = model
        self.model_part: Kratos.ModelPart = None

        evaluated_model_parts = [model[model_part_name] for model_part_name in parameters["evaluated_model_part_names"].GetStringArray()]
        if len(evaluated_model_parts) == 0:
            raise RuntimeError(f"No model parts were provided for LinearStrainEnergyResponseFunction. [ response name = \"{self.GetName()}\"]")
        self.model_part = ModelPartUtilities.GetOperatingModelPart(ModelPartUtilities.OperationType.UNION, f"response_{self.GetName()}", evaluated_model_parts, False)

        self.perturbation_size = parameters["perturbation_size"].GetDouble()

        self.measurement_data_file = parameters["measurement_data_file"].GetString()
        self.measurement_std = parameters["measurement_standard_deviation"].GetDouble()

        self.adjoint_analysis_file_name = parameters["adjoint_analysis_settings_file_name"].GetString()

        self.primal_analysis_execution_policy_decorator: ExecutionPolicyDecorator = optimization_problem.GetExecutionPolicy(parameters["primal_analysis_name"].GetString())

    def GetImplementedPhysicalKratosVariables(self) -> list[SupportedSensitivityFieldVariableTypes]:
        return [Kratos.YOUNG_MODULUS]

    def Initialize(self) -> None:
        if not KratosOA.ModelPartUtils.CheckModelPartStatus(self.model_part, "element_specific_properties_created"):
            KratosOA.OptimizationUtils.CreateEntitySpecificPropertiesForContainer(self.model_part, self.model_part.Elements)
            KratosOA.ModelPartUtils.LogModelPartStatus(self.model_part, "element_specific_properties_created")

        with open(self.adjoint_analysis_file_name, 'r') as parameter_file:
            self.adjoint_parameters = Kratos.Parameters(parameter_file.read())

        file = open(self.measurement_data_file)
        self.measurement_data = json.load(file)

        ModelPartUtilities.ExecuteOperationOnModelPart(self.model_part)

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
        KratosOA.PropertiesVariableExpressionIO.Check(Kratos.Expression.ElementExpression(self.model_part), Kratos.YOUNG_MODULUS)

    def Finalize(self) -> None:
        pass

    def GetEvaluatedModelPart(self) -> Kratos.ModelPart:
        if self.model_part is None:
            raise RuntimeError("Please call MeasurementLikelihoodResponseFunction::Initialize first.")
        return self.model_part

    def GetAnalysisModelPart(self) -> Kratos.ModelPart:
        return self.primal_analysis_execution_policy_decorator.GetAnalysisModelPart()

    def CalculateValue(self) -> float:

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

    def CalculateGradient(self, physical_variable_collective_expressions: "dict[SupportedSensitivityFieldVariableTypes, KratosOA.CollectiveExpression]") -> None:
        # first merge all the model parts
        merged_model_part_map = ModelPartUtilities.GetMergedMap(physical_variable_collective_expressions, False)

        # now get the intersected model parts
        intersected_model_part_map = ModelPartUtilities.GetIntersectedMap(self.model_part, merged_model_part_map, True)

        for variable, collective_expression in physical_variable_collective_expressions.items():

            self.adjoint_parameters["solver_settings"]["sensitivity_settings"]["element_data_value_sensitivity_variables"].SetStringArray([variable.Name()])

            self.adjoint_model = Kratos.Model()
            adjoint_analysis = structural_mechanics_analysis.StructuralMechanicsAnalysis(self.adjoint_model, self.adjoint_parameters)
            adjoint_model_part = adjoint_analysis._GetSolver().GetComputingModelPart()

            adjoint_analysis.Initialize()

            # create element specific properties in adjoint model part
            KratosOA.OptimizationUtils.CreateEntitySpecificPropertiesForContainer(adjoint_model_part, adjoint_model_part.Elements)

            # read primal E
            primal_youngs_modulus = Kratos.Expression.ElementExpression(self.model_part)
            KratosOA.PropertiesVariableExpressionIO.Read(primal_youngs_modulus, Kratos.YOUNG_MODULUS)
            # print("primal_youngs_modulus", primal_youngs_modulus.Evaluate())

            # assign primal E to adjoint
            adjoint_young_modulus = Kratos.Expression.ElementExpression(adjoint_model_part)
            adjoint_young_modulus.SetExpression(primal_youngs_modulus.GetExpression())
            KratosOA.PropertiesVariableExpressionIO.Write(adjoint_young_modulus, Kratos.YOUNG_MODULUS)

            # def modify_parameters_in_analysis_stage():
            #     for i, element in enumerate(self.model.GetModelPart(root_model_part_name).Elements):
            #         element.Properties[Kratos.YOUNG_MODULUS] = self.GetAnalysisModelPart().GetElement(element.Id).Properties[Kratos.YOUNG_MODULUS]

            # adjoint_analysis.ModifyInitialProperties = modify_parameters_in_analysis_stage  # override responding function, it's a hack!!!!!

            adjoint_analysis.RunSolutionLoop()
            adjoint_analysis.Finalize()

            # for element in self.adjoint_model[self.adjoint_model.GetModelPartNames()[0]].Elements:
            # print("getvalue", element.GetValue(KratosOA.YOUNG_MODULUS_SENSITIVITY))
            # print("propertiesvalue", element.Properties[KratosOA.YOUNG_MODULUS_SENSITIVITY])
            # print("propertiesvalue", element.Properties[Kratos.YOUNG_MODULUS])

            # raise RuntimeError(1)

            # Kratos.VariableUtils().CopyModelPartElementalVar(KratosOA.YOUNG_MODULUS_SENSITIVITY, self.adjoint_model[root_model_part_name], self.model[root_model_part_name])

            # now fill the collective expressions

            for expression in collective_expression.GetContainerExpressions():
                if isinstance(expression, Kratos.Expression.ElementExpression):
                    adjoint_young_modulus_sensitivity = Kratos.Expression.ElementExpression(adjoint_model_part)
                    Kratos.Expression.VariableExpressionIO.Read(adjoint_young_modulus_sensitivity, KratosOA.YOUNG_MODULUS_SENSITIVITY)
                    expression.SetExpression(adjoint_young_modulus_sensitivity.GetExpression())
                    # print("expression", expression.Evaluate())
                    # raise RuntimeError(1)
                else:
                    KratosOA.PropertiesVariableExpressionIO.Read(expression, Kratos.KratosGlobals.GetVariable(variable.Name() + "_SENSITIVITY"))

            # print("adjoint_young_modulus_sensitivity", adjoint_young_modulus_sensitivity.Evaluate())
            # print("adjoint_young_modulus", adjoint_young_modulus.Evaluate())

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
        intersected_model_part_map = ModelPartUtilities.GetIntersectedMap(self.model_part, merged_model_part_map, True)

        # with open("./primal_parameters.json", 'r') as parameter_file:
        #     primal_parameters = Kratos.Parameters(parameter_file.read())

        # primal_analysis = structural_mechanics_analysis.StructuralMechanicsAnalysis(self.model, primal_parameters)

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

                self.primal_analysis_execution_policy_decorator.Execute()

                objective_value = self.CalculateValue()

                gradient.append((reference_value-objective_value)/perturbation)
                element.Properties[Kratos.YOUNG_MODULUS] = old_stiffness

            return gradient

    def CalculateGradientWithFiniteDifferencing(self, physical_variable_collective_expressions: dict[SupportedSensitivityFieldVariableTypes, KratosOA.ContainerExpression.CollectiveExpressions]) -> None:
        merged_model_part_map = ModelPartUtilities.GetMergedMap(self.model_part, physical_variable_collective_expressions, False)
        intersected_model_part_map = ModelPartUtilities.GetIntersectedMap(self.GetAnalysisModelPart(), merged_model_part_map, True)

        perturbation = 1e-6

        with open("./primal_parameters_import_mdpa.json", 'r') as parameter_file:
            primal_parameters = Kratos.Parameters(parameter_file.read())

        for variable, collective_expression in physical_variable_collective_expressions.items():

            reference_value = self.CalculateValue()

            gradient = []
            num_elements = len(self.model_part.GetElements())

            for i in range(num_elements):

                model_primal = Kratos.Model()
                primal_analysis = structural_mechanics_analysis.StructuralMechanicsAnalysis(model_primal, primal_parameters)
                model_part = model_primal.GetModelPart("Structure")
                # print(model_part)
                primal_analysis.Initialize()

                element: Kratos.Element = model_part.GetElements()[i]

                # model_part.AddProperties(Kratos.Properties.GetTable)

                # for p in model_part.GetProperties():
                #     print(p.GetValue(Kratos.YOUNG_MODULUS))
                #     print(str(p))

                print(element.Properties)
                element.SetValue(Kratos.YOUNG_MODULUS, 0.0)

                stiffness = element.Properties[Kratos.YOUNG_MODULUS]
                element.Properties[Kratos.YOUNG_MODULUS] += perturbation
                # element.SetValue(Kratos.YOUNG_MODULUS, stiffness + perturbation)
                # print("E: ", stiffness, stiffness + perturbation)

                primal_analysis.Run()
                primal_analysis.Finalize()

                objective_value = 0
                for sensor in self.measurement_data["load_cases"][0]["sensors_infos"]:
                    measured_displacement = sensor["measured_value"]
                    measured_direction = Kratos.Array3(sensor["measurement_direction_normal"])

                    mesh_node = model_part.GetNode(sensor["mesh_node_id"])
                    node_displacement = mesh_node.GetSolutionStepValue(Kratos.DISPLACEMENT)
                    in_measurement_direction_projected_vector = (measured_direction[0]*node_displacement[0])+(measured_direction[1]*node_displacement[1])+(measured_direction[2]*node_displacement[2])

                    objective_value += (measured_displacement-in_measurement_direction_projected_vector)**2

                objective_value *= 0.5 * 1/self.measurement_std

                gradient.append((reference_value-objective_value)/perturbation)

                # element.SetValue(Kratos.YOUNG_MODULUS, stiffness)
                element.Properties[Kratos.YOUNG_MODULUS] = stiffness

                # e1: Kratos.Element = element
                # stiffness = e1.GetValue(Kratos.YOUNG_MODULUS)
                # print("E: ", stiffness)
                # e1.SetValue(Kratos.YOUNG_MODULUS, stiffness + perturbation)

                # self.primal_analysis_execution_policy_decorator.Execute()

                # gradient.append(self.CalculateValue()/perturbation)

                # e1.SetValue(Kratos.YOUNG_MODULUS, stiffness)

            print(gradient)

            # self.adjoint_parameters["solver_settings"]["sensitivity_settings"]["element_data_value_sensitivity_variables"].SetStringArray([variable.Name()])

            # self.adjoint_model = Kratos.Model()

            # adjoint_analysis = structural_mechanics_analysis.StructuralMechanicsAnalysis(self.adjoint_model, self.adjoint_parameters)
            # adjoint_analysis.Run()

            # root_model_part_name = self.adjoint_model.GetModelPartNames()[0]
            # Kratos.VariableUtils().CopyModelPartElementalVar(KratosOA.YOUNG_MODULUS_SENSITIVITY, self.adjoint_model[root_model_part_name], self.model[root_model_part_name])

            # # now fill the collective expressions

            # for expression in collective_expression.GetContainerExpressions():
            #     if isinstance(expression, KratosOA.ContainerExpression.ElementPropertiesExpression):
            #         element_data_expression = Kratos.ContainerExpression.ElementNonHistoricalExpression(expression.GetModelPart())
            #         element_data_expression.Read(StructuralMechanicsApplication.YOUNG_MODULUS_SENSITIVITY)
            #         expression.CopyFrom(element_data_expression)
            #     else:
            #         expression.Read(Kratos.KratosGlobals.GetVariable(variable.Name() + "_SENSITIVITY"))

    def __str__(self) -> str:
        return f"Response [type = {self.__class__.__name__}, name = {self.GetName()}, model part name = {self.model_part.FullName()}]"

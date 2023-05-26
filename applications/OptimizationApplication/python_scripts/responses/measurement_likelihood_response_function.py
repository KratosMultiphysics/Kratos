from __future__ import annotations
import json

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
            "measurement_standard_deviation"  : 0.001,
            "perturbation_size"              : 1e-8,
            "evaluated_model_part_names"     : [
                "PLEASE_PROVIDE_A_MODEL_PART_NAME"
            ]
        }""")
        parameters.ValidateAndAssignDefaults(default_settings)

        self.model_part_names = parameters["evaluated_model_part_names"].GetStringArray()
        self.perturbation_size = parameters["perturbation_size"].GetDouble()

        self.measurement_data_file = parameters["measurement_data_file"].GetString()
        self.measurement_std = parameters["measurement_standard_deviation"].GetDouble()

        self.model = model
        self.model_part: Kratos.ModelPart = None
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

        for sensor in self.measurement_data["description_of_sensors"]:

            point = Kratos.Point(sensor["position_of_mesh_node"])
            search_configuration = Kratos.Configuration.Initial
            search_tolerance = 1e-6

            point_locator = Kratos.BruteForcePointLocator(self.model_part)

            found_node_id = point_locator.FindNode(
                point, search_configuration, search_tolerance
            )
            # node = self.model_part.GetNode(found_node_id)
            # print(node)
            sensor["mesh_node_id"] = found_node_id

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
        for sensor in self.measurement_data["description_of_sensors"]:
            measured_displacement = Kratos.Array3(sensor["measured_value"])

            mesh_node = self.model_part.GetNode(sensor["mesh_node_id"])
            node_displacement = mesh_node.GetSolutionStepValue(Kratos.DISPLACEMENT)

            objective_value += (measured_displacement - node_displacement).norm_2()

        objective_value *= 0.5 * 1/self.measurement_std

        # return KratosOA.ResponseUtils.LinearStrainEnergyResponseUtils.CalculateValue(self.model_part)
        return objective_value

    def CalculateGradient(self, physical_variable_collective_expressions: dict[SupportedSensitivityFieldVariableTypes, KratosOA.ContainerExpression.CollectiveExpressions]) -> None:
        # first calculate the gradients
        merged_model_part_map = ModelPartUtilities.GetMergedMap(self.model_part, physical_variable_collective_expressions, False)
        intersected_model_part_map = ModelPartUtilities.GetIntersectedMap(self.GetAnalysisModelPart(), merged_model_part_map, True)

        ###################################################################

        self.response_settings = Kratos.Parameters(
            """{
                "response_type": "adjoint_nodal_displacement",
                "primal_settings": "ProjectParameters.json",
                "adjoint_settings": "auto",
                "primal_data_transfer_with_python": true
                }"""
        )

        self.adjoint_file_name = "/home/fmeister/GitkrakenRepos/Kratos/scripts/dev_scripts/linear_strain_energy_test/linear_shell_test_nodal_disp_adjoint_parameters.json"

        # Reading the ProjectParameters
        with open(self.adjoint_file_name, 'r') as parameter_file:
            self.adjoint_parameters = Kratos.Parameters(parameter_file.read())

        self.model_part_name = self.adjoint_parameters["solver_settings"]["model_part_name"].GetString()

        # create adjoint analysis
        model_adjoint = Kratos.Model()
        self.adjoint_analysis = structural_mechanics_analysis.StructuralMechanicsAnalysis(model_adjoint, self.adjoint_parameters)
        # self.adjoint_analysis.Run()
        # self.adjoint_analysis.time = 0.0
        # self.adjoint_analysis.end_time = 1.0
        # adjoint_response_function = StructuralMechanicsApplication.AdjointNodalDisplacementResponseFunction(self.model_part, self.response_settings)
        # adjoint_response_function.CalculateGradient()

        adjoint_model_part = self.adjoint_analysis.model.GetModelPart(self.model_part_name)
        for element in adjoint_model_part.GetElements():
            print(f"Sensi: {element.GetValue(StructuralMechanicsApplication.YOUNG_MODULUS_SENSITIVITY)}")

        #####################################################################

        # Likelihood Gradient
        # covariance = 0.001
        # S = -(-1/covariance) * b@A
        # Gradient = error @ d_disp/d_Youngs_modulus

        # TODO remove after testing
        KratosOA.ResponseUtils.LinearStrainEnergyResponseUtils.CalculateGradient(list(merged_model_part_map.keys()), list(
            merged_model_part_map.values()), list(intersected_model_part_map.values()), self.perturbation_size)

        # now fill the collective expressions
        for variable, collective_expression in physical_variable_collective_expressions.items():
            collective_expression.Read(Kratos.KratosGlobals.GetVariable(variable.Name() + "_SENSITIVITY"))

    def __str__(self) -> str:
        return f"Response [type = {self.__class__.__name__}, name = {self.GetName()}, model part name = {self.model_part.FullName()}]"

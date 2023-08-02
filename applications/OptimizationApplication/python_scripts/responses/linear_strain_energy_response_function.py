import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.execution_policies.execution_policy_decorator import ExecutionPolicyDecorator
from KratosMultiphysics.OptimizationApplication.responses.response_function import ResponseFunction
from KratosMultiphysics.OptimizationApplication.responses.response_function import SupportedSensitivityFieldVariableTypes
from KratosMultiphysics.OptimizationApplication.utilities.model_part_utilities import ModelPartUtilities

def Factory(model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem) -> ResponseFunction:
    if not parameters.Has("name"):
        raise RuntimeError(f"LinearStrainEnergyResponseFunction instantiation requires a \"name\" in parameters [ parameters = {parameters}].")
    if not parameters.Has("settings"):
        raise RuntimeError(f"LinearStrainEnergyResponseFunction instantiation requires a \"settings\" in parameters [ parameters = {parameters}].")

    return LinearStrainEnergyResponseFunction(parameters["name"].GetString(), model, parameters["settings"], optimization_problem)

class LinearStrainEnergyResponseFunction(ResponseFunction):
    def __init__(self, name: str, model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
        super().__init__(name)

        default_settings = Kratos.Parameters("""{
            "primal_analysis_name"           : "",
            "perturbation_size"              : 1e-8,
            "evaluated_model_part_names"     : [
                "PLEASE_PROVIDE_A_MODEL_PART_NAME"
            ]
        }""")
        parameters.ValidateAndAssignDefaults(default_settings)

        self.model_part_names = parameters["evaluated_model_part_names"].GetStringArray()
        self.perturbation_size = parameters["perturbation_size"].GetDouble()

        self.model = model
        self.model_part: Kratos.ModelPart = None
        self.primal_analysis_execution_policy_decorator: ExecutionPolicyDecorator = optimization_problem.GetExecutionPolicy(parameters["primal_analysis_name"].GetString())

        if len(self.model_part_names) == 0:
            raise RuntimeError("No model parts were provided for LinearStrainEnergyResponseFunction.")

    def GetImplementedPhysicalKratosVariables(self) -> list[SupportedSensitivityFieldVariableTypes]:
        return [KratosOA.SHAPE, Kratos.YOUNG_MODULUS, Kratos.THICKNESS, Kratos.POISSON_RATIO]

    def Initialize(self) -> None:
        model_parts_list = [self.model[model_part_name] for model_part_name in self.model_part_names]
        root_model_part = model_parts_list[0].GetRootModelPart()
        _, self.model_part = ModelPartUtilities.UnionModelParts(root_model_part, model_parts_list, False)

    def Check(self) -> None:
        pass

    def Finalize(self) -> None:
        pass

    def GetEvaluatedModelPart(self) -> Kratos.ModelPart:
        if self.model_part is None:
            raise RuntimeError("Please call LinearStrainEnergyResponseFunction::Initialize first.")
        return self.model_part

    def GetAnalysisModelPart(self) -> Kratos.ModelPart:
        return self.primal_analysis_execution_policy_decorator.GetExecutionPolicy().GetAnalysisModelPart()

    def CalculateValue(self) -> float:
        self.primal_analysis_execution_policy_decorator.GetExecutionPolicy().Execute()
        return KratosOA.ResponseUtils.LinearStrainEnergyResponseUtils.CalculateValue(self.model_part)

    def CalculateGradient(self, physical_variable_collective_expressions: dict[SupportedSensitivityFieldVariableTypes, KratosOA.ContainerExpression.CollectiveExpressions]) -> None:
        # first calculate the gradients
        merged_model_part_map = ModelPartUtilities.GetMergedMap(self.model_part, physical_variable_collective_expressions, False)
        intersected_model_part_map = ModelPartUtilities.GetIntersectedMap(self.GetAnalysisModelPart(), merged_model_part_map, True)
        KratosOA.ResponseUtils.LinearStrainEnergyResponseUtils.CalculateGradient(list(merged_model_part_map.keys()), list(merged_model_part_map.values()), list(intersected_model_part_map.values()), self.perturbation_size)

        # now fill the collective expressions
        for variable, collective_expression in physical_variable_collective_expressions.items():
            collective_expression.Read(Kratos.KratosGlobals.GetVariable(variable.Name() + "_SENSITIVITY"))

    def __str__(self) -> str:
        return f"Response [type = {self.__class__.__name__}, name = {self.GetName()}, model part name = {self.model_part.FullName()}]"
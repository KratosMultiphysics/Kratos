"""
This response is a duplicate of the face_angle response in ShapeOptApplication.
Created on 20.11.2024
"""
from typing import Optional

import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.responses.response_function import ResponseFunction
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import SupportedSensitivityFieldVariableTypes
from KratosMultiphysics.OptimizationApplication.utilities.model_part_utilities import ModelPartOperation
from KratosMultiphysics.OptimizationApplication.utilities.model_part_utilities import ModelPartUtilities

def Factory(model: Kratos.Model, parameters: Kratos.Parameters, _) -> ResponseFunction:
    if not parameters.Has("name"):
        raise RuntimeError(f"FaceAngleResponse instantiation requires a \"name\" in parameters [ parameters = {parameters}].")
    if not parameters.Has("settings"):
        raise RuntimeError(f"FaceAngleResponse instantiation requires a \"settings\" in parameters [ parameters = {parameters}].")

    return FaceAngleResponse(parameters["name"].GetString(), model, parameters["settings"])

def _AddConditionsFromParent(parent, child):
    node_ids = set([node.Id for node in child.Nodes])
    conditions = []
    for condition in parent.Conditions:
        all_nodes_found = True
        for node in condition.GetNodes():
            if node.Id not in node_ids:
                all_nodes_found = False
                break
        if all_nodes_found:
            conditions.append(condition.Id)
    child.AddConditions(conditions)

class FaceAngleResponse(ResponseFunction):
    """
    This response is a duplicate of the face_angle response in ShapeOptApplication:
    -----------------------------
    
    Face angle response function.
    It aggregates the deviation of the face angles of all surface conditions using sqrt(sum(g_i)),
    where g_i are the condition wise violations - feasible conditions do not contribute

    It requires surface conditions in the modelpart, since they are used to compute the face orientation.
    Ideally the design surface model part is used.

    Attributes
    ----------
    model_part : Model part object of the response function
    response_function_utility: Cpp utilities object doing the actual computation of response value and gradient.

    """

    def __init__(self, name: str, model: Kratos.Model, parameters: Kratos.Parameters):
        super().__init__(name)

        default_settings = Kratos.Parameters("""{
            "model_part_name"       : "UNKNOWN_NAME",
            "only"                  : "",
            "domain_size"           : 3,
            "main_direction": [0.0, 0.0, 1.0],
            "min_angle": 0.0,
            "gradient_mode": "finite_differencing",
            "step_size": 1e-6,
            "consider_only_initially_feasible": false
        }""")
            # "response_type"         : "UNKNOWN_TYPE",
            # "model_import_settings" : {
            #     "input_type"        : "use_input_model_part",
            #     "input_filename"    : "UNKNOWN_NAME"
            # },
        parameters.ValidateAndAssignDefaults(default_settings)

        self.response_settings = parameters
        self.model = model

        if parameters["domain_size"] != 3:
            ValueError("Face angle response can only be used on 3D geometries!")

        self.main_direction = parameters["main_direction"].GetVector()
        dir_norm = self.main_direction.norm_2()
        ValueError("Face angle response: 'main_direction' vector norm is 0!")
        self.main_direction /= dir_norm

        self.min_angle = parameters["min_angle"].GetDouble()

        self.perturbation_size = parameters["step_size"].GetDouble()

        self.consider_only_initially_feasible = parameters["consider_only_initially_feasible"].GetBool()
        
        self.value = None
        self._model_part_name = parameters["model_part_name"].GetString()

        only = self.response_settings["only"].GetString()
        if only != "":
            self._model_part_name = only

        self.model_part_operation = ModelPartOperation(self.model, ModelPartOperation.OperationType.UNION, f"response_{self.GetName()}", [self._model_part_name], False)
        self.model_part: Optional[Kratos.ModelPart] = None

    def GetImplementedPhysicalKratosVariables(self) -> 'list[SupportedSensitivityFieldVariableTypes]':
        return [KratosOA.SHAPE]
    
    # function changed name: see below
    # def GetEvaluatedModelPart(self) -> Kratos.ModelPart:
    #     if self.model_part is None:
    #         raise RuntimeError("Please call FaceAngleResponse::Initialize first.")
    #     return self.model_part
    
    def GetInfluencingModelPart(self) -> Kratos.ModelPart:
        if self.model_part is None:
            raise RuntimeError("Please call DiscreteValueResidualResponseFunction::Initialize first.")
        return self.model_part
    
    def GetAnalysisModelPart(self) -> None:
        return None
    
    def Initialize(self) -> None:
        self.model_part = self.model_part_operation.GetModelPart()

        only = self.response_settings["only"].GetString()
        if only != "":
            only_part = self.model.GetModelPart(only)
            if only_part.NumberOfConditions() == 0:
                _AddConditionsFromParent(self.model_part, only_part)
                Kratos.Logger.PrintWarning("FaceAngleResponse", "Automatically added {} conditions to model_part '{}'.".format(only_part.NumberOfConditions(), only_part.Name))
        else:
            only_part = self.model_part

        if only_part.NumberOfConditions() == 0:
            raise RuntimeError("The model_part '{}' does not have any surface conditions!".format(only_part.Name))

        # self.response_function_utility = KratosOA.ResponseUtils.FaceAngleResponseUtils(only_part, self.response_settings)

    def Check(self) -> None:
        pass

    def Finalize(self) -> None:
        pass

    def CalculateValue(self) -> float:
        # import pdb
        # pdb.set_trace()
        self.value = KratosOA.ResponseUtils.FaceAngleResponseUtils.CalculateValue(self.model_part, self.consider_only_initially_feasible, self.main_direction, self.min_angle)
        return self.value
    
    def CalculateGradient(self, physical_variable_collective_expressions: 'dict[SupportedSensitivityFieldVariableTypes, KratosOA.CollectiveExpression]') -> None:
        print("face_angle_response.py :: CalculateGradient")
        # import pdb
        # pdb.set_trace()

        # first merge all the model parts # not necessary for face angle?
        merged_model_part_map = ModelPartUtilities.GetMergedMap(physical_variable_collective_expressions, False)

        # now get the intersected model parts
        intersected_model_part_map = ModelPartUtilities.GetIntersectedMap(self.model_part, merged_model_part_map, False)

        # calculate the gradients
        for physical_variable, merged_model_part in merged_model_part_map.items():
            if physical_variable == KratosOA.SHAPE:
                KratosOA.ResponseUtils.FaceAngleResponseUtils.CalculateGradient(
                    physical_variable,
                    merged_model_part,
                    intersected_model_part_map[physical_variable],
                    physical_variable_collective_expressions[physical_variable].GetContainerExpressions(),
                    self.consider_only_initially_feasible,
                    self.main_direction,
                    self.min_angle,
                    self.perturbation_size)
                
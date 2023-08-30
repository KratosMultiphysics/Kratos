import numpy as np

from typing import Optional

import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.controls.control import Control
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import ContainerExpressionTypes
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import SupportedSensitivityFieldVariableTypes
from KratosMultiphysics.OptimizationApplication.utilities.helper_utilities import IsSameContainerExpression
from KratosMultiphysics.OptimizationApplication.utilities.model_part_utilities import ModelPartOperation


def Factory(model: Kratos.Model, parameters: Kratos.Parameters, _) -> Control:
    if not parameters.Has("name"):
        raise RuntimeError(f"MaterialPropertiesControlSystemIdentification instantiation requires a \"name\" in parameters [ parameters = {parameters}].")
    if not parameters.Has("settings"):
        raise RuntimeError(f"MaterialPropertiesControlSystemIdentification instantiation requires a \"settings\" in parameters [ parameters = {parameters}].")

    return MaterialPropertiesControlSystemIdentification(parameters["name"].GetString(), model, parameters["settings"])


class MaterialPropertiesControlSystemIdentification(Control):
    """Material properties control

    This is a generic material properties control which creates a control
    for the specified control variable. This does not do any filtering.

    TODO: Extend with filtering techniques when they are implemented.

    """

    def __init__(self, name: str, model: Kratos.Model, parameters: Kratos.Parameters):
        super().__init__(name)

        default_settings = Kratos.Parameters("""{
            "model_part_names"     : [""],
            "control_variable_name": "YOUNG_MODULUS",
            "mapping_options":{
                "technique_sequence":[""],
                "smoothing_technique":"element_node_element",
                "num_smoothing_iterations" : 3,
                "set_sensitivity_to_zero_for_nodes": [],
                "exponential_mapping_exponent": 2

            }
        }""")

        # Specific setup
        parameters.ValidateAndAssignDefaults(default_settings)
        self.mapping_options = parameters["mapping_options"]
        self.excluded_nodes = Kratos.Vector()  # self.mapping_options["set_sensitivity_to_zero_for_nodes"].GetVector()

        # Standard setup
        self.model = model

        control_variable_name = parameters["control_variable_name"].GetString()
        control_variable_type = Kratos.KratosGlobals.GetVariableType(control_variable_name)
        if control_variable_type != "Double":
            raise RuntimeError(f"{control_variable_name} with {control_variable_type} type is not supported. Only supports double variables")
        self.controlled_physical_variable = Kratos.KratosGlobals.GetVariable(control_variable_name)

        controlled_model_part_names = parameters["model_part_names"].GetStringArray()
        if len(controlled_model_part_names) == 0:
            raise RuntimeError(f"No model parts were provided for MaterialPropertiesControl. [ control name = \"{self.GetName()}\"]")

        self.model_part_operation = ModelPartOperation(self.model, ModelPartOperation.OperationType.UNION, f"control_{self.GetName()}", controlled_model_part_names, False)
        self.model_part: Optional[Kratos.ModelPart] = None

    def Initialize(self) -> None:
        self.model_part = self.model_part_operation.GetModelPart()

        if not KratosOA.ModelPartUtils.CheckModelPartStatus(self.model_part, "element_specific_properties_created"):
            KratosOA.OptimizationUtils.CreateEntitySpecificPropertiesForContainer(self.model_part, self.model_part.Elements)
            KratosOA.ModelPartUtils.LogModelPartStatus(self.model_part, "element_specific_properties_created")

    def Check(self) -> None:
        pass

    def Finalize(self) -> None:
        pass

    def GetPhysicalKratosVariables(self) -> 'list[SupportedSensitivityFieldVariableTypes]':
        return [self.controlled_physical_variable]

    def GetEmptyField(self) -> ContainerExpressionTypes:
        field = Kratos.Expression.ElementExpression(self.model_part)
        Kratos.Expression.LiteralExpressionIO.SetData(field, 0.0)
        return field

    def GetControlField(self) -> ContainerExpressionTypes:
        field = self.GetEmptyField()
        KratosOA.PropertiesVariableExpressionIO.Read(field, self.controlled_physical_variable)
        return field

    def MapGradient(self, physical_gradient_variable_container_expression_map: 'dict[SupportedSensitivityFieldVariableTypes, ContainerExpressionTypes]') -> ContainerExpressionTypes:
        keys = physical_gradient_variable_container_expression_map.keys()
        kratos_variable = Kratos.KratosGlobals.GetVariable(list(keys)[0].Name() + "_SENSITIVITY")

        if len(keys) != 1:
            raise RuntimeError(f"Provided more than required gradient fields for control \"{self.GetName()}\". Following are the variables:\n\t" + "\n\t".join([k.Name() for k in keys]))
        if self.controlled_physical_variable not in keys:
            raise RuntimeError(
                f"The required gradient for control \"{self.GetName()}\" w.r.t. {self.controlled_physical_variable.Name()} not found. Followings are the variables:\n\t" + "\n\t".join([k.Name() for k in keys]))

        physical_gradient = physical_gradient_variable_container_expression_map[self.controlled_physical_variable]
        if not IsSameContainerExpression(physical_gradient, self.GetEmptyField()):
            raise RuntimeError(
                f"Gradients for the required element container not found for control \"{self.GetName()}\". [ required model part name: {self.model_part.FullName()}, given model part name: {physical_gradient.GetModelPart().FullName()} ]")

        cloned_container_expression = physical_gradient_variable_container_expression_map[self.controlled_physical_variable].Clone()

        mapping_techniques_list = self.mapping_options["technique_sequence"].GetStringArray()

        for technique in mapping_techniques_list:

            if technique == "smoothing":

                if self.mapping_options["smoothing_technique"].GetString() == "element_node_element":
                    self._ElementNodeElementSmoothingForExpressionValues(
                        container_expression=cloned_container_expression)
                else:
                    raise ValueError(self.mapping_options["smoothing_technique"], "Only 'element_node_element' is allowed as a smoothing_technique.")

            elif technique == "exponential_mapping":

                gradient = cloned_container_expression.Evaluate()
                gradient_abs = np.abs(gradient)
                gradient_mapped = (gradient_abs-np.min(gradient_abs))/(np.max(gradient_abs)-np.min(gradient_abs))
                n = self.mapping_options["exponential_mapping_exponent"].GetDouble()
                scaling = ((np.power(np.e, n))-1)
                gradient_mapped = np.exp(gradient_mapped)
                gradient_mapped = np.power(gradient_mapped, n)
                gradient_mapped = gradient_mapped-1
                gradient_mapped /= scaling
                gradient_abs = (gradient_mapped * (np.max(gradient_abs)-np.min(gradient_abs))) + np.min(gradient_abs)
                gradient = gradient_abs * np.sign(gradient)
                Kratos.Expression.CArrayExpressionIO.Read(cloned_container_expression, gradient)

            else:
                raise RuntimeError(f"Mapping technique '{technique}' is not defined. \n Currently only 'smoothing' and 'exponential_mapping' are available.")

        return cloned_container_expression

    def _SetSensitivitiesToZero(self, container_expression: "ContainerExpressionTypes", list_of_ids: Kratos.Vector):
        if list_of_ids.Size() > 0:
            sensitivities = container_expression.Evaluate()

            for ID in list_of_ids:
                sensitivities[int(ID)] = 0

            Kratos.Expression.CArrayExpressionIO.Read(container_expression, sensitivities)

    def _ElementNodeElementSmoothingForExpressionValues(self, container_expression: "ContainerExpressionTypes") -> None:

        neighbors_per_node = Kratos.Expression.NodalExpression(self.model_part)
        KratosOA.ExpressionUtils.ComputeNumberOfNeighbourElements(neighbors_per_node)

        node_sensitivity_expression = Kratos.Expression.NodalExpression(self.model_part)

        for i in range(self.mapping_options["num_smoothing_iterations"].GetInt()):
            KratosOA.ExpressionUtils.MapContainerVariableToNodalVariable(node_sensitivity_expression, container_expression, neighbors_per_node)
            self._SetSensitivitiesToZero(node_sensitivity_expression, self.excluded_nodes)
            KratosOA.ExpressionUtils.MapNodalVariableToContainerVariable(container_expression, node_sensitivity_expression)

    def Update(self, control_field: ContainerExpressionTypes) -> bool:
        if not IsSameContainerExpression(control_field, self.GetEmptyField()):
            raise RuntimeError(
                f"Updates for the required element container not found for control \"{self.GetName()}\". [ required model part name: {self.model_part.FullName()}, given model part name: {control_field.GetModelPart().FullName()} ]")

        # TODO: Implement inverse filtering mechanisms here
        # since no filtering is implemented yet, we are checking the unfiltered updates with the filtered updates. This needs to be changed once the
        # filtering mechanisms are implemented.

        KratosOA.PropertiesVariableExpressionIO.Write(control_field, self.controlled_physical_variable)
        return True

    def __str__(self) -> str:
        return f"Control [type = {self.__class__.__name__}, name = {self.GetName()}, model part name = {self.model_part.FullName()}, control variable = {self.controlled_physical_variable.Name()}"

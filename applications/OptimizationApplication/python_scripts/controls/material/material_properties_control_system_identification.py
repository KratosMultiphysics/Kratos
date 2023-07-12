import numpy as np

import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.controls.control import Control
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import ContainerExpressionTypes
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import SupportedSensitivityFieldVariableTypes
from KratosMultiphysics.OptimizationApplication.utilities.helper_utilities import IsSameContainerExpression
from KratosMultiphysics.OptimizationApplication.utilities.model_part_utilities import ModelPartUtilities


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
                "use_smoothing": false,
                "smoothing_technique":"element_node_element",
                "num_smoothing_iterations" : 3
            }
        }""")
        parameters.ValidateAndAssignDefaults(default_settings)
        self.mapping_options = parameters["mapping_options"]

        self.model = model

        control_variable_name = parameters["control_variable_name"].GetString()
        control_variable_type = Kratos.KratosGlobals.GetVariableType(control_variable_name)
        if control_variable_type != "Double":
            raise RuntimeError(f"{control_variable_name} with {control_variable_type} type is not supported. Only supports double variables")
        self.controlled_physical_variable = Kratos.KratosGlobals.GetVariable(control_variable_name)

        controlled_model_parts = [model[model_part_name] for model_part_name in parameters["model_part_names"].GetStringArray()]
        if len(controlled_model_parts) == 0:
            raise RuntimeError(f"No model parts were provided for MaterialPropertiesControl. [ control name = \"{self.GetName()}\"]")

        self.model_part = ModelPartUtilities.GetOperatingModelPart(ModelPartUtilities.OperationType.UNION, f"control_{self.GetName()}", controlled_model_parts, False)

    def Initialize(self) -> None:
        ModelPartUtilities.ExecuteOperationOnModelPart(self.model_part)

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

        if self.mapping_options["use_smoothing"].GetBool():
            if self.mapping_options["smoothing_technique"].GetString() == "element_node_element":
                self._ElementNodeElementSmoothingForExpressionValues(
                    container_expression=cloned_container_expression,
                    kratos_variable=kratos_variable)
            else:
                ValueError(self.mapping_options["smoothing_technique"], "Only 'element_node_element' is allowed as a smoothing_technique.")

        return cloned_container_expression

    def _ElementNodeElementSmoothingForExpressionValues(self, container_expression: "ContainerExpressionTypes", kratos_variable: "SupportedSensitivityFieldVariableTypes") -> None:

        # num_of_nodes = len(self.model_part.Nodes)
        # elements_per_node_list = np.ndarray(num_of_nodes)

        neighbors_per_node = Kratos.Expression.NodalExpression(self.model_part)
        KratosOA.ExpressionUtils.ComputeNumberOfNeighbourElements(neighbors_per_node)

        nodal_sensitivity_expression = Kratos.Expression.NodalExpression(self.model_part)

        # Kratos.Expression.VariableExpressionIO.Write(container_expression, KratosOA.YOUNG_MODULUS_SENSITIVITY)

        # for i, element in enumerate(container_expression.GetModelPart().Elements):
        #     for node in element.GetNodes():
        #         if not hasattr(node, str(kratos_variable)):
        #             node.AddDof(kratos_variable)

        for i in range(self.mapping_options["num_smoothing_iterations"].GetInt()):
            KratosOA.ExpressionUtils.MapContainerVariableToNodalVariable(nodal_sensitivity_expression, container_expression, neighbors_per_node)
            KratosOA.ExpressionUtils.MapNodalVariableToContainerVariable(container_expression, nodal_sensitivity_expression)
            # elements_per_node_list *= 0.0

            # for i, element in enumerate(container_expression.GetModelPart().Elements):
            #     for node in element.GetNodes():
            #         if elements_per_node_list[node.Id-1] == 0:
            #             node[kratos_variable] = 0.0
            #         node[kratos_variable] += element.GetValue(KratosOA.YOUNG_MODULUS_SENSITIVITY)
            #         elements_per_node_list[node.Id-1] += 1

            # for element in container_expression.GetModelPart().Elements:
            #     new_stiffness = 0.0
            #     for node in element.GetNodes():
            #         new_stiffness += node[kratos_variable] / elements_per_node_list[node.Id-1]
            #     new_stiffness /= len(element.GetNodes())
            #     element.SetValue(KratosOA.YOUNG_MODULUS_SENSITIVITY, new_stiffness)

        # smoothed_sensitivity_expression = Kratos.Expression.ElementExpression(container_expression.GetModelPart())
        # Kratos.Expression.VariableExpressionIO.Read(smoothed_sensitivity_expression, KratosOA.YOUNG_MODULUS_SENSITIVITY)
        # container_expression.SetExpression(smoothed_sensitivity_expression.GetExpression())

    def Update(self, control_field: ContainerExpressionTypes) -> bool:
        if not IsSameContainerExpression(control_field, self.GetEmptyField()):
            raise RuntimeError(
                f"Updates for the required element container not found for control \"{self.GetName()}\". [ required model part name: {self.model_part.FullName()}, given model part name: {control_field.GetModelPart().FullName()} ]")

        # TODO: Implement inverse filtering mechanisms here
        # since no filtering is implemented yet, we are checking the unfiltered updates with the filtered updates. This needs to be changed once the
        # filtering mechanisms are implemented.

        # get the current unfiltered control field
        # unfiltered_control_field = self.GetControlField()

        # if KratosOA.ExpressionUtils.NormL2(unfiltered_control_field - control_field) > 1e-9:
        KratosOA.PropertiesVariableExpressionIO.Write(control_field, self.controlled_physical_variable)
        return True

        # return False

    def __str__(self) -> str:
        return f"Control [type = {self.__class__.__name__}, name = {self.GetName()}, model part name = {self.model_part.FullName()}, control variable = {self.controlled_physical_variable.Name()}"

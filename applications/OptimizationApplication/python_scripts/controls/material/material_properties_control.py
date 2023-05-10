import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.controls.control import Control
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import ContainerExpressionTypes
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import SupportedSensitivityFieldVariableTypes
from KratosMultiphysics.OptimizationApplication.utilities.helper_utilities import IsSameContainerExpression
from KratosMultiphysics.OptimizationApplication.utilities.model_part_utilities import ModelPartUtilities

def Factory(model: Kratos.Model, parameters: Kratos.Parameters, _) -> Control:
    if not parameters.Has("name"):
        raise RuntimeError(f"MaterialPropertiesControl instantiation requires a \"name\" in parameters [ parameters = {parameters}].")
    if not parameters.Has("settings"):
        raise RuntimeError(f"MaterialPropertiesControl instantiation requires a \"settings\" in parameters [ parameters = {parameters}].")

    return MaterialPropertiesControl(parameters["name"].GetString(), model, parameters["settings"])

class MaterialPropertiesControl(Control):
    """Material properties control

    This is a generic material properties control which creates a control
    for the specified control variable. This does not do any filtering.

    TODO: Extend with filtering techniques when they are implemented.

    """
    def __init__(self, name: str, model: Kratos.Model, parameters: Kratos.Parameters):
        super().__init__(name)

        default_settings = Kratos.Parameters("""{
            "model_part_names"     : [""],
            "control_variable_name": ""
        }""")
        parameters.ValidateAndAssignDefaults(default_settings)

        self.model_part_names = parameters["model_part_names"].GetStringArray()
        self.model = model

        control_variable_name = parameters["control_variable_name"].GetString()
        control_variable_type = Kratos.KratosGlobals.GetVariableType(control_variable_name)

        if control_variable_type != "Double":
            raise RuntimeError(f"{control_variable_name} with {control_variable_type} type is not supported. Only supports double variables")
        if len(self.model_part_names) == 0:
            raise RuntimeError(f"\"{self.GetName()}\" control has not specified any model part names.")

        self.controlled_physical_variable = Kratos.KratosGlobals.GetVariable(control_variable_name)

    def Initialize(self) -> None:
        model_parts_list = [self.model[model_part_name] for model_part_name in self.model_part_names]
        root_model_part = model_parts_list[0].GetRootModelPart()
        is_new_model_part, self.model_part = ModelPartUtilities.UnionModelParts(root_model_part, model_parts_list, False)

        if is_new_model_part:
            # now create entity specific properties for the merged model part which is used for the control.
            KratosOA.OptimizationUtils.CreateEntitySpecificPropertiesForContainer(self.model_part, self.model_part.Elements)

    def Check(self) -> None:
        pass

    def Finalize(self) -> None:
        pass

    def GetPhysicalKratosVariables(self) -> list[SupportedSensitivityFieldVariableTypes]:
        return [self.controlled_physical_variable]

    def GetEmptyField(self) -> ContainerExpressionTypes:
        field = KratosOA.ContainerExpression.ElementPropertiesExpression(self.model_part)
        field.SetData(0.0)
        return field

    def GetControlField(self) -> ContainerExpressionTypes:
        field = self.GetEmptyField()
        field.Read(self.controlled_physical_variable)
        return field

    def MapGradient(self, physical_gradient_variable_container_expression_map: dict[SupportedSensitivityFieldVariableTypes, ContainerExpressionTypes]) -> ContainerExpressionTypes:
        keys = physical_gradient_variable_container_expression_map.keys()
        if len(keys) != 1:
            raise RuntimeError(f"Provided more than required gradient fields for control \"{self.GetName()}\". Following are the variables:\n\t" + "\n\t".join([k.Name() for k in keys]))
        if self.controlled_physical_variable not in keys:
            raise RuntimeError(f"The required gradient for control \"{self.GetName()}\" w.r.t. {self.controlled_physical_variable.Name()} not found. Followings are the variables:\n\t" + "\n\t".join([k.Name() for k in keys]))

        physical_gradient = physical_gradient_variable_container_expression_map[self.controlled_physical_variable]
        if not IsSameContainerExpression(physical_gradient, self.GetEmptyField()):
            raise RuntimeError(f"Gradients for the required element container not found for control \"{self.GetName()}\". [ required model part name: {self.model_part.FullName()}, given model part name: {physical_gradient.GetModelPart().FullName()} ]")

        # TODO: Implement filtering mechanisms here
        return physical_gradient_variable_container_expression_map[self.controlled_physical_variable].Clone()

    def Update(self, control_field: ContainerExpressionTypes) -> bool:
        if not IsSameContainerExpression(control_field, self.GetEmptyField()):
            raise RuntimeError(f"Updates for the required element container not found for control \"{self.GetName()}\". [ required model part name: {self.model_part.FullName()}, given model part name: {control_field.GetModelPart().FullName()} ]")

        # TODO: Implement inverse filtering mechanisms here
        # since no filtering is implemented yet, we are checking the unfiltered updates with the filtered updates. This needs to be changed once the
        # filtering mechanisms are implemented.

        # get the current unfiltered control field
        unfiltered_control_field = self.GetEmptyField()
        unfiltered_control_field.Read(self.controlled_physical_variable)

        if KratosOA.ContainerExpressionUtils.NormL2(unfiltered_control_field - control_field) > 1e-9:
            control_field.Evaluate(self.controlled_physical_variable)
            return True

        return False

    def __str__(self) -> str:
        return f"Control [type = {self.__class__.__name__}, name = {self.GetName()}, model part name = {self.model_part.FullName()}, control variable = {self.controlled_physical_variable.Name()}"

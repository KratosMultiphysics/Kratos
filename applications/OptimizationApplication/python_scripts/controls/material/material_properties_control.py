import numpy
from typing import Optional

import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.controls.control import Control
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import SupportedSensitivityFieldVariableTypes
from KratosMultiphysics.OptimizationApplication.utilities.model_part_utilities import ModelPartOperation

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
            "model_part_names"                  : [""],
            "control_variable_name"             : "",
            "consider_recursive_property_update": false
        }""")
        parameters.ValidateAndAssignDefaults(default_settings)

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

        self.consider_recursive_property_update = parameters["consider_recursive_property_update"].GetBool()

    def Initialize(self) -> None:
        self.model_part = self.model_part_operation.GetModelPart()

        if not KratosOA.OptAppModelPartUtils.CheckModelPartStatus(self.model_part, "element_specific_properties_created"):
            KratosOA.OptimizationUtils.CreateEntitySpecificPropertiesForContainer(self.model_part, self.model_part.Elements, self.consider_recursive_property_update)
            KratosOA.OptAppModelPartUtils.LogModelPartStatus(self.model_part, "element_specific_properties_created")

    def Check(self) -> None:
        pass

    def Finalize(self) -> None:
        pass

    def GetPhysicalKratosVariables(self) -> 'list[SupportedSensitivityFieldVariableTypes]':
        return [self.controlled_physical_variable]

    def GetEmptyField(self) -> Kratos.TensorAdaptors.DoubleTensorAdaptor:
        field = KratosOA.TensorAdaptors.PropertiesVariableTensorAdaptor(self.model_part.Elements, self.controlled_physical_variable)
        field.data[:] = 0.0
        return Kratos.TensorAdaptors.DoubleTensorAdaptor(field, copy=False)

    def GetControlField(self) -> Kratos.TensorAdaptors.DoubleTensorAdaptor:
        field = KratosOA.TensorAdaptors.PropertiesVariableTensorAdaptor(self.model_part.Elements, self.controlled_physical_variable)
        field.CollectData()
        return Kratos.TensorAdaptors.DoubleTensorAdaptor(field, copy=False)

    def MapGradient(self, physical_gradient_variable_container_expression_map: 'dict[SupportedSensitivityFieldVariableTypes, Kratos.TensorAdaptors.DoubleTensorAdaptor]') -> Kratos.TensorAdaptors.DoubleTensorAdaptor:
        keys = physical_gradient_variable_container_expression_map.keys()
        if len(keys) != 1:
            raise RuntimeError(f"Provided more than required gradient fields for control \"{self.GetName()}\". Following are the variables:\n\t" + "\n\t".join([k.Name() for k in keys]))
        if self.controlled_physical_variable not in keys:
            raise RuntimeError(f"The required gradient for control \"{self.GetName()}\" w.r.t. {self.controlled_physical_variable.Name()} not found. Followings are the variables:\n\t" + "\n\t".join([k.Name() for k in keys]))

        physical_gradient = physical_gradient_variable_container_expression_map[self.controlled_physical_variable]
        if physical_gradient.GetContainer() != self.GetEmptyField().GetContainer():
            raise RuntimeError(f"Gradients for the required element container not found for control \"{self.GetName()}\". [ required model part name: {self.model_part.FullName()} ]")

        # TODO: Implement filtering mechanisms here
        return Kratos.TensorAdaptors.DoubleTensorAdaptor(physical_gradient_variable_container_expression_map[self.controlled_physical_variable])

    def Update(self, control_field: Kratos.TensorAdaptors.DoubleTensorAdaptor) -> bool:
        if control_field.GetContainer() != self.GetEmptyField().GetContainer():
            raise RuntimeError(f"Updates for the required element container not found for control \"{self.GetName()}\". [ required model part name: {self.model_part.FullName()} ]")

        # TODO: Implement inverse filtering mechanisms here
        # since no filtering is implemented yet, we are checking the unfiltered updates with the filtered updates. This needs to be changed once the
        # filtering mechanisms are implemented.

        # get the current unfiltered control field
        unfiltered_control_field = self.GetControlField()

        if numpy.linalg.norm(unfiltered_control_field.data - control_field.data) > 1e-9:
            KratosOA.TensorAdaptors.PropertiesVariableTensorAdaptor(control_field, self.controlled_physical_variable, copy=False).StoreData()

            if self.consider_recursive_property_update:
                KratosOA.OptimizationUtils.UpdatePropertiesVariableWithRootValueRecursively(control_field.GetContainer(), self.controlled_physical_variable)
            return True

        return False

    def __str__(self) -> str:
        return f"Control [type = {self.__class__.__name__}, name = {self.GetName()}, model part name = {self.model_part.FullName()}, control variable = {self.controlled_physical_variable.Name()}"

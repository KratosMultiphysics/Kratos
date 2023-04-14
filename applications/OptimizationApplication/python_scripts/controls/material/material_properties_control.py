import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.controls.control import Control
from KratosMultiphysics.OptimizationApplication.utilities.optimization_info import OptimizationInfo
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import SupportedSensitivityFieldVariableTypes

class MaterialPropertiesControl(Control):
    """Material properties control

    This is a generic material properties control which creates a control
    for the specified control variable. This does not do any filtering.

    TODO: Extend with filtering techniques when they are implemented.

    """
    def __init__(self, model: Kratos.Model, parameters: Kratos.Parameters, optimization_info: OptimizationInfo):
        super().__init__()

        default_settings = Kratos.Parameters("""{
            "model_part_names"            : [""],
            "control_variable_name"       : ""
        }""")
        parameters.ValidateAndAssignDefaults(default_settings)

        self.optimization_info = optimization_info
        self.model_parts = [model[model_part_name] for model_part_name in parameters["model_part_names"].GetStringArray()]

        control_variable_name = parameters["control_variable_name"].GetString()
        control_variable_type = Kratos.KratosGlobals.GetVariableType(control_variable_name)
        if control_variable_type != "Double":
            raise RuntimeError(f"{control_variable_name} with {control_variable_type} type is not supported. Only supports double variables")

        sensitivity_variable_name = control_variable_name + "_SENSITIVITY"
        if not Kratos.KratosGlobals.HasVariable(sensitivity_variable_name):
            raise RuntimeError(f"Sensitivity variable {sensitivity_variable_name} for control variable {control_variable_name} not found.")
        elif Kratos.KratosGlobals.GetVariableType(sensitivity_variable_name) != "Double":
            raise RuntimeError(f"Sensitivity variable {sensitivity_variable_name} type is not double.")

        self.control_variable = Kratos.KratosGlobals.GetVariable(control_variable_name)
        self.sensitivity_variable: SupportedSensitivityFieldVariableTypes = Kratos.KratosGlobals.GetVariable(sensitivity_variable_name)

    def Initialize(self) -> None:
        # creates element specific properties
        problem_data_container = self.optimization_info.GetProblemDataContainer()
        if not problem_data_container.HasValue("model_parts_with_element_specific_properties"):
            problem_data_container["model_parts_with_element_specific_properties"] = []

        for model_part in self.model_parts:
            # check if element specific properties are created for the model part by a different process.
            if not f"{model_part.FullName()}.Elements" in problem_data_container["model_parts_with_element_specific_properties"]:
                KratosOA.OptimizationUtils.CreateEntitySpecificPropertiesForContainer(model_part, model_part.Elements)
                problem_data_container["model_parts_with_element_specific_properties"].append(f"{model_part.FullName()}.Elements")

    def GetRequiredSensitivityFieldVariables(self) -> 'list[SupportedSensitivityFieldVariableTypes]':
        return [self.sensitivity_variable]

    def GetCollectiveExpression(self) -> KratosOA.ContainerExpression.CollectiveExpressions:
        control_collective_expression = KratosOA.ContainerExpression.CollectiveExpressions()
        for model_part in self.model_parts:
            control_collective_expression.Add(KratosOA.ContainerExpression.ElementPropertiesExpression(model_part))

        return control_collective_expression

    def MapFirstDerivative(self, sensitivity_variable_first_derivative_map: 'dict[SupportedSensitivityFieldVariableTypes, KratosOA.ContainerExpression.CollectiveExpressions]', mapped_first_derivative: KratosOA.ContainerExpression.CollectiveExpressions):
        # first check the input is correct
        self._CheckSensitivityCollectiveExpressions(sensitivity_variable_first_derivative_map)

        # clear the container
        mapped_first_derivative.Clear()

        for provided_container_expression in sensitivity_variable_first_derivative_map[self.sensitivity_variable].GetContainerExpressions():
            # we are cloning in here because, it is not sure what the algorithm will
            # do with the sensitivity container expressions. This is just a move operation
            # of the pointers the expression is holding, hence not at all expensive operation.
            # it does not also clone the memory space.
            mapped_first_derivative.Add(provided_container_expression.Clone())

    def Update(self, update: KratosOA.ContainerExpression.CollectiveExpressions) -> None:
        # first check the input is correct
        self._CheckUpdateCollectiveExpressions(update)

        for provided_container_expression in update.GetContainerExpressions():
            # first read the existing values from the element properties expression
            current_values_container = KratosOA.ContainerExpression.ElementPropertiesExpression(provided_container_expression.GetModelPart())
            current_values_container.Read(self.control_variable)

            # now update the properties
            new_values_container: KratosOA.ContainerExpression.ElementPropertiesExpression = current_values_container + provided_container_expression
            new_values_container.Evaluate(self.control_variable)
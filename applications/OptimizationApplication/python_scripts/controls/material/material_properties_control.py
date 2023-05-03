import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.controls.control import Control
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import ContainerExpressionTypes
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import SupportedSensitivityFieldVariableTypes
from KratosMultiphysics.OptimizationApplication.utilities.helper_utilities import IsSameContainerExpression

def Factory(model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
    return MaterialPropertiesControl(model, parameters, optimization_problem)

class MaterialPropertiesControl(Control):
    """Material properties control

    This is a generic material properties control which creates a control
    for the specified control variable. This does not do any filtering.

    TODO: Extend with filtering techniques when they are implemented.

    """
    def __init__(self, model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
        super().__init__()

        default_settings = Kratos.Parameters("""{
            "combined_output_model_part_name": "<CONTROL_NAME>_combined_no_neighbours_control",
            "model_part_names"               : [""],
            "control_variable_name"          : ""
        }""")
        parameters.ValidateAndAssignDefaults(default_settings)

        self.output_model_part_name = parameters["combined_output_model_part_name"].GetString()
        self.model_part_names = parameters["model_part_names"].GetStringArray()
        self.optimization_problem = optimization_problem
        self.model = model

        control_variable_name = parameters["control_variable_name"].GetString()
        control_variable_type = Kratos.KratosGlobals.GetVariableType(control_variable_name)
        if control_variable_type != "Double":
            raise RuntimeError(f"{control_variable_name} with {control_variable_type} type is not supported. Only supports double variables")

        self.control_variable = Kratos.KratosGlobals.GetVariable(control_variable_name)

    def Initialize(self) -> None:
        # get the model part name
        self.control_name = self.optimization_problem.GetComponentName(self)
        output_model_part_name = self.output_model_part_name.replace("<CONTROL_NAME>", self.control_name)

        # get root model part
        model_parts_list = [self.model[model_part_name] for model_part_name in self.model_part_names]
        root_model_part = model_parts_list[0].GetRootModelPart()

        # create the combined model part
        if not root_model_part.HasSubModelPart(output_model_part_name):
            self.model_part: Kratos.ModelPart = Kratos.ModelPartOperationUtilities.Merge(output_model_part_name, root_model_part, model_parts_list, False)
            # now create entity specific properties for the merged model part which is used for the control.
            KratosOA.OptimizationUtils.CreateEntitySpecificPropertiesForContainer(self.model_part, self.model_part.Elements)
        else:
            self.model_part = root_model_part.GetSubModelPart(output_model_part_name)

    def Check(self) -> None:
        pass

    def Finalize(self) -> None:
        pass

    def GetPhysicalKratosVariables(self) -> list[SupportedSensitivityFieldVariableTypes]:
        return [self.control_variable]

    def GetEmptyControlField(self) -> ContainerExpressionTypes:
        return KratosOA.ContainerExpression.ElementPropertiesExpression(self.model_part)

    def MapGradient(self, physical_gradient_variable_container_expression_map: dict[SupportedSensitivityFieldVariableTypes, ContainerExpressionTypes]) -> ContainerExpressionTypes:
        keys = physical_gradient_variable_container_expression_map.keys()
        if len(keys) != 1:
            raise RuntimeError(f"Provided more than required gradient fields for control \"{self.control_name}\". Following are the variables:\n\t" + "\n\t".join([k.Name() for k in keys]))
        if self.control_variable not in keys:
            raise RuntimeError(f"The required gradient for control \"{self.control_name}\" w.r.t. {self.control_variable.Name()} not found. Followings are the variables:\n\t" + "\n\t".join([k.Name() for k in keys]))

        physical_gradient = physical_gradient_variable_container_expression_map[self.control_variable]
        if not IsSameContainerExpression(physical_gradient, self.GetEmptyControlField()):
            raise RuntimeError(f"Gradients for the required element container not found for control \"{self.control_name}\". [ required model part name: {self.model_part.FullName()}, given model part name: {physical_gradient.GetModelPart().FullName()} ]")

        # TODO: Implement filtering mechanisms here
        return physical_gradient_variable_container_expression_map[self.control_variable].Clone()

    def Update(self, control_field: ContainerExpressionTypes) -> bool:
        if not IsSameContainerExpression(control_field, self.GetEmptyControlField()):
            raise RuntimeError(f"Updates for the required element container not found for control \"{self.control_name}\". [ required model part name: {self.model_part.FullName()}, given model part name: {control_field.GetModelPart().FullName()} ]")

        # TODO: Implement inverse filtering mechanisms here
        # since no filtering is implemented yet, we are checking the unfiltered updates with the filtered updates. This needs to be changed once the
        # filtering mechanisms are implemented.

        # get the current unfiltered control field
        unfiltered_control_field = self.GetEmptyControlField()
        unfiltered_control_field.Read(self.control_variable)

        if KratosOA.ContainerExpressionUtils.NormL2(unfiltered_control_field - control_field) > 1e-9:
            control_field.Evaluate(self.control_variable)
            return True

        return False
import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.controls.control import Control
from KratosMultiphysics.OptimizationApplication.optimization_info import OptimizationInfo
from KratosMultiphysics.OptimizationApplication.utilities.helper_utils import ContainerEnum
from KratosMultiphysics.OptimizationApplication.execution_policies.execution_policy_wrapper import ExecutionPolicyWrapper

class ShapeControl(Control):
    def __init__(self, model: Kratos.Model, parameters: Kratos.Parameters, optimization_info: OptimizationInfo):
        super().__init__(model, parameters, optimization_info)

        default_settings = Kratos.Parameters("""{
            "model_part_name"             : "",
            "control_update_variable_name": "VECTOR_CONTROL_UPDATE",
            "mesh_moving_analysis_name"   : ""
        }""")
        parameters.ValidateAndAssignDefaults(default_settings)

        self.model_part = self.model[parameters["model_part_name"].GetString()]

        control_update_variable_name = parameters["control_update_variable_name"].GetString()
        control_update_variable_type = Kratos.KratosGlobals.GetVariableType(control_update_variable_name)
        if control_update_variable_type != "Array":
            raise RuntimeError(f"{control_update_variable_name} with {control_update_variable_type} type is not supported. Only supports array variables")

        self.control_update_variable = Kratos.KratosGlobals.GetVariable(control_update_variable_name)
        self.mesh_moving_execution_policy_wrapper: ExecutionPolicyWrapper = self.optimization_info.GetOptimizationRoutine("ExecutionPolicyWrapper", self.parameters["mesh_moving_analysis_name"].GetString())

    def UpdateControls(self, control_values: Kratos.Vector):
        KratosOA.OptimizationUtils.AssignVectorToHistoricalContainer(self.model_part, self.model_part.ProcessInfo[Kratos.DOMAIN_SIZE], Kratos.MESH_DISPLACEMENT, control_values)
        self.mesh_moving_execution_policy_wrapper.Execute()

    def GetNewControlValuesVector(self) -> Kratos.Vector:
        return self.GetControlUpdatesVector()

    def GetModelPart(self) -> Kratos.ModelPart:
        return self.model_part

    def GetContainerType(self) -> ContainerEnum:
        return ContainerEnum.NODES

    def GetControlSensitivityVariable(self) -> any:
        return Kratos.SHAPE_SENSITIVITY

    def GetControlUpdateVariable(self) -> any:
        return self.control_update_variable


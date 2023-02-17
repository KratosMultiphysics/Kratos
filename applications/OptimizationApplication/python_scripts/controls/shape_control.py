import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.controls.control import Control
from KratosMultiphysics.OptimizationApplication.optimization_info import OptimizationInfo
from KratosMultiphysics.OptimizationApplication.utilities.execution_policy_decorator import ExecutionPolicyDecorator
from KratosMultiphysics.OptimizationApplication.utilities.response_function_implementor import ResponseFunctionImplementor

class ShapeControl(Control):
    def __init__(self, model: Kratos.Model, parameters: Kratos.Parameters, optimization_info: OptimizationInfo):
        super().__init__()

        default_settings = Kratos.Parameters("""{
            "model_part_names"            : [""],
            "mesh_moving_analysis_name"   : ""
        }""")
        parameters.ValidateAndAssignDefaults(default_settings)

        self.model_parts = [model[model_part_name] for model_part_name in parameters["model_part_names"].GetStringArray()]
        self.mesh_moving_execution_policy_wrapper: ExecutionPolicyDecorator = optimization_info.GetOptimizationProcess(ExecutionPolicyDecorator, parameters["mesh_moving_analysis_name"].GetString())

    def CalculateSensitivity(self, response_function: ResponseFunctionImplementor, output_sensitivities: KratosOA.CollectiveVariableDataHolder):
        # clear the container
        output_sensitivities.ClearVariableDataHolders()

        for model_part in self.model_parts:
            d_j_d_x = KratosOA.NodalContainerVariableDataHolder(model_part)
            response_function.GetStandardizedSensitivity(Kratos.SHAPE_SENSITIVITY, d_j_d_x)

            output_sensitivities.AddVariableDataHolder(d_j_d_x)

    def UpdateControl(self, update: KratosOA.CollectiveVariableDataHolder):
        for current_update in update.GetVariableDataHolders():
            model_part: Kratos.ModelPart = current_update.GetModelPart()

            if model_part not in self.model_parts:
                raise RuntimeError(f"Trying to update {model_part.FullName()} which is not controlled by {self.GetName()}.")

            historical_container = KratosOA.HistoricalContainerVariableDataHolder(current_update)
            historical_container.AssignDataToContainerVariable(Kratos.MESH_DISPLACEMENT)

        # run the mesh motion solver
        self.mesh_moving_execution_policy_wrapper.Execute()

        # update initial mesh to the current mesh
        for model_part in self.model_parts:
            Kratos.VariableUtils().UpdateCurrentToInitialConfiguration(model_part.Nodes)

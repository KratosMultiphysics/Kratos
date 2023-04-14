import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.controls.control import Control
from KratosMultiphysics.OptimizationApplication.utilities.optimization_info import OptimizationInfo
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import SupportedSensitivityFieldVariableTypes
from KratosMultiphysics.OptimizationApplication.execution_policies.stepping_analysis_execution_policy import SteppingAnalysisExecutionPolicy
from KratosMultiphysics.OptimizationApplication.execution_policies.independent_analysis_execution_policy import IndependentAnalysisExecutionPolicy

class ShapeControlNew(Control):
    """Shape control

    This class is planned to replace existing "shape_control.py" when the filtering
    techniques are introduced to the new OptimizationApplication framework.

    TODO: Update the original shape_control.py when filtering techniques are done.
          Add support for filtering techniques.

    """
    def __init__(self, model: Kratos.Model, parameters: Kratos.Parameters, _: OptimizationInfo):
        super().__init__()

        default_settings = Kratos.Parameters("""{
            "model_part_names"                 : [""],
            "mesh_motion_execution_policy_type": "SteppingAnalysisExecutionPolicy",
            "mesh_motion_module"               : "KratosMultiphysics.MeshMovingApplication",
            "mesh_motion_analysis_type"        : "MeshMovingAnalysis",
            "mesh_motion_settings"             : {}
        }""")
        parameters.ValidateAndAssignDefaults(default_settings)

        # create stepping analysis parameters
        execution_policy_parameters = Kratos.Parameters()
        execution_policy_parameters.AddString("analysis_module", parameters["mesh_motion_module"].GetString())
        execution_policy_parameters.AddString("analysis_type", parameters["mesh_motion_analysis_type"].GetString())
        execution_policy_parameters.AddValue("analysis_settings", parameters["mesh_motion_settings"])

        self.model_parts = [model[model_part_name] for model_part_name in parameters["model_part_names"].GetStringArray()]

        mesh_motion_execution_policy_type = parameters["mesh_motion_execution_policy_type"].GetString()
        if mesh_motion_execution_policy_type == "SteppingAnalysisExecutionPolicy":
            execution_policy_parameters.AddStringArray("model_part_names", parameters["model_part_names"].GetStringArray())
            self.mesh_moving_execution_policy = SteppingAnalysisExecutionPolicy(model, execution_policy_parameters, None)
        elif mesh_motion_execution_policy_type == "IndependentAnalysisExecutionPolicy":
            self.mesh_moving_execution_policy = IndependentAnalysisExecutionPolicy(model, execution_policy_parameters, None)
        else:
            raise RuntimeError(f"Only supports \"SteppingAnalysisExecutionPolicy\" or \"IndependentAnalysisExecutionPolicy\" type for \"mesh_motion_execution_policy_type\" [ Provided \"mesh_motion_execution_policy_type\" = \"{mesh_motion_execution_policy_type}\" ].")

    def GetRequiredSensitivityFieldVariables(self) -> 'list[SupportedSensitivityFieldVariableTypes]':
        return [Kratos.SHAPE_SENSITIVITY]

    def GetCollectiveExpression(self) -> KratosOA.ContainerExpression.CollectiveExpressions:
        control_collective_expression = KratosOA.ContainerExpression.CollectiveExpressions()
        for model_part in self.model_parts:
            control_collective_expression.Add(Kratos.ContainerExpression.NodalNonHistoricalExpression(model_part))

        return control_collective_expression

    def MapFirstDerivative(self, sensitivity_variable_first_derivative_map: 'dict[SupportedSensitivityFieldVariableTypes, KratosOA.ContainerExpression.CollectiveExpressions]', mapped_first_derivative: KratosOA.ContainerExpression.CollectiveExpressions):
        # first check the input is correct
        self._CheckSensitivityCollectiveExpressions(sensitivity_variable_first_derivative_map)

        # clear the container
        mapped_first_derivative.Clear()

        for provided_container_expression in sensitivity_variable_first_derivative_map[Kratos.SHAPE_SENSITIVITY].GetContainerExpressions():
            # we are cloning in here because, it is not sure what the algorithm will
            # do with the sensitivity container expressions. This is just a move operation
            # of the pointers the expression is holding, hence not at all expensive operation.
            # it does not also clone the memory space.
            mapped_first_derivative.Add(provided_container_expression.Clone())

    def Update(self, update: KratosOA.ContainerExpression.CollectiveExpressions) -> None:
        # first check the input is correct
        self._CheckUpdateCollectiveExpressions(update)

        for provided_container_expression in update.GetContainerExpressions():
            model_part: Kratos.ModelPart = provided_container_expression.GetModelPart()

            # transfer and assign non-historical data to historical variable
            historical_container = Kratos.ContainerExpression.HistoricalExpression(provided_container_expression)
            historical_container.Evaluate(Kratos.MESH_DISPLACEMENT)

        # run the mesh motion solver
        self.mesh_moving_execution_policy.Execute()

        # update initial mesh to the current mesh
        for model_part in self.model_parts:
            Kratos.VariableUtils().UpdateCurrentToInitialConfiguration(model_part.Nodes)

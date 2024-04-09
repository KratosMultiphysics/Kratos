import typing
import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
import KratosMultiphysics.DigitalTwinApplication as KratosDT
from KratosMultiphysics.OptimizationApplication.controls.control import Control
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import ContainerExpressionTypes
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import SupportedSensitivityFieldVariableTypes
from KratosMultiphysics.OptimizationApplication.utilities.model_part_utilities import ModelPartOperation
from KratosMultiphysics.OptimizationApplication.utilities.helper_utilities import IsSameContainerExpression
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView
from KratosMultiphysics.OptimizationApplication.filtering.filter import Factory as FilterFactory

def Factory(model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem) -> Control:
    if not parameters.Has("name"):
        raise RuntimeError(f"SensorStatusControl instantiation requires a \"name\" in parameters [ parameters = {parameters}].")
    if not parameters.Has("settings"):
        raise RuntimeError(f"SensorStatusControl instantiation requires a \"settings\" in parameters [ parameters = {parameters}].")
    return SensorStatusControl(parameters["name"].GetString(), model, parameters["settings"], optimization_problem)

class SensorStatusControl(Control):
    def __init__(self, name: str, model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem) -> None:
        super().__init__(name)

        self.parameters = parameters
        self.model = model
        self.optimization_problem = optimization_problem

        default_settings = Kratos.Parameters("""{
            "controlled_model_part_names": [""],
            "initial_sensor_status"      : 0.5,
            "filter_settings"            : {},
            "beta_settings": {
                "initial_value": 5,
                "max_value"    : 30,
                "adaptive"     : true,
                "increase_fac" : 1.05,
                "update_period": 50
            }
        }""")
        parameters.ValidateAndAssignDefaults(default_settings)

        self.controlled_physical_variables = [KratosDT.SENSOR_STATUS]
        self.initial_sensor_status = parameters["initial_sensor_status"].GetDouble()

        # beta settings
        beta_settings = parameters["beta_settings"]
        beta_settings.ValidateAndAssignDefaults(default_settings["beta_settings"])
        self.beta = beta_settings["initial_value"].GetDouble()
        self.beta_max = beta_settings["max_value"].GetDouble()
        self.beta_adaptive = beta_settings["adaptive"].GetBool()
        self.beta_increase_frac = beta_settings["increase_fac"].GetDouble()
        self.beta_update_period = beta_settings["update_period"].GetInt()
        self.beta_computed_step = 1

        controlled_model_names_parts = parameters["controlled_model_part_names"].GetStringArray()
        if len(controlled_model_names_parts) == 0:
            raise RuntimeError(f"No model parts are provided for SensorStatusControl. [ control name = \"{self.GetName()}\"]")
        self.model_part_operation = ModelPartOperation(self.model, ModelPartOperation.OperationType.UNION, f"control_{self.GetName()}", controlled_model_names_parts, False)
        self.model_part: 'typing.Optional[Kratos.ModelPart]' = None

        # filtering settings
        self.filter = FilterFactory(self.model, self.model_part_operation.GetModelPartFullName(), KratosDT.SENSOR_STATUS, Kratos.Globals.DataLocation.NodeNonHistorical, parameters["filter_settings"])

    def Initialize(self) -> None:
        self.un_buffered_data = ComponentDataView(self, self.optimization_problem).GetUnBufferedData()

        # init model parts
        self.model_part = self.model_part_operation.GetModelPart()

        # calculate phi from existing density
        sensor_status = Kratos.Expression.NodalExpression(self.model_part)
        Kratos.Expression.LiteralExpressionIO.SetData(sensor_status, self.initial_sensor_status)
        self.phi = KratosOA.ControlUtils.SigmoidalProjectionUtils.ProjectBackward(sensor_status, [0, 1], [0, 1], self.beta, 1)
        self.un_buffered_data.SetValue("phi", self.phi, overwrite=True)
        self.__ApplySensorStatus()

        self.filter.SetComponentDataView(ComponentDataView(self, self.optimization_problem))
        self.filter.Initialize()

    def Check(self) -> None:
        self.filter.Check()

    def Finalize(self) -> None:
        self.filter.Finalize()

    def GetPhysicalKratosVariables(self) -> 'list[SupportedSensitivityFieldVariableTypes]':
        return self.controlled_physical_variables

    def GetEmptyField(self) -> ContainerExpressionTypes:
        field = Kratos.Expression.NodalExpression(self.model_part)
        Kratos.Expression.LiteralExpressionIO.SetData(field, 0.0)
        return field

    def GetControlField(self) -> ContainerExpressionTypes:
        return self.phi.Clone()

    def MapGradient(self, physical_gradient_variable_container_expression_map: 'dict[SupportedSensitivityFieldVariableTypes, ContainerExpressionTypes]') -> ContainerExpressionTypes:
        keys = physical_gradient_variable_container_expression_map.keys()
        if len(keys) != 1:
            raise RuntimeError(f"Not provided required gradient fields for control \"{self.GetName()}\". Following are the variables:\n\t" + "\n\t".join([k.Name() for k in keys]))
        if  KratosDT.SENSOR_STATUS not in keys:
            raise RuntimeError(f"The required gradient for control \"{self.GetName()}\" w.r.t. SENSOR_STATUS not found. Followings are the variables:\n\t" + "\n\t".join([k.Name() for k in keys]))

        # first calculate the sensor status partial sensitivity of the response function
        d_j_d_sensor_status= physical_gradient_variable_container_expression_map[KratosDT.SENSOR_STATUS]

        # now calculate the total sensitivities of sensor_status w.r.t. phi
        d_sensor_status_d_phi = KratosOA.ControlUtils.SigmoidalProjectionUtils.CalculateForwardProjectionGradient(self.phi, [0, 1], [0, 1], self.beta, 1)

        # now compute response function total sensitivity w.r.t. phi
        d_j_d_phi = d_j_d_sensor_status * d_sensor_status_d_phi

        return self.filter.BackwardFilterField(d_j_d_phi)

    def Update(self, control_field: ContainerExpressionTypes) -> bool:
        if not IsSameContainerExpression(control_field, self.GetEmptyField()):
            raise RuntimeError(f"Updates for the required element container not found for control \"{self.GetName()}\". [ required model part name: {self.model_part.FullName()}, given model part name: {control_field.GetModelPart().FullName()} ]")

        if Kratos.Expression.Utils.NormL2(self.phi - control_field) > 1e-15:
            self.phi = self.filter.ForwardFilterField(control_field)
            self.un_buffered_data.SetValue("filtered_phi", self.phi.Clone(), overwrite=True)

            self.__UpdateBeta()

            self.__ApplySensorStatus()

    def __UpdateBeta(self) -> None:
        if self.beta_adaptive:
            step = self.optimization_problem.GetStep()
            if step % self.beta_update_period == 0 and self.beta_computed_step != step:
                self.beta_computed_step = step
                self.beta = min(self.beta * self.beta_increase_frac, self.beta_max)
                Kratos.Logger.PrintInfo(f"::{self.GetName()}::", f"Increased beta to {self.beta}.")

    def __ApplySensorStatus(self) -> None:
        sensor_status =  KratosOA.ControlUtils.SigmoidalProjectionUtils.ProjectForward(self.phi, [0, 1], [0, 1], self.beta, 1)
        Kratos.Expression.VariableExpressionIO.Write(sensor_status, KratosDT.SENSOR_STATUS, False)
        self.un_buffered_data.SetValue("SENSOR_STATUS", sensor_status.Clone(), overwrite=True)

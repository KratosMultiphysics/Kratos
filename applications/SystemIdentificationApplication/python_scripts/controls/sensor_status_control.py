import typing
import KratosMultiphysics as Kratos
import KratosMultiphysics.SystemIdentificationApplication as KratosSI
from KratosMultiphysics.OptimizationApplication.controls.control import Control
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import ContainerExpressionTypes
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import SupportedSensitivityFieldVariableTypes
from KratosMultiphysics.OptimizationApplication.utilities.model_part_utilities import ModelPartOperation
from KratosMultiphysics.OptimizationApplication.utilities.logger_utilities import TimeLogger
from KratosMultiphysics.OptimizationApplication.utilities.helper_utilities import IsSameContainerExpression
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView
from KratosMultiphysics.OptimizationApplication.utilities.opt_projection import CreateProjection
from KratosMultiphysics.SystemIdentificationApplication.utilities.sensor_utils import UpdateSensorStatusControlUpdaters

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
            "output_all_fields"          : false,
            "projection_settings"        : {}
        }""")
        parameters.ValidateAndAssignDefaults(default_settings)

        self.output_all_fields = parameters["output_all_fields"].GetBool()

        controlled_model_names_parts = parameters["controlled_model_part_names"].GetStringArray()
        if len(controlled_model_names_parts) == 0:
            raise RuntimeError(f"No model parts are provided for SensorStatusControl. [ control name = \"{self.GetName()}\"]")
        self.model_part_operation = ModelPartOperation(self.model, ModelPartOperation.OperationType.UNION, f"control_{self.GetName()}", controlled_model_names_parts, False)
        self.model_part: 'typing.Optional[Kratos.ModelPart]' = None

        self.projection = CreateProjection(parameters["projection_settings"], optimization_problem)

    def Initialize(self) -> None:
        self.model_part = self.model_part_operation.GetModelPart()

        ComponentDataView(self, self.optimization_problem).SetDataBuffer(1)
        self.un_buffered_data = ComponentDataView(self, self.optimization_problem).GetUnBufferedData()

        self.projection.SetProjectionSpaces([0, 1], [0, 1])

        sensor_status = Kratos.Expression.NodalExpression(self.model_part)
        Kratos.Expression.VariableExpressionIO.Read(sensor_status, KratosSI.SENSOR_STATUS, False)

        # project backward the uniform physical control field and assign it to the control field
        self.physical_phi_field = self.projection.ProjectBackward(sensor_status)

        # now update the physical thickness field
        self._UpdateAndOutputFields(self.GetEmptyField())

    def Check(self) -> None:
        pass

    def Finalize(self) -> None:
        pass

    def GetPhysicalKratosVariables(self) -> 'list[SupportedSensitivityFieldVariableTypes]':
        return [KratosSI.SENSOR_STATUS]

    def GetEmptyField(self) -> ContainerExpressionTypes:
        field = Kratos.Expression.NodalExpression(self.model_part)
        Kratos.Expression.LiteralExpressionIO.SetData(field, 0.0)
        return field

    def GetControlField(self) -> ContainerExpressionTypes:
        return self.physical_phi_field

    def MapGradient(self, physical_gradient_variable_container_expression_map: 'dict[SupportedSensitivityFieldVariableTypes, ContainerExpressionTypes]') -> ContainerExpressionTypes:
        keys = physical_gradient_variable_container_expression_map.keys()
        if len(keys) != 1:
            raise RuntimeError(f"Not provided required gradient fields for control \"{self.GetName()}\". Following are the variables:\n\t" + "\n\t".join([k.Name() for k in keys]))
        if  KratosSI.SENSOR_STATUS not in keys:
            raise RuntimeError(f"The required gradient for control \"{self.GetName()}\" w.r.t. SENSOR_STATUS not found. Followings are the variables:\n\t" + "\n\t".join([k.Name() for k in keys]))

        # first calculate the sensor status partial sensitivity of the response function
        physical_gradient = physical_gradient_variable_container_expression_map[KratosSI.SENSOR_STATUS]

        # multiply the physical sensitivity field with projection derivatives
        return physical_gradient * self.projection_derivative_field

    def Update(self, new_control_field: ContainerExpressionTypes) -> bool:
        if not IsSameContainerExpression(new_control_field, self.GetEmptyField()):
            raise RuntimeError(f"Updates for the required element container not found for control \"{self.GetName()}\". [ required model part name: {self.model_part.FullName()}, given model part name: {new_control_field.GetModelPart().FullName()} ]")

        update = new_control_field - self.physical_phi_field
        if Kratos.Expression.Utils.NormL2(update) > 1e-15:
            with TimeLogger(self.__class__.__name__, f"Updating {self.GetName()}...", f"Finished updating of {self.GetName()}.",False):
                # now update the physical field
                self._UpdateAndOutputFields(update)

                self.projection.Update()
                return True

        self.projection.Update()
        return False

    def _UpdateAndOutputFields(self, physical_phi_update: ContainerExpressionTypes) -> None:
        self.physical_phi_field = Kratos.Expression.Utils.Collapse(self.physical_phi_field + physical_phi_update)

        # project forward the filtered thickness field
        projected_sensor_field = self.projection.ProjectForward(self.physical_phi_field)

        # now update physical field
        Kratos.Expression.VariableExpressionIO.Write(projected_sensor_field, KratosSI.SENSOR_STATUS, False)

        # now update the sensor status control updaters
        UpdateSensorStatusControlUpdaters(self.optimization_problem)

        # compute and stroe projection derivatives for consistent filtering of the sensitivities
        self.projection_derivative_field = self.projection.ForwardProjectionGradient(self.physical_phi_field)

        # now output the fields
        un_buffered_data = ComponentDataView(self, self.optimization_problem).GetUnBufferedData()
        un_buffered_data.SetValue("sensor_status", projected_sensor_field.Clone(),overwrite=True)
        if self.output_all_fields:
            un_buffered_data.SetValue("sensor_status_physical_phi", self.physical_phi_field.Clone(), overwrite=True)
            un_buffered_data.SetValue("sensor_status_physical_phi_update", physical_phi_update.Clone(), overwrite=True)
            un_buffered_data.SetValue("sensor_status_projection_derivative", self.projection_derivative_field.Clone(), overwrite=True)

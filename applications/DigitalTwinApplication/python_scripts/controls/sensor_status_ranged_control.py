import typing
import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
import KratosMultiphysics.DigitalTwinApplication as KratosDT
from KratosMultiphysics.OptimizationApplication.controls.control import Control
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import ContainerExpressionTypes
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import SupportedSensitivityFieldVariableTypes
from KratosMultiphysics.OptimizationApplication.utilities.model_part_utilities import ModelPartOperation
from KratosMultiphysics.OptimizationApplication.utilities.logger_utilities import TimeLogger
from KratosMultiphysics.OptimizationApplication.utilities.helper_utilities import IsSameContainerExpression
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView

def Factory(model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem) -> Control:
    if not parameters.Has("name"):
        raise RuntimeError(f"SensorStatusRangedControl instantiation requires a \"name\" in parameters [ parameters = {parameters}].")
    if not parameters.Has("settings"):
        raise RuntimeError(f"SensorStatusRangedControl instantiation requires a \"settings\" in parameters [ parameters = {parameters}].")
    return SensorStatusRangedControl(parameters["name"].GetString(), model, parameters["settings"], optimization_problem)

class SensorStatusRangedControl(Control):
    def __init__(self, name: str, model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem) -> None:
        super().__init__(name)

        self.parameters = parameters
        self.model = model
        self.optimization_problem = optimization_problem

        default_settings = Kratos.Parameters("""{
            "controlled_model_part_names": [""],
            "mask_expression_name"       : "YOUNG_MODULUS_SENSITIVITY_mask",
            "penalty_factor"             : 3,
            "output_all_fields"          : false,
            "beta_settings": [
                {
                    "end_iteration": 5e+3,
                    "beta"         : 30.0
                },
                {
                    "end_iteration": 1e+4,
                    "beta"         : 60.0
                },
                {
                    "end_iteration": 1.5e+4,
                    "beta"         : 300.0
                }
            ]
        }""")
        parameters.ValidateAndAssignDefaults(default_settings)

        self.controlled_physical_variables = [KratosDT.SENSOR_STATUS]
        self.output_all_fields = parameters["output_all_fields"].GetBool()
        self.penalty_factor = parameters["penalty_factor"].GetInt()
        self.mask_expression_name = parameters["mask_expression_name"].GetString()

        # beta settings
        beta_settings = parameters["beta_settings"]
        self.beta_values: 'list[tuple[int, float]]' = []
        for values in beta_settings.values():
            self.beta_values.append((values["end_iteration"].GetInt(), values["beta"].GetDouble()))

        self.beta = self.beta_values[0][1]

        controlled_model_names_parts = parameters["controlled_model_part_names"].GetStringArray()
        if len(controlled_model_names_parts) == 0:
            raise RuntimeError(f"No model parts are provided for SensorStatusRangedControl. [ control name = \"{self.GetName()}\"]")
        self.model_part_operation = ModelPartOperation(self.model, ModelPartOperation.OperationType.UNION, f"control_{self.GetName()}", controlled_model_names_parts, False)
        self.model_part: 'typing.Optional[Kratos.ModelPart]' = None

    def Initialize(self) -> None:
        self.model_part = self.model_part_operation.GetModelPart()

        self.un_buffered_data = ComponentDataView(self, self.optimization_problem).GetUnBufferedData()

        sensor_status = Kratos.Expression.NodalExpression(self.model_part)
        Kratos.Expression.VariableExpressionIO.Read(sensor_status, KratosDT.SENSOR_STATUS, False)

        # project backward the uniform physical control field and assign it to the control field
        self.physical_phi_field = KratosOA.ControlUtils.SigmoidalProjectionUtils.ProjectBackward(
                                                sensor_status,
                                                [0, 1],
                                                [0, 1],
                                                self.beta,
                                                self.penalty_factor)

        sensor_data = ComponentDataView("sensor", self.optimization_problem)
        sensors_list: 'list[KratosDT.Sensors.Sensor]' = sensor_data.GetUnBufferedData().GetValue("list_of_sensors")
        if len(sensors_list) == 0:
            raise RuntimeError("Empty sensors list.")

        first_sensor = sensors_list[0]
        if self.mask_expression_name in first_sensor.GetNodalExpressionsMap().keys():
            self.sensor_mask_status = KratosDT.MaskUtils.SensorNodalMaskStatus(self.model_part, [sensor.GetNodalExpression(self.mask_expression_name) for sensor in sensors_list])
            self.sensor_mask_status_kd_tree = KratosDT.MaskUtils.SensorNodalMaskStatusKDTree(self.sensor_mask_status, 4)
        elif self.mask_expression_name in first_sensor.GetConditionExpressionsMap().keys():
            self.sensor_mask_status = KratosDT.MaskUtils.SensorConditionMaskStatus(self.model_part, [sensor.GetConditionExpression(self.mask_expression_name) for sensor in sensors_list])
            self.sensor_mask_status_kd_tree = KratosDT.MaskUtils.SensorConditionMaskStatusKDTree(self.sensor_mask_status, 4)
        elif self.mask_expression_name in first_sensor.GetElementExpressionsMap().keys():
            self.sensor_mask_status = KratosDT.MaskUtils.SensorElementMaskStatus(self.model_part, [sensor.GetElementExpression(self.mask_expression_name) for sensor in sensors_list])
            self.sensor_mask_status_kd_tree = KratosDT.MaskUtils.SensorElementMaskStatusKDTree(self.sensor_mask_status, 4)
        else:
            raise RuntimeError(f"The sensor mask expression \"{self.mask_expression_name}\" not found in the list of sensors.")

        sensor_data.GetUnBufferedData().SetValue("sensor_mask_status", self.sensor_mask_status)
        sensor_data.GetUnBufferedData().SetValue("sensor_mask_status_kd_tree", self.sensor_mask_status_kd_tree)

        # now update the physical thickness field
        self._UpdateAndOutputFields(self.GetEmptyField())

    def Check(self) -> None:
        pass

    def Finalize(self) -> None:
        pass

    def GetPhysicalKratosVariables(self) -> 'list[SupportedSensitivityFieldVariableTypes]':
        return self.controlled_physical_variables

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
        if  KratosDT.SENSOR_STATUS not in keys:
            raise RuntimeError(f"The required gradient for control \"{self.GetName()}\" w.r.t. SENSOR_STATUS not found. Followings are the variables:\n\t" + "\n\t".join([k.Name() for k in keys]))

        # first calculate the sensor status partial sensitivity of the response function
        physical_gradient = physical_gradient_variable_container_expression_map[KratosDT.SENSOR_STATUS]

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

                self.__UpdateBeta()
                return True

        self.__UpdateBeta()
        return False

    def __UpdateBeta(self) -> None:
        step = self.optimization_problem.GetStep()

        for end_iteration, beta in self.beta_values:
            if step < end_iteration:
                if self.beta != beta:
                    self.beta = beta
                    # compute and stroe projection derivatives for consistent filtering of the sensitivities
                    self.projection_derivative_field = KratosOA.ControlUtils.SigmoidalProjectionUtils.CalculateForwardProjectionGradient(
                                                            self.physical_phi_field,
                                                            [0, 1],
                                                            [0, 1],
                                                            self.beta,
                                                            self.penalty_factor)
                    Kratos.Logger.PrintInfo(self.__class__.__name__, f"Changed beta to {self.beta}...")
                break

    def _UpdateAndOutputFields(self, physical_phi_update: ContainerExpressionTypes) -> None:
        self.physical_phi_field = Kratos.Expression.Utils.Collapse(self.physical_phi_field + physical_phi_update)

        # project forward the filtered thickness field
        projected_sensor_field = KratosOA.ControlUtils.SigmoidalProjectionUtils.ProjectForward(
                                                self.physical_phi_field,
                                                [0, 1],
                                                [0, 1],
                                                self.beta,
                                                self.penalty_factor)
        # now update physical field
        Kratos.Expression.VariableExpressionIO.Write(projected_sensor_field, KratosDT.SENSOR_STATUS, False)
        list_of_sensors: 'list[KratosDT.Sensors.Sensor]' = ComponentDataView("sensor", self.optimization_problem).GetUnBufferedData().GetValue("list_of_sensors")
        for i, node in enumerate(self.model_part.Nodes):
            list_of_sensors[i].SetValue(KratosDT.SENSOR_STATUS, node.GetValue(KratosDT.SENSOR_STATUS))

        self.sensor_mask_status.Update()
        self.sensor_mask_status_kd_tree.Update()

        # compute and stroe projection derivatives for consistent filtering of the sensitivities
        self.projection_derivative_field = KratosOA.ControlUtils.SigmoidalProjectionUtils.CalculateForwardProjectionGradient(
                                                self.physical_phi_field,
                                                [0, 1],
                                                [0, 1],
                                                self.beta,
                                                self.penalty_factor)

        # now output the fields
        un_buffered_data = ComponentDataView(self, self.optimization_problem).GetUnBufferedData()
        un_buffered_data.SetValue("sensor_status", projected_sensor_field.Clone(),overwrite=True)
        if self.output_all_fields:
            un_buffered_data.SetValue("sensor_status_physical_phi", self.physical_phi_field.Clone(), overwrite=True)
            un_buffered_data.SetValue("sensor_status_physical_phi_update", physical_phi_update.Clone(), overwrite=True)
            un_buffered_data.SetValue("sensor_status_projection_derivative", self.projection_derivative_field.Clone(), overwrite=True)

from pathlib import Path
import typing
import KratosMultiphysics as Kratos
import KratosMultiphysics.SystemIdentificationApplication as KratosSI

from KratosMultiphysics.SystemIdentificationApplication.utilities.sensor_utils import CreateSensors
from KratosMultiphysics.SystemIdentificationApplication.utilities.sensor_utils import PrintSensorListToCSV
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView
from KratosMultiphysics.SystemIdentificationApplication.utilities.sensor_utils import GetSensors

def Factory(arg1: Kratos.Parameters, arg2: Kratos.Model, optimization_problem: 'typing.Optional[OptimizationProblem]' = None):
    if isinstance(arg1, Kratos.Parameters) and isinstance(arg2, Kratos.Model):
        # this is to use SensorOutputProcess as a normal kratos process
        return SensorOutputProcess(arg2, arg1["Parameters"], optimization_problem)
    elif isinstance(arg1, Kratos.Model) and isinstance(arg2, Kratos.Parameters):
        # this is to use SensorOutputProcess as a process for optimization application where optimization_problem cannot be none
        return SensorOutputProcess(arg1, arg2["settings"], optimization_problem)
    else:
        raise RuntimeError("Argument mismatch.")

class SensorOutputProcess(Kratos.OutputProcess):
    def __init__(self, model: Kratos.Model, settings: Kratos.Parameters, optimization_problem: 'typing.Optional[OptimizationProblem]' = None):
        Kratos.OutputProcess.__init__(self)

        default_settings = Kratos.Parameters("""{
            "model_part_name" : "PLEASE_PROVIDE_A_MODEL_PART_NAME",
            "sensor_model_part_name": "SensorModelPart",
            "sensor_group_name": "sensors",
            "output_file_name": "",
            "write_only_active_sensors": false,
            "properties_list": [
                "type",
                "name",
                "location",
                "value"
            ],
            "list_of_sensors": []
        }""")

        if optimization_problem is not None:
            default_settings.RemoveValue("list_of_sensors")
            default_settings.AddEmptyValue("response_name")
            default_settings["response_name"].SetString("")
        settings.ValidateAndAssignDefaults(default_settings)

        self.model_part = model[settings["model_part_name"].GetString()]
        self.properties_list: 'list[str]' = settings["properties_list"].GetStringArray()
        self.output_file_name = Path(settings["output_file_name"].GetString())
        self.optimization_problem = optimization_problem
        self.write_only_active_sensors = settings["write_only_active_sensors"].GetBool()
        self.sensor_group_name = settings["sensor_group_name"].GetString()

        self.list_of_sensors: 'list[KratosSI.Sensors.Sensor]' = []
        if optimization_problem is None:
            self.initialized_sensors = True
            self.compute_sensor_value = True
            sensor_model_part_name = settings["sensor_model_part_name"].GetString()
            if not model.HasModelPart(sensor_model_part_name):
                model.CreateModelPart(sensor_model_part_name)
            self.list_of_sensors = CreateSensors(model[sensor_model_part_name], self.model_part, settings["list_of_sensors"].values())
        else:
            self.initialized_sensors = False
            self.compute_sensor_value = False
            self.response_name = settings["response_name"].GetString()

    def PrintOutput(self) -> None:
        s_name = str(self.output_file_name)

        if not self.initialized_sensors:
            self.list_of_sensors = GetSensors(ComponentDataView(self.sensor_group_name, self.optimization_problem))
            self.initialized_sensors = True
            s_name = s_name.replace("<step>", f"{self.optimization_problem.GetStep():06d}")
        else:
            s_name = s_name.replace("<step>", f"{self.model_part.ProcessInfo[Kratos.STEP]:06d}")

        # calculate and assign values for each sensor
        if self.compute_sensor_value:
            for sensor in self.list_of_sensors:
                sensor.SetSensorValue(sensor.CalculateValue(self.model_part))

        if self.write_only_active_sensors:
            PrintSensorListToCSV(Path(s_name), [sensor for sensor in self.list_of_sensors if sensor.IsActive()], self.properties_list)
        else:
            PrintSensorListToCSV(Path(s_name), self.list_of_sensors, self.properties_list)






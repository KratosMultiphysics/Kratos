from pathlib import Path
import typing
import KratosMultiphysics as Kratos
import KratosMultiphysics.SystemIdentificationApplication as KratosSI

from KratosMultiphysics.SystemIdentificationApplication.utilities.sensor_utils import GetSensors
from KratosMultiphysics.SystemIdentificationApplication.utilities.sensor_utils import PrintSensorListToCSV
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView
from KratosMultiphysics.OptimizationApplication.responses.response_function import ResponseFunction

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
            "output_file_name": "",
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

        self.list_of_sensors: 'list[KratosSI.Sensors.Sensor]' = []
        if optimization_problem is None:
            self.initialized_sensors = True
            self.compute_sensor_value = True
            self.list_of_sensors = GetSensors(self.model_part, settings["list_of_sensors"].values())
        else:
            self.initialized_sensors = False
            self.compute_sensor_value = False
            self.response_name = settings["response_name"].GetString()

    def PrintOutput(self) -> None:
        s_name = str(self.output_file_name)

        if not self.initialized_sensors:
            self.list_of_sensors = ComponentDataView(self.optimization_problem.GetComponent(self.response_name, ResponseFunction), self.optimization_problem).GetUnBufferedData().GetValue("sensors")
            self.initialized_sensors = True
            s_name = s_name.replace("<step>", f"{self.optimization_problem.GetStep():06d}")
        else:
            s_name = s_name.replace("<step>", f"{self.model_part.ProcessInfo[Kratos.STEP]:06d}")

        # calculate and assign values for each sensor
        if self.compute_sensor_value:
            for sensor in self.list_of_sensors:
                if "eigenvector" in sensor.GetName():            
                    print(sensor.CalculateValueVector(self.model_part))
                    sensor.SetSensorValueVector(sensor.CalculateValueVector(self.model_part))
                else:
                    sensor.SetSensorValue(sensor.CalculateValue(self.model_part))

        PrintSensorListToCSV(Path(s_name), self.list_of_sensors, self.properties_list)






from pathlib import Path
import KratosMultiphysics as Kratos
import KratosMultiphysics.DigitalTwinApplication as KratosDT

from KratosMultiphysics.DigitalTwinApplication.utilities.sensor_utils import GetSensors
from KratosMultiphysics.DigitalTwinApplication.utilities.sensor_utils import PrintSensorViewsListToCSV

def Factory(settings: Kratos.Parameters, model: Kratos.Model):
    if (not isinstance(model, Kratos.Model)):
        raise RuntimeError("expected input shall be a model object, encapsulating a json string")
    if (not isinstance(settings, Kratos.Parameters)):
        raise RuntimeError("expected input shall be a Parameters object, encapsulating a json string")
    return SensorOutputProcess(model, settings["Parameters"])

class SensorOutputProcess(Kratos.OutputProcess):
    def __init__(self, model: Kratos.Model, settings: Kratos.Parameters):
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
        settings.ValidateAndAssignDefaults(default_settings)

        self.model_part = model[settings["model_part_name"].GetString()]
        self.list_of_sensors: 'list[KratosDT.Sensors.Sensor]' = GetSensors(self.model_part, settings["list_of_sensors"].values())
        self.properties_list: 'list[str]' = settings["properties_list"].GetStringArray()
        self.output_file_name = Path(settings["output_file_name"].GetString())

    def PrintOutput(self) -> None:
        s_name = str(self.output_file_name)
        s_name = s_name.replace("<step>", str(self.model_part.ProcessInfo[Kratos.STEP]))
        s_name = s_name.replace("<time>", str(self.model_part.ProcessInfo[Kratos.TIME]))
        # calculate and assign values for each sensor
        for sensor in self.list_of_sensors:
            sensor.SetSensorValue(sensor.CalculateValue(self.model_part))
        PrintSensorViewsListToCSV(Path(s_name), self.list_of_sensors, self.properties_list)






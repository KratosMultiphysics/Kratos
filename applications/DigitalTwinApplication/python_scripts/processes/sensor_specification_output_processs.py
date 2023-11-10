from pathlib import Path
import KratosMultiphysics as Kratos
import KratosMultiphysics.DigitalTwinApplication as KratosDT

from KratosMultiphysics.DigitalTwinApplication.utilities.sensor_specification_utils import GetSpecifications

def Factory(settings: Kratos.Parameters, model: Kratos.Model):
    if (not isinstance(model, Kratos.Model)):
        raise Exception(
            "expected input shall be a model object, encapsulating a json string"
        )
    if (not isinstance(settings, Kratos.Parameters)):
        raise Exception(
            "expected input shall be a Parameters object, encapsulating a json string"
        )
    return SensorSpecificationOutputProcess(model, settings["Parameters"])


class SensorSpecificationOutputProcess(Kratos.OutputProcess):
    def __init__(self, model: Kratos.Model, settings: Kratos.Parameters):
        Kratos.OutputProcess.__init__(self)

        default_settings = Kratos.Parameters("""{
            "model_part_name" : "PLEASE_PROVIDE_A_MODEL_PART_NAME",
            "output_file_name": "",
            "sensor_settings" : [
                {
                    "sensor_type"           : "PLEASE_SPECIFY_SENSOR_TYPE",
                    "base_settings"         : {},
                    "list_of_specifications": []
                }
            ]
        }""")
        settings.ValidateAndAssignDefaults(default_settings)

        self.model_part = model[settings["model_part_name"].GetString()]

        # create the list of specifications
        self.dict_of_sensor_specifications: 'dict[KratosDT.Sensors.AdjointSensor, list[KratosDT.Sensors.SensorSpecification]]' = {}
        sensor_setting: Kratos.Parameters
        for sensor_setting in settings["sensor_settings"].values():
            sensor_setting.ValidateAndAssignDefaults(default_settings["sensor_settings"][0])

            sensor_type = sensor_setting["sensor_type"].GetString()

            if sensor_type == "adjoint_displacement_sensor":
                adjoint_sensor = KratosDT.Sensors.AdjointDisplacementSensor(model, sensor_setting["base_settings"])
                self.dict_of_sensor_specifications[adjoint_sensor] = GetSpecifications(self.model_part, sensor_setting["list_of_specifications"].values())
            else:
                raise RuntimeError(f"Unsupported sensor_type = \"{sensor_type}\".")

        self.output_file_name = Path(settings["output_file_name"].GetString())

    def PrintOutput(self) -> None:
        self.output_file_name.parent.mkdir(exist_ok=True, parents=True)
        with open(str(self.output_file_name), "w") as file_output:
            file_output.write("#; type; value; location_x; location_y; location_z\n")
            for adjoint_sensor, list_of_specifications in self.dict_of_sensor_specifications.items():
                for spec in list_of_specifications:
                    adjoint_sensor.SetSensorSpecification(spec)
                    spec.SetSensorValue(adjoint_sensor.CalculateValue(self.model_part))
                    loc = spec.GetLocation()
                    file_output.write(f"{spec.Id}; {spec.GetType()}; {spec.GetSensorValue()}; {loc[0]}; {loc[1]}; {loc[2]}\n")






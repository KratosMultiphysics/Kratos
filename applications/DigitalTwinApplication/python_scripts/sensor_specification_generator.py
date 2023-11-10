import typing
import KratosMultiphysics as Kratos
import json

class SensorSpecificationGenerator:
    def __init__(self, model: Kratos.Model, parameters: Kratos.Parameters) -> None:
        default_settings = Kratos.Parameters("""{
            "model_part_name"       : "PLEASE_SPECIFY_MODEL_PART_NAME",
            "mdpa_file_name"        : "PlEASE_SPECIFY_FILE_NAME",
            "output_file_name"      : "sensor_specification_data.json",
            "specification_settings": {
                "model_part_name": "PLEASE_SPECIFY_MODEL_PART_NAME",
                "container_type" : "elements",
                "common_settings": {}
            }
        }""")

        self.parameters = parameters
        self.parameters.ValidateAndAssignDefaults(default_settings)
        self.parameters["specification_settings"].ValidateAndAssignDefaults(default_settings["specification_settings"])

        self.model = model
        self.model_part = self.model.CreateModelPart(self.parameters["model_part_name"].GetString())

    def Run(self) -> None:
        Kratos.ModelPartIO(self.parameters["mdpa_file_name"].GetString(), Kratos.ModelPartIO.READ | Kratos.ModelPartIO.MESH_ONLY).ReadModelPart(self.model_part)

        spec_settings = self.parameters["specification_settings"]
        model_part = self.model[spec_settings["model_part_name"].GetString()]

        container_type = spec_settings["container_type"].GetString()
        container: 'typing.Union[Kratos.ElementsArray, Kratos.ConditionsArray]'
        if container_type == "elements":
            container = model_part.Elements
        elif container_type == "conditions":
            container = model_part.Conditions
        else:
            raise RuntimeError(f"Unsupported container_type = \"{container_type}\".")

        specifications = {
            "list_of_specifications": []
        }

        prototype_spec = {}
        SensorSpecificationGenerator.__FillSpecPrototype(prototype_spec, spec_settings["common_settings"])

        for entity in container:
            spec = dict(prototype_spec)
            spec["id"] = entity.Id
            loc = entity.GetGeometry().Center()
            spec["location"] = [loc[0], loc[1], loc[2]]
            specifications["list_of_specifications"].append(spec)

        with open(self.parameters["output_file_name"].GetString(), "w") as file_output:
            file_output.write(json.dumps(specifications, indent=4))

    @staticmethod
    def __FillSpecPrototype(prototype, params: Kratos.Parameters):
        value: Kratos.Parameters
        for key, value in params.items():
            if value.IsString():
                prototype[key] = value.GetString()
            elif value.IsInt():
                prototype[key] = value.GetInt()
            elif value.IsDouble():
                prototype[key] = value.GetDouble()
            elif value.IsBool():
                prototype[key] = value.GetBool()
            elif value.IsVector():
                vec = value.GetVector()
                list_of_values = [v for v in vec]
                prototype[key] = list_of_values
            elif value.IsSubParameter():
                prototype[key] = {}
                SensorSpecificationGenerator.__FillSpecPrototype(prototype[key], value)
            else:
                raise RuntimeError(f"Unsupported value at key {key} [ value = {value} ]")






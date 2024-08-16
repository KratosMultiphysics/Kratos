import typing
import json

import KratosMultiphysics as Kratos

class SensorGeneratorAnalysis:
    def __init__(self, model: Kratos.Model, parameters: Kratos.Parameters) -> None:
        default_settings = Kratos.Parameters("""{
            "model_part_name" : "PLEASE_SPECIFY_MODEL_PART_NAME",
            "mdpa_file_name"  : "PlEASE_SPECIFY_FILE_NAME",
            "output_file_name": "sensor_data.json",
            "sensor_groups"   : [
                {
                    "generation_method": "mesh_based",
                    "container_type"   : "elements",
                    "location"         : "center",
                    "sensor_settings"  : {}
                }
            ]
        }""")

        self.parameters = parameters
        self.parameters.ValidateAndAssignDefaults(default_settings)
        self.model = model
        self.model_part = self.model.CreateModelPart(self.parameters["model_part_name"].GetString())

    def Run(self) -> None:
        Kratos.ModelPartIO(self.parameters["mdpa_file_name"].GetString(), Kratos.ModelPartIO.READ | Kratos.ModelPartIO.MESH_ONLY).ReadModelPart(self.model_part)

        sensors_dict = {"list_of_sensors": []}
        for sensor_group in self.parameters["sensor_groups"].values():
            generation_method = sensor_group["generation_method"].GetString()
            if generation_method == "mesh_based":
                self.__GenerateMeshBased(sensors_dict, sensor_group)
            else:
                raise RuntimeError("Unsupported generation method")

        with open(self.parameters["output_file_name"].GetString(), "w") as file_output:
            file_output.write(json.dumps(sensors_dict, indent=4))

    def __GenerateMeshBased(self, sensors_dict, sensor_group_params: Kratos.Parameters):
        container_type = sensor_group_params["container_type"].GetString()
        container: 'typing.Union[Kratos.ElementsArray, Kratos.ConditionsArray]'
        if container_type == "nodes":
            container = self.model_part.Nodes
        elif container_type == "elements":
            container = self.model_part.Elements
        elif container_type == "conditions":
            container = self.model_part.Conditions
        else:
            raise RuntimeError(f"Unsupported container_type = \"{container_type}\".")

        sensor_params = sensor_group_params["sensor_settings"]
        if not sensor_params.Has("location"):
            sensor_params.AddVector("location", Kratos.Vector())

        if container_type == "nodes":
            location_retrieval = lambda x: [x.X, x.Y, x.Z]
        elif location in ["elements", "conditions"]:
            if location == "center":
                location = sensor_group_params["location"].GetString()  
                location_retrieval = lambda x: x.GetGeometry().Center()
            else:
                raise RuntimeError(f"Unsupported location = \"{location}\".")

        for entity in container:
            current_params = sensor_params.Clone()
            current_params["location"].SetVector(location_retrieval(entity))
            json_dict = json.loads(current_params.WriteJsonString().replace("<ENTITY_ID>", str(entity.Id)))
            sensors_dict["list_of_sensors"].append(json_dict)

import typing
import json

import KratosMultiphysics as Kratos
import KratosMultiphysics.SystemIdentificationApplication as KratosSI

class SensorGeneratorAnalysis:
    def __init__(self, model: Kratos.Model, parameters: Kratos.Parameters) -> None:
        default_settings = Kratos.Parameters("""{
            "model_part_name"     : "PLEASE_SPECIFY_MODEL_PART_NAME",
            "mdpa_file_name"      : "PlEASE_SPECIFY_FILE_NAME",
            "output_file_name"    : "sensor_data.json",
            "csv_output_file_name": "sensor_data.csv",
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

        model_part_name = self.parameters["model_part_name"].GetString()
        root_model_part = model_part_name.split(".")[0]
        self.root_model_part = self.model.CreateModelPart(root_model_part)

    def Run(self) -> None:
        Kratos.ModelPartIO(self.parameters["mdpa_file_name"].GetString(), Kratos.ModelPartIO.READ | Kratos.ModelPartIO.MESH_ONLY).ReadModelPart(self.root_model_part)
        self.model_part = self.model[self.parameters["model_part_name"].GetString()]

        sensors_dict = {"list_of_sensors": []}
        for sensor_group in self.parameters["sensor_groups"].values():
            generation_method = sensor_group["generation_method"].GetString()
            if generation_method == "mesh_based":
                self.__GenerateMeshBased(sensors_dict, sensor_group)
            elif generation_method == "bounding_surface":
                self.__GenerationBoundingSurfaceBased(sensors_dict, sensor_group)
            else:
                raise RuntimeError("Unsupported generation method")

        with open(self.parameters["output_file_name"].GetString(), "w") as file_output:
            file_output.write(json.dumps(sensors_dict, indent=4))

        with open(self.parameters["csv_output_file_name"].GetString(), "w") as file_output:
            file_output.write("#,name,location_0,location_1,location_2\n")
            for i, sensor_info in enumerate(sensors_dict["list_of_sensors"]):
                name = sensor_info["name"]
                location = sensor_info["location"]
                file_output.write(f"{i+1},{name},{location[0]},{location[1]},{location[2]}\n")

    def __GenerationBoundingSurfaceBased(self, sensors_dict, sensor_group_params: Kratos.Parameters):
        defaults = Kratos.Parameters("""{
                "generation_method"        : "bounding_surface",
                "bounding_surface_corner_1": [0.0, 0.0, 0.0],
                "bounding_surface_corner_2": [0.0, 0.0, 0.0],
                "number_of_sensors"        : [1, 1, 1],
                "sensor_settings"  : {}
        }""")
        sensor_group_params.ValidateAndAssignDefaults(defaults)

        point_locator = Kratos.BruteForcePointLocator(self.model_part)

        corner_1 = sensor_group_params["bounding_surface_corner_1"].GetVector()
        corner_2 = sensor_group_params["bounding_surface_corner_2"].GetVector()
        number_of_sensors = sensor_group_params["number_of_sensors"].GetVector()
        number_of_sensors = [int(v) for v in number_of_sensors]

        sensor_params = sensor_group_params["sensor_settings"]
        if not sensor_params.Has("location"):
            sensor_params.AddVector("location", Kratos.Vector())

        distance = corner_2 - corner_1
        index = 1
        for i_x in range(number_of_sensors[0]):
            if distance[0] > 0.0:
                x_coord = corner_1[0] + distance[0] * (i_x) / (number_of_sensors[0] - 1)
            else:
                x_coord = corner_1[0]
            for i_y in range(number_of_sensors[1]):
                if distance[1] > 0.0:
                    y_coord = corner_1[1] + distance[1] * (i_y) / (number_of_sensors[1] - 1)
                else:
                    y_coord = corner_1[1]
                for i_z in range(number_of_sensors[2]):
                    if distance[2] > 0.0:
                        z_coord = corner_1[2] + distance[2] * (i_z) / (number_of_sensors[2] - 1)
                    else:
                        z_coord = corner_1[2]
                    loc = Kratos.Point(x_coord, y_coord, z_coord)
                    shape_funcs = Kratos.Vector()
                    elem_id = point_locator.FindElement(loc, shape_funcs, Kratos.Configuration.Initial, 1e-8)
                    if elem_id != -1 and KratosSI.SensorUtils.IsPointInGeometry(loc, self.model_part.GetElement(elem_id).GetGeometry()):
                        current_params = sensor_params.Clone()
                        current_params["location"].SetVector(loc)
                        json_dict = json.loads(current_params.WriteJsonString().replace("<ENTITY_ID>", str(index)))
                        index += 1
                        sensors_dict["list_of_sensors"].append(json_dict)

    def __GenerateMeshBased(self, sensors_dict, sensor_group_params: Kratos.Parameters):
        defaults = Kratos.Parameters("""{
                "generation_method": "mesh_based",
                "container_type"   : "elements",
                "location"         : "center",
                "sensor_settings"  : {}
        }""")
        sensor_group_params.ValidateAndAssignDefaults(defaults)
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

        location = sensor_group_params["location"].GetString()
        if container_type == "nodes":
            location_retrieval = lambda x: [x.X, x.Y, x.Z]
        elif location == "center":
            location_retrieval = lambda x: x.GetGeometry().Center()
        else:
            raise RuntimeError(f"Unsupported location = \"{location}\".")

        for entity in container:
            current_params = sensor_params.Clone()
            current_params["location"].SetVector(location_retrieval(entity))
            json_dict = json.loads(current_params.WriteJsonString().replace("<ENTITY_ID>", str(entity.Id)))
            sensors_dict["list_of_sensors"].append(json_dict)

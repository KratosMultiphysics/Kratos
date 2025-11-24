import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.model_part_controllers.model_part_controller import ModelPartController

def Factory(model: Kratos.Model, parameters: Kratos.Parameters, _) -> ModelPartController:
    if not parameters.Has("settings"):
        raise RuntimeError(f"ExternalResponseFunctionModelPartController instantiation requires a \"settings\" in parameters [ parameters = {parameters}].")
    return ExternalResponseFunctionModelPartController(model, parameters["settings"])

class LocationAgnosticDesignPoints:
    def __init__(self, model_part: Kratos.ModelPart, parameters: Kratos.Parameters):
        default_settings = Kratos.Parameters("""{
            "type"                              : "location_agnostic_design_points",
            "number_of_design_points"           : 1
        }""")
        parameters.ValidateAndAssignDefaults(default_settings)
        self.model_part = model_part

        self.number_of_nodes = parameters["number_of_design_points"].GetInt()

    def Execute(self):
        for i in range(self.number_of_nodes):
            node: Kratos.Node = self.model_part.CreateNewNode(i + 1, i + 1, 0, 0)
            node.SetValue(KratosOA.CUSTOM_DESIGN_VARIABLE, 0.0)

class CSVBasedDesignPoints:
    def __init__(self, model_part: Kratos.ModelPart, parameters: Kratos.Parameters):
        default_settings = Kratos.Parameters("""{
            "type"                              : "csv_based_design_points",
            "csv_file_name"                     : "",
            "location_column_headers"           : [""]
        }""")
        parameters.ValidateAndAssignDefaults(default_settings)
        self.model_part = model_part

        self.csv_file_name = parameters["csv_file_name"].GetString()
        self.location_column_headers = parameters["location_column_headers"].GetStringArray()

        if len(self.location_column_headers) > 3:
            raise RuntimeError("Number of location column headers cannot exceed 3.")

    def Execute(self):
        with open(self.csv_file_name, "r") as file_input:
            lines = file_input.readlines()

        if len(lines) == 0:
            raise RuntimeError(f"The provided csv file = \"{self.csv_file_name}\" is empty.")

        headers = lines[0].split(",")
        location_column_indices = [-1, -1, -1]

        for j, location_header in enumerate(self.location_column_headers):
            for i, header in enumerate(headers):
                if header.strip() == location_header:
                    location_column_indices[j] = i
                    break

        for i, line in enumerate(lines[1:]):
            data = [float(v.strip()) for v in line[:-1].split(',')]
            location = [
                data[location_column_indices[0]] if location_column_indices[0] != -1 else 0.0,
                data[location_column_indices[1]] if location_column_indices[1] != -1 else 0.0,
                data[location_column_indices[2]] if location_column_indices[2] != -1 else 0.0
            ]
            node: Kratos.Node = self.model_part.CreateNewNode(i + 1, location[0], location[1], location[2])
            node.SetValue(KratosOA.CUSTOM_DESIGN_VARIABLE, 0.0)

class ExternalResponseFunctionModelPartController(ModelPartController):
    def __init__(self, model: Kratos.Model, parameters: Kratos.Parameters):
        default_settings = Kratos.Parameters("""{
            "model_part_name"       : "",
            "design_point_settings" : {}
        }""")

        parameters.ValidateAndAssignDefaults(default_settings)

        model_part_name = parameters["model_part_name"].GetString()
        if model_part_name == "":
            raise RuntimeError("Empty \"model_part_name\" is not allowed which is given with following parameters:\n" + str(parameters))

        self.model_part = model.CreateModelPart(model_part_name)

        design_point_settings = parameters["design_point_settings"]
        if not design_point_settings.Has("type"):
            raise RuntimeError(f"No design point creation method was specified. \n{design_point_settings}")

        design_point_type = design_point_settings["type"].GetString()
        if design_point_type == "location_agnostic_design_points":
            self.design_point_creation_method = LocationAgnosticDesignPoints(self.model_part, design_point_settings)
        elif design_point_type == "csv_based_design_points":
            self.design_point_creation_method = CSVBasedDesignPoints(self.model_part, design_point_settings)
        else:
            raise RuntimeError(f"Unsupported design point creation method type = \"{design_point_type}\". Followings are supported:\n\tlocation_agnostic_design_points\n\tcsv_based_design_points")

    def ImportModelPart(self) -> None:
        self.design_point_creation_method.Execute()

    def GetModelPart(self) -> Kratos.ModelPart:
        return self.model_part

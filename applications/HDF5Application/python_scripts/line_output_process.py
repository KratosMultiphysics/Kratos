import KratosMultiphysics
from KratosMultiphysics.HDF5Application.point_set_output_process import PointSetOutputProcess


def Factory(parameters: KratosMultiphysics.Parameters,
            model: KratosMultiphysics.Model):
    if not isinstance(parameters, KratosMultiphysics.Parameters):
        raise Exception("expecting input parameters of type KratosMultiphysics.Parameters, but got {}".format(type(parameters)))
    return LineOutputProcess(model, parameters["Parameters"])

class LineOutputProcess(KratosMultiphysics.OutputProcess):
    """Convenience process wrapping PointSetOutputProcess with points along a line."""
    @classmethod
    def GetDefaultParameters(cls) -> KratosMultiphysics.Parameters:
        return KratosMultiphysics.Parameters("""
            {
                "model_part_name": "PLEASE_SPECIFY_MODEL_PART_NAME",
                "file_settings": {
                    "file_name"        : "",
                    "time_format"      : "0.4f",
                    "file_access_mode" : "truncate",
                    "max_files_to_keep": "unlimited",
                    "echo_level"       :  0
                },
                "output_time_settings": {
                    "output_control_type": "step",
                    "output_interval"    : 1.0,
                    "interval"           : [0.0, "End"]
                },
                "line_output_settings": {
                    "prefix"              : "/VertexData/Coordinates",
                    "time_format"         : "0.4f",
                    "start_point"         : [0.0, 0.0, 0.0],
                    "end_point"           : [0.0, 0.0, 0.0],
                    "number_of_points"    : 2,
                    "search_configuration": "initial",
                    "search_tolerance"    : 1e-6,
                    "custom_attributes"   : {}
                },
                "nodal_solution_step_data_settings": {
                    "prefix"           : "/VertexData/VertexSolutionStepData/",
                    "list_of_variables": [],
                    "time_format"      : "0.4f",
                    "custom_attributes": {}
                },
                "nodal_data_value_settings": {
                    "prefix"           : "/VertexData/VertexDataValues/",
                    "list_of_variables": [],
                    "time_format"      : "0.4f",
                    "custom_attributes": {}
                }
            }""")

    def __init__(self, model: KratosMultiphysics.Model, parameters: KratosMultiphysics.Parameters):
        KratosMultiphysics.OutputProcess.__init__(self)
        self.wrapped_process = PointSetOutputProcess(model, self.__ConvertToPointSetOutputParameters(parameters))

    def IsOutputStep(self) -> bool:
        return self.wrapped_process.IsOutputStep()

    def PrintOutput(self) -> None:
        return self.wrapped_process.PrintOutput()

    @staticmethod
    def __ConvertToPointSetOutputParameters(parameters: KratosMultiphysics.Parameters):
        defaults = LineOutputProcess.GetDefaultParameters()
        parameters.ValidateAndAssignDefaults(defaults)

        line_output_settings =  parameters["line_output_settings"]
        line_output_settings.ValidateAndAssignDefaults(defaults["line_output_settings"])

        start_point = line_output_settings["start_point"].GetVector()
        end_point = line_output_settings["end_point"].GetVector()
        number_of_points = line_output_settings["number_of_points"].GetInt()

        if number_of_points < 2:
            raise RuntimeError("LineOutputProcess requires at least 2 output points")

        step_vector = (end_point - start_point) / (number_of_points - 1)

        point_set_output_parameters = parameters.Clone()
        point_set_output_parameters.RemoveValue("line_output_settings")
        point_set_output_parameters.AddValue("point_output_settings", KratosMultiphysics.Parameters("""{
            "prefix"              : "/VertexData/Coordinates",
            "time_format"         : "0.4f",
            "positions"           : [],
            "search_configuration": "initial",
            "search_tolerance"    : 1e-6
        }"""))

        point_output_settings = point_set_output_parameters["point_output_settings"]
        point_output_settings["prefix"].SetString(line_output_settings["prefix"].GetString())
        point_output_settings["time_format"].SetString(line_output_settings["time_format"].GetString())
        point_output_settings["search_configuration"].SetString(line_output_settings["search_configuration"].GetString())
        point_output_settings["search_tolerance"].SetInt(line_output_settings["search_tolerance"].GetInt())
        point_output_settings.AddValue("custom_attributes", line_output_settings["custom_attributes"])

        for i_point in range(number_of_points):
            point_output_settings["positions"].Append(start_point + i_point * step_vector)

        return point_set_output_parameters
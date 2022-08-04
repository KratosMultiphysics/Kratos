from os import stat
import KratosMultiphysics
import KratosMultiphysics.HDF5Application as HDF5Application
from KratosMultiphysics.HDF5Application.point_set_output_process import PointSetOutputProcess


def Factory(parameters: KratosMultiphysics.Parameters,
            model: KratosMultiphysics.Model):
    if not isinstance(parameters, KratosMultiphysics.Parameters):
        raise Exception("expecting input parameters of type KratosMultiphysics.Parameters, but got {}".format(type(parameters)))
    return LineOutputProcess(model, parameters["Parameters"])


class LineOutputProcess(KratosMultiphysics.OutputProcess):
    """Convenience process wrapping PointSetOutputProcess with points along a line."""

    def __init__(self, model: KratosMultiphysics.Model, parameters: KratosMultiphysics.Parameters):
        KratosMultiphysics.OutputProcess.__init__(self)
        self.wrapped_process = PointSetOutputProcess(model, self.__ConvertToPointSetOutputParameters(parameters))


    def IsOutputStep(self):
        return self.wrapped_process.IsOutputStep()


    def ExecuteInitialize(self):
        self.wrapped_process.ExecuteInitialize()


    def PrintOutput(self):
        self.wrapped_process.PrintOutput()


    @staticmethod
    def GetDefaultParameters():
        return KratosMultiphysics.Parameters("""{
            "model_part_name"       : "",
            "interval"              : [0.0, "End"],
            "output_frequency"      : 1,
            "start_point"           : [0.0, 0.0, 0.0],
            "end_point"             : [0.0, 0.0, 0.0],
            "number_of_points"      : 2,
            "output_variables"      : [],
            "historical_value"      : true,
            "search_configuration"  : "initial",
            "search_tolerance"      : 1e-6,
            "coordinates_prefix"    : "/<model_part_name>_point_set_output",
            "variables_prefix"      : "/<model_part_name>_point_set_output/step_<step>",
            "file_parameters"       : {
                "file_name"         : "",
                "file_access_mode"  : "read_write",
                "echo_level"        : 0
            }
        }""")


    @staticmethod
    def __ConvertToPointSetOutputParameters(parameters: KratosMultiphysics.Parameters):
        parameters.ValidateAndAssignDefaults(LineOutputProcess.GetDefaultParameters())
        start_point = parameters["start_point"].GetVector()
        end_point = parameters["end_point"].GetVector()
        number_of_points = parameters["number_of_points"].GetInt()

        if number_of_points < 2:
            raise RuntimeError("LineOutputProcess requires at least 2 output points")

        step_vector = (end_point - start_point) / (number_of_points - 1)

        point_set_output_parameters = parameters.Clone()
        point_set_output_parameters.RemoveValue("start_point")
        point_set_output_parameters.RemoveValue("end_point")
        point_set_output_parameters.RemoveValue("number_of_points")
        point_set_output_parameters.AddValue("positions", KratosMultiphysics.Parameters("[]"))

        for i_point in range(number_of_points):
            point_set_output_parameters["positions"].Append(start_point + i_point * step_vector)

        return point_set_output_parameters
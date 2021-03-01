"""
This module contains an interface to the available response functions
"""
import time as timer
import KratosMultiphysics as Kratos
from KratosMultiphysics.response_functions.response_function_interface import ResponseFunctionInterface

class AverageLocationDeviationResponseFunction(ResponseFunctionInterface):
    def __init__(self, identifier, response_settings, model):
        self.identifier = identifier
        self.response_settings = response_settings

        default_parameters = Kratos.Parameters("""
        {
            "response_type"            : "average_location_deviation_response",
            "model_part_name"          : "PLEASE_SPECIFY_MODEL_PART_NAME"
        }
        """)

        self.response_settings.ValidateAndAssignDefaults(default_parameters)
        self.model = model

    def Initialize(self):
        self.model_part = self.model[self.response_settings["model_part_name"].GetString()]

        self.model_part_center = Kratos.Array3(0.0)
        self.number_of_nodes = self.model_part.GetCommunicator().GlobalNumberOfNodes()

        for node in self.model_part.Nodes:
            self.model_part_center[0] += node.X
            self.model_part_center[1] += node.Y
            self.model_part_center[2] += node.Z

        self.model_part_center = self.model_part.GetCommunicator().GetDataCommunicator().SumAll(self.model_part_center)
        self.model_part_center /= self.number_of_nodes

    def CalculateValue(self):
        start_time = timer.time()

        average_location = Kratos.Array3(0.0)
        for node in self.model_part.Nodes:
            average_location[0] += node.X
            average_location[1] += node.Y
            average_location[2] += node.Z

        average_location = self.model_part.GetCommunicator().GetDataCommunicator().SumAll(average_location)

        self.value_array = (average_location / self.number_of_nodes  - self.model_part_center)
        self.value = self.value_array[0] ** 2 + self.value_array[1] ** 2 + self.value_array[2] ** 2

        Kratos.Logger.PrintInfo(self._GetLabel(), "Time needed for calculating the response value = ",round(timer.time() - start_time,2),"s")

    def CalculateGradient(self):
        pass

    def GetValue(self):
        return self.value

    def GetNodalGradient(self, variable):
        if variable != Kratos.SHAPE_SENSITIVITY:
            raise RuntimeError("GetNodalGradient: No gradient for {}!".format(variable.Name))

        gradient = {}
        for node in self.model_part.Nodes:
            gradient[node.Id] = 2.0 * self.value_array / self.number_of_nodes

        return gradient

    @staticmethod
    def _GetLabel():
        return "AverageLocationDeviationResponseFunction"



import KratosMultiphysics as KM
from KratosMultiphysics.NeuralNetworkApplication.neural_network_process import NeuralNetworkProcess
from KratosMultiphysics.NeuralNetworkApplication.data_loading_utilities import ImportDataFromFile
import numpy as np


def Factory(settings):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a parameters object, encapsulating a json string")
    return InvertTransformationsProcess(settings["parameters"])

class InvertTransformationsProcess(NeuralNetworkProcess):

    def __init__(self,parameters):

        super().__init__()


        default_settings = KM.Parameters("""{
            "echo"                : 0,
            "input_file"          : "",
            "predictions_file"    : "",
            "input_file_name"     : "",
            "output_file_name"    : "",
            "save_format"         : "ascii"            
        }""")
        parameters.ValidateAndAssignDefaults(default_settings)

        self.echo = parameters["echo"].GetInt()
        self.input_file = parameters["input_file"].GetString()
        self.predictions_file = parameters["predictions_file"].GetString()
        self.input_file_name = parameters["input_file_name"].GetString()
        self.output_file_name = parameters["output_file_name"].GetString()
        self.save_format = parameters["save_format"].GetString()
        if  not (self.save_format == "ascii") and  not (self.save_format == "npy"):
            raise Exception("Saving format not supported. Supported formats are ascii and npy.")

    def TransformPredictions(self, processes):
        # Load the data
        self.data_in = ImportDataFromFile(self.input_file, "InputData")
        self.data_out = ImportDataFromFile(self.predictions_file, "OutputData")

        # Invert the transformations
        for process in reversed(processes):
            try:
                if not process.load_from_log:
                    inversion_settings = process.Invert(self.data_in, self.data_out)
                    if not inversion_settings is None:
                        self.data_in = inversion_settings[0]
                        self.data_out = inversion_settings[1]
            except AttributeError:
                continue

        # Save the data
        if self.save_format == "ascii":
            extension  = ".csv"
            self.input_file_name = self.input_file_name + extension
            self.output_file_name = self.output_file_name + extension
            np.savetxt(self.input_file_name, self.data_in)
            np.savetxt(self.output_file_name,self.data_out)

        elif self.save_format == "npy":
            extension = ".npy"
            self.input_file_name = self.input_file_name + extension
            self.output_file_name = self.output_file_name + extension
            np.save(self.input_file_name, self.data_in)
            np.save(self.output_file_name,self.data_out)
        else:
            raise Exception("Output save format not supported. Supported formats are ascii and npy.")



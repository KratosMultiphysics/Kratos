import KratosMultiphysics as KM
from KratosMultiphysics.NeuralNetworkApplication.neural_network_process import NeuralNetworkProcess
from KratosMultiphysics.NeuralNetworkApplication.data_loading_utilities import ImportDataFromFile
import numpy as np
import matplotlib.pyplot as plt
import keras.models


def Factory(settings):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a parameters object, encapsulating a json string")
    return PredictFromModelProcess(settings["parameters"])

class PredictFromModelProcess(NeuralNetworkProcess):

    def __init__(self,parameters):

        super().__init__()

        default_settings = KM.Parameters("""{
            "model"               : "",
            "echo"                : 0,
            "write_output"        : false,
            "predictions_file"    : "",
            "save_format"         : "ascii"            
        }""")
        parameters.ValidateAndAssignDefaults(default_settings)

        self.model_file = parameters["model"].GetString()
        self.echo = parameters["echo"].GetInt()
        self.write_output = parameters["write_output"].GetBool()
        self.predictions_file = parameters["predictions_file"].GetString()
        self.save_format = parameters["save_format"].GetString()
        if  not (self.save_format == "ascii") and  not (self.save_format == "npy"):
            raise Exception("Saving format not supported. Supported formats are ascii and npy.")

    def Predict(self, model, data_in):

        predictions = model.predict(data_in.ExportAsArray())

        if self.write_output:

            if self.save_format == "ascii":
                self.predictions_file = self.predictions_file + '.csv'
                np.savetxt(self.predictions_file, predictions)
            
            if self.save_format == "npy":
                self.predictions_file = self.predictions_file + '.npy'
                np.save(self.predictions_file, predictions)

        return predictions
                



import KratosMultiphysics as KM
from KratosMultiphysics.NeuralNetworkApplication.neural_network_process import NeuralNetworkProcess
from KratosMultiphysics.NeuralNetworkApplication.data_loading_utilities import ImportDataFromFile
import numpy as np
import matplotlib.pyplot as plt
import keras.models


def Factory(settings):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a parameters object, encapsulating a json string")
    return TestingProcess(settings["parameters"])

class TestingProcess(NeuralNetworkProcess):

    def __init__(self,parameters):

        super().__init__()


        default_settings = KM.Parameters("""{
            "model"               : "",
            "input_file"          : "",
            "target_file"         : "",
            "echo"                : 0,
            "predictions_file"    : "",
            "save_format"         : "ascii"            
        }""")
        parameters.ValidateAndAssignDefaults(default_settings)

        self.model = parameters["model"].GetString()
        self.input_file = parameters["input_file"].GetString()
        self.target_file = parameters["target_file"].GetString()
        self.echo = parameters["echo"].GetInt()
        self.predictions_file = parameters["predictions_file"].GetString()
        self.save_format = parameters["save_format"].GetString()
        if  not (self.save_format == "ascii") and  not (self.save_format == "npy"):
            raise Exception("Saving format not supported. Supported formats are ascii and npy.")

    def ExecuteFinalize(self):

        print("Testing process started... \n")
        input = ImportDataFromFile(self.input_file, "InputData").ExportAsArray()
        target = ImportDataFromFile(self.target_file, "OutputData").ExportAsArray()
        # The compiling options of the model do not matter, as it will not be trained
        model = keras.models.load_model(self.model)
        model.compile()
        predictions = model.predict(input)

        if self.echo > 0:
            timestep = 0
            for entry in predictions:
                print('Predicted: ' + str(entry))
                print('Ground truth: ' + str(target[timestep]))
                timestep = timestep + 1
        
        if self.save_format == "ascii":
            self.predictions_file = self.predictions_file + '.csv'
            np.savetxt(self.predictions_file, predictions)
        

        if self.save_format == "npy":
            self.predictions_file = self.predictions_file + '.npy'
            np.save(self.predictions_file, predictions)
            



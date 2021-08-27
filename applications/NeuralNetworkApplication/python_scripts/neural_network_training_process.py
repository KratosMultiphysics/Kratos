from importlib import import_module
import KratosMultiphysics as KM
from KratosMultiphysics.NeuralNetworkApplication.neural_network_process import NeuralNetworkProcess
import KratosMultiphysics.NeuralNetworkApplication.data_loading_utilities as DataLoadingUtilities
import numpy as np
import h5py

def Factory(settings):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a parameters object, encapsulating a json string")
    return NeuralNetworkTrainingProcess(settings["parameters"])


class NeuralNetworkTrainingProcess(NeuralNetworkProcess):
    """
    This class is the parent for the classes that deal with the training of the neural network.

    It is similar to the normal Kratos Process but it does not receive a model, only settings.
    """
    def __init__(self, parameters):
        super().__init__()
        self.input_file = parameters["data_input"].GetString()
        self.output_file = parameters["data_output"].GetString()
        self.validation = False
        self.validation_split = 0.0
        if parameters.Has("val_input"):
            self.validation = True
            self.val_input_file = parameters["val_input"].GetString()
        if parameters.Has("val_output"):
            self.validation = True
            self.val_output_file = parameters["val_output"].GetString()
        if parameters.Has("validation_split"):
            self.validation_split = parameters["validation_split"].GetDouble()
        if parameters.Has("recurrent"):
            self.recurrent = parameters["recurrent"].GetBool()
        else:
            self.recurrent = False

    def ExecuteBeforeTraining(self):
        # Input
        
        self.test_input = DataLoadingUtilities.ImportDataFromFile(self.input_file,"InputData", lookback=self.recurrent)

        # Output
        
        self.test_output = DataLoadingUtilities.ImportDataFromFile(self.output_file,"OutputData")
    

        # Validation files
        
        if self.validation:
            # Input
        
            self.val_test_input = DataLoadingUtilities.ImportDataFromFile(self.val_input_file,"InputData", lookback=self.recurrent)

            # Output
            
            self.val_test_output = DataLoadingUtilities.ImportDataFromFile(self.val_output_file,"OutputData")

        



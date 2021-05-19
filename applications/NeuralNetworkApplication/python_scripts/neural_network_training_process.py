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

    def ExecuteBeforeTraining(self):
        # Input
        
        # Data loading for h5
        if self.input_file.endswith('.h5'):
            self.test_input = DataLoadingUtilities.ImportH5(self.input_file, "InputData")
        # Data loading for dat
        elif self.input_file.endswith('.dat'):
            self.test_input = DataLoadingUtilities.ImportAscii(self.input_file)
        # Exception for non-supported formats
        else:
            raise Exception("Input data format not supported. Supported formats are .dat and .h5")

        # Output
        
        # Data loading for h5
        if self.output_file.endswith('.h5'):
            self.test_output = DataLoadingUtilities.ImportH5(self.output_file,"OutputData")
  
        # Data loading for dat
        elif self.output_file.endswith('.dat'):
            self.test_output = DataLoadingUtilities.ImportAscii(self.output_file)
        # Exception for non-supported formats
        else:
            raise Exception("Output data format not supported. Supported formats are .dat and .h5")

        # Validation files
        
        if self.validation:
            # Data loading for h5
            if self.val_input_file.endswith('.h5'):
                self.val_input = DataLoadingUtilities.ImportH5(self.val_input_file, "InputData")
            # Data loading for dat
            elif self.val_input_file.endswith('.dat'):
                self.val_input = DataLoadingUtilities.ImportAscii(self.val_input_file)
            # Exception for non-supported formats
            else:
                raise Exception("Input validation data format not supported. Supported formats are .dat and .h5")
            # Data loading for h5
            if self.val_output_file.endswith('.h5'):
                self.val_output = DataLoadingUtilities.ImportH5(self.val_output_file, "OutputData")
            # Data loading for dat
            elif self.val_output_file.endswith('.dat'):
                self.val_output = DataLoadingUtilities.ImportAscii(self.val_output_file)
            # Exception for non-supported formats
            else:
                raise Exception("Output validation data format not supported. Supported formats are .dat and .h5")

        



import numpy as np
import warnings
import KratosMultiphysics as KM
from  KratosMultiphysics.NeuralNetworkApplication.neural_network_process import NeuralNetworkProcess


def Factory(settings):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return PreprocessingProcess(settings["parameters"])

class PreprocessingProcess(NeuralNetworkProcess):

    def __init__(self, settings):
        super().__init__()
        """ The default constructor of the class

        Keyword arguments:
        self -- It signifies an instance of a class.
        model -- the container of the different model parts.
        settings -- Kratos parameters containing process settings.
        """
        if settings.Has("load_from_log"): 
            self.load_from_log = settings["load_from_log"].GetBool()
        else:
            self.load_from_log = False
        if settings.Has("input_log_name"):
            self.input_log_name = settings["input_log_name"].GetString()
        else:
            if self.load_from_log:
                warnings.warn("No input_log_file is specified. Data cannot be loaded.")
            else:
                warnings.warn("No input_log_file is specified. The transformation will not be logged.")

        if settings.Has("output_log_name"):
            self.output_log_name = settings["output_log_name"].GetString()
        else:
            if self.load_from_log:
                warnings.warn("No output_log_file is specified. Data cannot be loaded.")
            else:
                warnings.warn("No output_log_file is specified. The transformation will not be logged.")

            

    def Preprocess(self, data_in, data_out):
        """Preprocessing of the data."""
        return [data_in, data_out]
    
    def Invert(self, data_in, data_out):
        """Inverts de transformation."""
        return [data_in, data_out]

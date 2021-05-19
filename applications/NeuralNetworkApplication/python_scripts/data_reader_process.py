import KratosMultiphysics as KM
from KratosMultiphysics.NeuralNetworkApplication.preprocessing_process import PreprocessingProcess
import KratosMultiphysics.NeuralNetworkApplication.data_loading_utilities as DataLoadingUtilities
import numpy as np
import json
from importlib import import_module

def Factory(settings):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return DataReaderProcess(settings["parameters"])

class DataReaderProcess(PreprocessingProcess):

    def __init__(self, parameters):
        super().__init__(parameters)
        # Input
        self.input_file = parameters["data_input"].GetString()
        # Data loading for h5
        if self.input_file.endswith('.h5'):
            self.input = DataLoadingUtilities.ImportH5(self.input_file, "InputData")
        # Data loading for dat
        elif self.input_file.endswith('.dat'):
            self.input = DataLoadingUtilities.ImportAscii(self.input_file)
        # Exception for non-supported formats
        else:
            raise Exception("Input data format not supported. Supported formats are .dat and .h5")

        # Output
        self.output_file = parameters["data_output"].GetString()
        # Data loading for h5
        if self.output_file.endswith('.h5'):
            self.output = DataLoadingUtilities.ImportH5(self.output_file,"OutputData")

        # Data loading for dat
        elif self.output_file.endswith('.dat'):
            self.output = DataLoadingUtilities.ImportAscii(self.output_file)
        # Exception for non-supported formats
        else:
            raise Exception("Output data format not supported. Supported formats are .dat and .h5")

    def Preprocess(self, data_in, data_out):
        if not self.load_from_log:
            dic_in = {"data_input" : self.input_file}
            with open(self.input_log_name,'w') as f:
                json.dump(dic_in,f)

            dic_out = {"data_output" : self.input_file}
            with open(self.output_log_name,'w') as f:
                json.dump(dic_out,f)

        # Return the data
        return [self.input, self.output]
import KratosMultiphysics as KM
from KratosMultiphysics.NeuralNetworkApplication.preprocessing_process import PreprocessingProcess
import KratosMultiphysics.NeuralNetworkApplication.data_loading_utilities as DataLoadingUtilities
import json


def Factory(settings):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return DataReaderProcess(settings["parameters"])

class DataReaderProcess(PreprocessingProcess):

    def __init__(self, parameters):
        super().__init__(parameters)
        # Input
        self.input_file = parameters["data_input"].GetString()
        self.input = DataLoadingUtilities.ImportDataFromFile(self.input_file, "InputData")
        # Output
        self.output_file = parameters["data_output"].GetString()
        self.output = DataLoadingUtilities.ImportDataFromFile(self.output_file, "OutputData")

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
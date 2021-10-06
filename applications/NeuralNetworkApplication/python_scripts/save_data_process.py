import KratosMultiphysics as KM
from KratosMultiphysics.NeuralNetworkApplication.preprocessing_process import PreprocessingProcess
import KratosMultiphysics.NeuralNetworkApplication.data_loading_utilities as DataLoadingUtilities
import numpy as np
import pickle

def Factory(settings):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return SaveDataProcess(settings["parameters"])

class SaveDataProcess(PreprocessingProcess):

    def __init__(self, parameters):
        """ The default constructor of the class

        Keyword arguments:
        self -- It signifies an instance of a class.
        model -- the container of the different model parts.
        settings -- Kratos parameters containing process settings.
        """
        super().__init__(parameters)
        self.input_file_name = parameters["input_file_name"].GetString()
        self.output_file_name = parameters["output_file_name"].GetString()
        self.format = parameters["format"].GetString()

    def Preprocess(self, data_structure_in, data_structure_out):

        data_in = data_structure_in.ExportAsArray()
        data_out = data_structure_out.ExportAsArray()

        # Save the data
        if self.format == "ascii":
            extension  = ".dat"
            self.input_file_name = self.input_file_name + extension
            self.output_file_name = self.output_file_name + extension
            np.savetxt(self.input_file_name, data_in)
            np.savetxt(self.output_file_name,data_out)

        elif self.format == "npy":
            extension = ".npy"
            self.input_file_name = self.input_file_name + extension
            self.output_file_name = self.output_file_name + extension
            np.save(self.input_file_name, data_in)
            np.save(self.output_file_name,data_out)

        elif self.format == "pkl":
            # This is the prefered method for complex data structures.
            extension = ".pkl"
            self.input_file_name = self.input_file_name + extension
            self.output_file_name = self.output_file_name + extension
            with open(self.input_file_name, 'wb') as in_file:
                pickle.dump(data_structure_in, in_file, pickle.HIGHEST_PROTOCOL)
            with open(self.output_file_name, 'wb') as out_file:
                pickle.dump(data_structure_out, out_file, pickle.HIGHEST_PROTOCOL)
        else:
            raise Exception("Output save format not supported. Supported formats are ascii, pkl and npy.")

        if not self.load_from_log:
            input_log = DataLoadingUtilities.ImportDictionaryFromText(self.input_log_name)
            output_log = DataLoadingUtilities.ImportDictionaryFromText(self.output_log_name)

            input_log.update({ "saved_file_name" : self.input_file_name })
            output_log.update({ "saved_file_name" : self.output_file_name })

            DataLoadingUtilities.UpdateDictionaryJson(self.input_log_name, input_log)
            DataLoadingUtilities.UpdateDictionaryJson(self.output_log_name, output_log)
            
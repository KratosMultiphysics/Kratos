import numpy as np
import KratosMultiphysics as KM
from  KratosMultiphysics.NeuralNetworkApplication.preprocessing_process import PreprocessingProcess
from KratosMultiphysics.NeuralNetworkApplication.data_loading_utilities import ImportDictionaryFromText, UpdateDictionaryJson


def Factory(settings):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return CenteringProcess(settings["parameters"])

class CenteringProcess(PreprocessingProcess):

    def __init__(self, settings):
        super().__init__(settings)
        """ The default constructor of the class

        Keyword arguments:
        self -- It signifies an instance of a class.
        model -- the container of the different model parts.
        settings -- Kratos parameters containing process settings.
        """
       
        if self.load_from_log:
            self.center = "file"
        elif settings.Has("center"):
            self.center = settings["center"].GetString()
        else:
            self.center = "mean"
        self.objective = settings["objective"].GetString()
       

    def Preprocess(self, data_in, data_out):
        try:
            input_log = ImportDictionaryFromText(self.input_log_name)
            output_log = ImportDictionaryFromText(self.output_log_name)
        except AttributeError:
            print("No logging.")
            input_log = {}
            output_log = {}

        # Centering from the mean
        if self.center == "mean":
            if self.objective == "input":
                mean_in = np.mean(data_in, axis = 0)
                data_in = data_in - mean_in
                input_log.update({"centering" : mean_in.tolist()})
            if self.objective == "output":
                mean_out = np.mean(data_out, axis = 0)
                data_out = data_out - mean_out
                output_log.update({"centering" : mean_out.tolist()})
            if self.objective == "all":
                mean_in = np.mean(data_in, axis = 0)
                data_in = data_in - mean_in
                input_log.update({"centering" : mean_in.tolist()})
                mean_out = np.mean(data_out, axis = 0)
                data_out = data_out - mean_out
                output_log.update({"centering" : mean_out.tolist()})

        # Centering from the minimum
        if self.center == "min":
            if self.objective == "input":
                min = data_in.min(axis = 0)
                data_in = data_in - min
                input_log.update({"centering" : min.tolist()})
            if self.objective == "output":
                min = data_out.min(axis = 0)
                data_out = data_out - min
                output_log.update({"centering" : min.tolist()})
            if self.objective == "all":
                min = data_in.min(axis = 0)
                data_in = data_in - min
                input_log.update({"centering" : min.tolist()})
                min = data_out.min(axis = 0)
                data_out = data_out - min
                output_log.update({"centering" : min.tolist()})

        # Centering from file log
        if self.center == "file":
            if self.objective == "input":
                data_in = data_in - input_log.get("centering")
            if self.objective == "output":
                data_out = data_out - output_log.get("centering")
            if self.objective == "all":
                data_in = data_in - input_log.get("centering")
                data_out = data_out - output_log.get("centering")

        # Updating the file log
        if not self.load_from_log:
            try:
                UpdateDictionaryJson(self.input_log_name, input_log)
                UpdateDictionaryJson(self.output_log_name, output_log)
            except AttributeError:
                pass

        return [data_in, data_out]

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
        try:
            self.log_denominator = settings["log_denominator"].GetString()
        except AttributeError:
            self.log_denominator = "centering"
       

    def Preprocess(self, data_structure_in, data_structure_out):

        if len(data_structure_in)>0:
            data_in = data_structure_in.ExportAsArray()
        if len(data_structure_out)>0:
            data_out = data_structure_out.ExportAsArray()

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
                input_log.update({self.log_denominator : mean_in.tolist()})
            if self.objective == "output":
                mean_out = np.mean(data_out, axis = 0)
                data_out = data_out - mean_out
                output_log.update({self.log_denominator : mean_out.tolist()})
            if self.objective == "all":
                mean_in = np.mean(data_in, axis = 0)
                data_in = data_in - mean_in
                input_log.update({self.log_denominator : mean_in.tolist()})
                mean_out = np.mean(data_out, axis = 0)
                data_out = data_out - mean_out
                output_log.update({self.log_denominator : mean_out.tolist()})

        # Centering from the minimum
        if self.center == "min":
            if self.objective == "input":
                min = data_in.min(axis = 0)
                data_in = data_in - min
                input_log.update({self.log_denominator : min.tolist()})
            if self.objective == "output":
                min = data_out.min(axis = 0)
                data_out = data_out - min
                output_log.update({self.log_denominator : min.tolist()})
            if self.objective == "all":
                min = data_in.min(axis = 0)
                data_in = data_in - min
                input_log.update({self.log_denominator : min.tolist()})
                min = data_out.min(axis = 0)
                data_out = data_out - min
                output_log.update({self.log_denominator : min.tolist()})

         # Centering for avoiding the 0 in soft_maxmin normalization
        if self.center == "soft_minmax":
            if self.objective == "input":
                gap = (1.0 - data_in.max(axis = 0))/2.0
                data_in = data_in + gap
                input_log.update({self.log_denominator : (-gap).tolist()})
            if self.objective == "output":
                gap = (1.0 - data_out.max(axis = 0))/2.0
                data_out = data_out + gap
                output_log.update({self.log_denominator : (-gap).tolist()})
            if self.objective == "all":
                gap = (1.0 - data_in.max(axis = 0))/2.0
                data_in = data_in + gap
                input_log.update({self.log_denominator : (-gap).tolist()})
                gap = (1.0 - data_out.max(axis = 0))/2.0
                data_out = data_out + gap
                output_log.update({self.log_denominator: (-gap).tolist()})

        # Centering from file log
        if self.center == "file":
            if self.objective == "input":
                data_in = data_in - input_log.get(self.log_denominator)
            if self.objective == "output":
                data_out = data_out - output_log.get(self.log_denominator)
            if self.objective == "all":
                data_in = data_in - input_log.get(self.log_denominator)
                data_out = data_out - output_log.get(self.log_denominator)

        # Updating the file log
        if not self.load_from_log:
            try:
                UpdateDictionaryJson(self.input_log_name, input_log)
                UpdateDictionaryJson(self.output_log_name, output_log)
            except AttributeError:
                pass
        
        # Updating the data structure
        if len(data_structure_in)>0:
            data_structure_in.UpdateData(data_in)
        if len(data_structure_out)>0:
            data_structure_out.UpdateData(data_out)

        return [data_structure_in, data_structure_out]

    def Invert(self,data_structure_in,data_structure_out):
        
        data_in = data_structure_in.ExportAsArray()
        data_out = data_structure_out.ExportAsArray()

        input_log = ImportDictionaryFromText(self.input_log_name)
        output_log = ImportDictionaryFromText(self.output_log_name)

        if self.objective == "input" or self.objective == "predict_input":
            data_in = data_in + input_log.get(self.log_denominator)
        if self.objective == "output" or self.objective == "predict_output":
            data_out = data_out + output_log.get(self.log_denominator)
        if self.objective == "all" or self.objective == "predict_all":
            data_in = data_in + input_log.get(self.log_denominator)
            data_out = data_out + output_log.get(self.log_denominator)

        data_structure_in.UpdateData(data_in)
        data_structure_out.UpdateData(data_out)

        return[data_structure_in,data_structure_out]

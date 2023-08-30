import numpy as np
import KratosMultiphysics as KM
from  KratosMultiphysics.NeuralNetworkApplication.preprocessing_process import PreprocessingProcess
from KratosMultiphysics.NeuralNetworkApplication.data_loading_utilities import ImportDictionaryFromText, UpdateDictionaryJson

def Factory(settings):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ScalingProcess(settings["parameters"])

class ScalingProcess(PreprocessingProcess):

    def __init__(self, settings):
        super().__init__(settings)
        """ The default constructor of the class

        Keyword arguments:
        self -- It signifies an instance of a class.
        model -- the container of the different model parts.
        settings -- Kratos parameters containing process settings.
        """

        if self.load_from_log:
            self.scale = "file"
        elif settings.Has("scale"):
            self.scale = settings["scale"].GetString()
        else:
            self.scale = "std"
        self.objective = settings["objective"].GetString()
        try:
            self.log_denominator = settings["log_denominator"].GetString()
        except RuntimeError:
            self.log_denominator = "scaling"
       

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
        # Scaling from the standard deviation
        if self.scale == "std":
            if self.objective == "input":
                std_in = np.std(data_in, axis = 0)
                std_in = self.__check_for_zero(std_in)
                data_in = data_in / std_in
                input_log.update({self.log_denominator : std_in.tolist()})
            if self.objective == "output":
                std_out = np.std(data_out, axis = 0)
                std_out = self.__check_for_zero(std_out)
                data_out = data_out / std_out
                output_log.update({self.log_denominator : std_out.tolist()})
            if self.objective == "all":
                std_in = np.std(data_in, axis = 0)
                std_in = self.__check_for_zero(std_in)
                data_in = data_in / std_in
                input_log.update({self.log_denominator : std_in.tolist()})
                std_out = np.std(data_out, axis = 0)
                std_out = self.__check_for_zero(std_out)
                data_out = data_out / std_out
                output_log.update({self.log_denominator : std_out.tolist()})

        # Scaling from the min-max
        if self.scale == "minmax":
            if self.objective == "input":
                range = data_in.max(axis = 0) - data_in.min(axis = 0)
                range = self.__check_for_zero(range)
                data_in = data_in/ range
                input_log.update({self.log_denominator : range.tolist()})
            if self.objective == "output":
                range = data_out.max(axis = 0) - data_out.min(axis = 0)
                range = self.__check_for_zero(range)
                data_out = data_out/ range
                output_log.update({self.log_denominator : range.tolist()})
            if self.objective == "all":
                range = data_in.max(axis = 0) - data_in.min(axis = 0)
                range = self.__check_for_zero(range)
                data_in = data_in/ range
                input_log.update({self.log_denominator : range.tolist()})
                range = data_out.max(axis = 0) - data_out.min(axis = 0)
                range = self.__check_for_zero(range)
                data_out = data_out/ range
                output_log.update({self.log_denominator : range.tolist()})

        # Scaling from the soft min-max
        if self.scale == "soft_minmax":
            if self.objective == "input":
                range = 1.25*(data_in.max(axis = 0) - data_in.min(axis = 0))
                range = self.__check_for_zero(range)
                data_in = data_in/ range
                input_log.update({self.log_denominator : range.tolist()})
            if self.objective == "output":
                range = 1.25*(data_out.max(axis = 0) - data_out.min(axis = 0))
                range = self.__check_for_zero(range)
                data_out = data_out/ range
                output_log.update({self.log_denominator : range.tolist()})
            if self.objective == "all":
                range = 1.25*(data_in.max(axis = 0) - data_in.min(axis = 0))
                range = self.__check_for_zero(range)
                data_in = data_in/ range
                input_log.update({self.log_denominator : range.tolist()})
                range = 1.25*(data_out.max(axis = 0) - data_out.min(axis = 0))
                range = self.__check_for_zero(range)
                data_out = data_out/ range
                output_log.update({self.log_denominator : range.tolist()})

        # Scaling from file log
        if self.scale == "file":
            if self.objective == "input":
                data_in = data_in / input_log.get(self.log_denominator)
            if self.objective == "output":
                data_out = data_out / output_log.get(self.log_denominator)
            if self.objective == "all":
                data_in = data_in / input_log.get(self.log_denominator)
                data_out = data_out / output_log.get(self.log_denominator)
                
        # Updating the file log
        if not self.load_from_log:
            try:
                UpdateDictionaryJson(self.input_log_name, input_log)
                UpdateDictionaryJson(self.output_log_name, output_log)
            except AttributeError:
                pass
        if len(data_structure_in)>0:
            data_structure_in.UpdateData(data_in)
        if len(data_structure_out)>0:
            data_structure_out.UpdateData(data_out)

        return [data_structure_in, data_structure_out]

    def Invert(self,data_structure_in,data_structure_out):
        
        data_in = data_structure_in.ExportDataOnly()
        data_out = data_structure_out.ExportDataOnly()

        input_log = ImportDictionaryFromText(self.input_log_name)
        output_log = ImportDictionaryFromText(self.output_log_name)

        if self.objective == "input" or self.objective == "predict_input":
            data_in = data_in * input_log.get(self.log_denominator)
        if self.objective == "output" or self.objective == "predict_output":
            data_out = data_out * output_log.get(self.log_denominator)
        if self.objective == "all" or self.objective == "predict_all":
            data_in = data_in * input_log.get(self.log_denominator)
            data_out = data_out * output_log.get(self.log_denominator)


        data_structure_in.UpdateData(data_in)
        data_structure_out.UpdateData(data_out)

        return[data_structure_in,data_structure_out]

    def __check_for_zero(self, scale_vector):

        try:
            scale_vector[np.nonzero(abs(scale_vector)<1e-8)] = 1.0
        except TypeError:
            pass
        
        return scale_vector



import numpy as np
import KratosMultiphysics as KM
from  KratosMultiphysics.NeuralNetworkApplication.preprocessing_process import PreprocessingProcess
from KratosMultiphysics.NeuralNetworkApplication.data_loading_utilities import ImportDictionaryFromText, UpdateDictionaryJson
from KratosMultiphysics.NeuralNetworkApplication.input_dataclasses import ListNeuralNetworkData, ListDataWithLookback

def Factory(settings):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return FlattenVectorProcess(settings["parameters"])

class FlattenVectorProcess(PreprocessingProcess):

    def __init__(self, settings):
        super().__init__(settings)
        """ The default constructor of the class

        Keyword arguments:
        self -- It signifies an instance of a class.
        model -- the container of the different model parts.
        settings -- Kratos parameters containing process settings.
        """
        self.objective = settings["objective"].GetString()
        try:
            self.log_denominator = settings["log_denominator"].GetString()
        except RuntimeError:
            self.log_denominator = "flatten"
       

    def Preprocess(self, data_structure_in, data_structure_out):
        
        data_in = data_structure_in.ExportAsArray()
        data_out = data_structure_out.ExportAsArray()

        # Flatten for input
        if self.objective == "input":
            input_log = ImportDictionaryFromText(self.input_log_name)
            try:
                new_data_in = data_in.reshape((data_in.shape[0],data_in.shape[1]*data_in.shape[2]))
            except IndexError:
                new_data_in = data_in.reshape((data_in.shape[0]*data_in.shape[1]))
            new_data_out = data_out
            input_log.update({self.log_denominator : list(data_in.shape[1:])})
        # Flatten for predict input
        if self.objective == "predict_input":
            input_log = ImportDictionaryFromText(self.input_log_name)
            try:
                new_data_in = data_in.reshape((data_in.shape[0],data_in.shape[1]*data_in.shape[2]))
            except IndexError:
                new_data_in = data_in.reshape((data_in.shape[0]*data_in.shape[1]))
            new_data_out = data_out
            input_log.update({self.log_denominator : list(data_in.shape)})
        # Flatten for output
        if self.objective == "output":
            output_log = ImportDictionaryFromText(self.output_log_name)
            new_data_in = data_in
            try:
                new_data_out = data_out.reshape((data_out.shape[0],data_out.shape[1]*data_out.shape[2]))
            except IndexError:
                new_data_out = data_out.reshape((data_out.shape[0]*data_out.shape[1]))

            output_log.update({self.log_denominator : list(data_out.shape[1:])})
        if self.objective == "output":
            output_log = ImportDictionaryFromText(self.output_log_name)
            new_data_in = data_in
            try:
                new_data_out = data_out.reshape((data_out.shape[0],data_out.shape[1]*data_out.shape[2]))
            except IndexError:
                new_data_out = data_out.reshape((data_out.shape[0]*data_out.shape[1]))

            output_log.update({self.log_denominator : list(data_out.shape)})

        # Updating the file log
        if not self.load_from_log:
            try:
                UpdateDictionaryJson(self.input_log_name, input_log)
            except AttributeError:
                pass
            try:
                UpdateDictionaryJson(self.output_log_name, output_log)
            except AttributeError:
                pass
        
        data_structure_in.UpdateData(new_data_in)
        data_structure_out.UpdateData(new_data_out)
            
        return [data_structure_in, data_structure_out]

    def Invert(self, data_structure_in, data_structure_out):

        data_in_array = data_structure_in.ExportDataOnly()
        data_out_array = data_structure_out.ExportDataOnly()
        i = 0
        for data_in in data_in_array:
            
            if self.objective == "input":
                input_log = ImportDictionaryFromText(self.input_log_name)
                input_shape = input_log.get(self.log_denominator)
                try:
                    new_data_in = data_in.reshape(input_shape[1],input_shape[2])
                except IndexError:
                    new_data_in = data_in.reshape(data_in.shape[0],input_shape[0])
            if self.objective == "predict_input":
                input_log = ImportDictionaryFromText(self.input_log_name)
                input_shape = input_log.get(self.log_denominator)
                try:
                    new_data_in = data_in.reshape(input_shape[0],input_shape[1])
                except IndexError:
                    new_data_in = data_in.reshape(input_shape[0])
            if self.objective == "output":
                new_data_in = data_in
            if self.objective == "predict_output":
                new_data_in = data_in
            if self.objective == "all":
                input_log = ImportDictionaryFromText(self.input_log_name)
                input_shape = input_log.get(self.log_denominator)
                try:
                    new_data_in = data_in.reshape(data_in.shape[0],input_shape[0],input_shape[1])
                except IndexError:
                    new_data_in = data_in.reshape(data_in.shape[0],input_shape[0])

            if self.objective == "predict_all":
                input_log = ImportDictionaryFromText(self.input_log_name)
                input_shape = input_log.get(self.log_denominator)
                try:
                    new_data_in = data_in.reshape(input_shape[0],input_shape[1])
                except IndexError:
                    new_data_in = data_in.reshape(input_shape[0])
            try:
                data_structure_in[i].UpdateData(new_data_in)
            except AttributeError:
                data_structure_in.UpdateData(new_data_in)

        for data_out in data_out_array:
            
            if self.objective == "input":
                new_data_out = data_out
            if self.objective == "predict_input":
                new_data_out = data_out
            if self.objective == "output":
                output_log = ImportDictionaryFromText(self.output_log_name)
                output_shape = output_log.get(self.log_denominator)
                try:
                    new_data_out = data_out.reshape(output_shape[1],output_shape[2])
                except IndexError:
                    new_data_out = data_out.reshape(data_out.shape[0],output_shape[0])
            if self.objective == "predict_output":
                output_log = ImportDictionaryFromText(self.output_log_name)
                output_shape = output_log.get(self.log_denominator)
                try:
                    new_data_out = data_out.reshape(output_shape[0],output_shape[1])
                except IndexError:
                    new_data_out = data_out.reshape(output_shape[0])
            if self.objective == "all":
                output_log = ImportDictionaryFromText(self.output_log_name)
                output_shape = output_log.get(self.log_denominator)
                try:
                    new_data_out = data_out.reshape(data_out.shape[0],output_shape[0],output_shape[1])
                except IndexError:
                    new_data_out = data_out.reshape(data_out.shape[0],output_shape[0])
            if self.objective == "predict_all":
                output_log = ImportDictionaryFromText(self.output_log_name)
                output_shape = output_log.get(self.log_denominator)
                try:
                    new_data_out = data_out.reshape(output_shape[0],output_shape[1])
                except IndexError:
                    new_data_out = data_out.reshape(output_shape[0])
            try:
                data_structure_out[i].UpdateData(new_data_out)
            except AttributeError:
                data_structure_out.UpdateData(new_data_out)
            i+=1

        return [data_structure_in, data_structure_out]
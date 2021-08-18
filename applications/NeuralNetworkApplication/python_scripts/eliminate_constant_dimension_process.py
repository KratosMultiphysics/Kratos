import KratosMultiphysics as KM
from KratosMultiphysics.NeuralNetworkApplication.preprocessing_process import PreprocessingProcess
from KratosMultiphysics.NeuralNetworkApplication.data_loading_utilities import ImportDictionaryFromText, UpdateDictionaryJson
import numpy as np


def Factory(settings):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return EliminateConstantDimensionProcess(settings["parameters"])

class EliminateConstantDimensionProcess(PreprocessingProcess):

    def __init__(self, settings):
        super().__init__(settings)
        """ The default constructor of the class

        Keyword arguments:
        self -- It signifies an instance of a class.
        settings -- Kratos parameters containing process settings.
        """

        self.objective = settings["objective"].GetString()
        try:
            self.log_denominator = settings["log_denominator"].GetString()
        except RuntimeError:
            self.log_denominator = "eliminate_constant_dimension"
        
    def Preprocess(self, data_in, data_out):

        # Setup for input
        if self.objective == 'input':
            # dimensions = range(len(data_in[0,0]))
            non_constant_dim = np.zeros(data_in.shape[2], dtype=np.bool)
            # dataset_length = range(len(data_in[:,0]))
            # Extract which variables are always 0
            for ix in range(data_in.shape[2]):
                non_constant_dim[ix] = (abs(data_in[...,ix]- data_in[0,0,ix]) > 1e-8).any()
            new_data_in = data_in[...,non_constant_dim]

            try:
                input_log = ImportDictionaryFromText(self.input_log_name)
                input_log.update({self.log_denominator: non_constant_dim.astype(int).tolist()})
                UpdateDictionaryJson(self.input_log_name, input_log)
            except AttributeError:
                print("Input not logged")

            return [new_data_in, data_out]
            
        # Setup for output
        if self.objective == 'output':
            # dimensions = range(len(data_out[0,0]))
            non_constant_dim = np.zeros(data_out.shape[2], dtype=np.bool)
            # dataset_length = range(len(data_out[:,0]))
            # Extract which variables are always 0
            for ix in range(data_out.shape[2]):
                non_constant_dim[ix] = (abs(data_out[...,ix]- data_out[0,0,ix]) > 1e-8).any()
            new_data_out = data_out[...,non_constant_dim]

            try:
                output_log = ImportDictionaryFromText(self.output_log_name)
                output_log.update({self.log_denominator: non_constant_dim.astype(int).tolist()})
                UpdateDictionaryJson(self.output_log_name, output_log)
            except AttributeError:
                print("Output not logged")

            return [data_in, new_data_out]
        else:
            raise Exception("Masking objective not supported. Supported objectives are input and output")

    def Invert(self, data_in, data_out):  

        # This only works if the constant dimensions are 0!!

        if self.objective == "input" or self.objective == "predict_input":
            input_log = ImportDictionaryFromText(self.input_log_name)
            non_constant_dim = list(map(bool,input_log.get(self.log_denominator)))
            new_data_in = np.zeros((data_in.shape[0],data_in.shape[1],len(non_constant_dim)))
            new_data_in[:,:,non_constant_dim]=data_in
            new_data_out = data_out

        if self.objective == "output" or self.objective == "predict_output":
            output_log = ImportDictionaryFromText(self.output_log_name)
            non_constant_dim = list(map(bool,output_log.get(self.log_denominator)))
            new_data_out = np.zeros((data_out.shape[0],data_out.shape[1],len(non_constant_dim)))
            new_data_out[:,:,non_constant_dim]=data_out
            new_data_in = data_in

        if self.objective == "all" or self.objective == "predict_all":
            input_log = ImportDictionaryFromText(self.input_log_name)
            non_constant_dim_in = list(map(bool,input_log.get(self.log_denominator)))
            new_data_in = np.zeros((data_in.shape[0],data_in.shape[1],len(non_constant_dim_in)))
            new_data_in[:,:,non_constant_dim_in]=data_in
            output_log = ImportDictionaryFromText(self.output_log_name)
            non_constant_dim_out = list(map(bool,output_log.get(self.log_denominator)))
            new_data_out = np.zeros((data_out.shape[0],data_out.shape[1],len(non_constant_dim_out)))
            new_data_out[:,:,non_constant_dim_out]=data_out

        return [new_data_in, new_data_out]


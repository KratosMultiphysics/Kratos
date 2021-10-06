import KratosMultiphysics as KM
from KratosMultiphysics.NeuralNetworkApplication.preprocessing_process import PreprocessingProcess
from KratosMultiphysics.NeuralNetworkApplication.masking_process import MaskingProcess
from KratosMultiphysics.NeuralNetworkApplication.data_loading_utilities import ImportDictionaryFromText

import numpy as np


def Factory(settings):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return MaskZerosProcess(settings["parameters"])

class MaskZerosProcess(PreprocessingProcess):

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
            self.log_denominator = "masking_zeros"
        
    def Preprocess(self, data_structure_in, data_structure_out):
        
        data_in = data_structure_in.ExportAsArray()
        data_out = data_structure_out.ExportAsArray()
        
        mask_parameters = KM.Parameters()
        masking_variables = []
        if not self.load_from_log:
            # Setup for input
            if self.objective == 'input':
                try:
                    dimension_variables = range(len(data_in[0]))
                    dataset_length = range(len(data_in[:,0]))
                except IndexError:
                    raise Exception("Trying to mask a one-dimensional array.")
                # Extract which variables are always 0
                for variable in dimension_variables:
                    if all(abs(data_in[i,variable]) < 1e-13 for i in dataset_length):
                        masking_variables.append(variable)
                try:
                    mask_parameters.AddEmptyValue("input_log_name")
                    mask_parameters["input_log_name"].SetString(self.input_log_name)
                except AttributeError:
                    print("No input_log_name defined.")
                    mask_parameters.RemoveValue("input_log_name")

            # Setup for output
            elif self.objective == 'output':
                try:
                    dimension_variables = range(len(data_out[0]))
                    dataset_length = range(len(data_out[:,0]))
                except:
                    raise Exception("Trying to mask one-dimensional array.")
                # Extract which variables are always 0
                for variable in dimension_variables:
                    if all(abs(data_out[i,variable]) < 1e-13 for i in dataset_length):
                        masking_variables.append(variable)
                try:
                    mask_parameters.AddEmptyValue("output_log_name")
                    mask_parameters["output_log_name"].SetString(self.output_log_name)
                except AttributeError:
                        print("No output_log_name defined.")
                        mask_parameters.RemoveValue("output_log_name")
            # else:
            #     raise Exception("Masking objective not supported. Supported objectives are input and output")

        if len(masking_variables) > 0 or self.load_from_log:
            # Load the parameters to the masking process
            mask_parameters.AddEmptyValue("objective")
            mask_parameters["objective"].SetString(self.objective)
            mask_parameters.AddEmptyValue("load_from_log")
            mask_parameters["load_from_log"].SetBool(self.load_from_log)
            mask_parameters.AddEmptyValue("variable_ids")
            mask_parameters["variable_ids"].SetVector(masking_variables)
            mask_parameters.AddEmptyValue("log_denominator")
            mask_parameters["log_denominator"].SetString(self.log_denominator)
            if hasattr(self, 'input_log_name'):
                mask_parameters.AddEmptyValue("input_log_name")
                mask_parameters["input_log_name"].SetString(self.input_log_name)
            if hasattr(self, 'output_log_name'):   
                mask_parameters.AddEmptyValue("output_log_name")
                mask_parameters["output_log_name"].SetString(self.output_log_name)

            # Mask the dataset
            self.mask_process = MaskingProcess(mask_parameters)
            [data_structure_in , data_structure_out] = self.mask_process.Preprocess(data_structure_in, data_structure_out)

        return [data_structure_in, data_structure_out]

    def Invert(self, data_structure_in, data_structure_out):

        data_in = data_structure_in.ExportDataOnly()
        data_out = data_structure_out.ExportDataOnly()

        if self.objective == "input" or self.objective == "predict_input":
            input_log = ImportDictionaryFromText(self.input_log_name)
            input_zero_mask = input_log.get(self.log_denominator)
            for index in input_zero_mask:
                data_in = np.insert(data_in, index, 0.0, axis = 1)
            
        if self.objective == "output" or self.objective == "predict_output":
            output_log = ImportDictionaryFromText(self.output_log_name)
            output_zero_mask = output_log.get(self.log_denominator)
            for index in output_zero_mask:
                data_out = np.insert(data_out, index, 0.0, axis = 1)

        if self.objective == "all" or self.objective == "predict_all":
            input_log = ImportDictionaryFromText(self.input_log_name)
            input_zero_mask = input_log.get(self.log_denominator)
            for index in input_zero_mask:
                data_in = np.insert(data_in, index, 0.0, axis = 1)
            output_log = ImportDictionaryFromText(self.output_log_name)
            output_zero_mask = output_log.get(self.log_denominator)
            for index in output_zero_mask:
                data_out = np.insert(data_out, index, 0.0, axis = 1)

        data_structure_in.UpdateData(data_in)
        data_structure_out.UpdateData(data_out)

        return [data_structure_in, data_structure_out]


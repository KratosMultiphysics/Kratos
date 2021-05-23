import numpy as np
import KratosMultiphysics as KM
from KratosMultiphysics.NeuralNetworkApplication.preprocessing_process import PreprocessingProcess
from KratosMultiphysics.NeuralNetworkApplication.data_loading_utilities import ImportDictionaryFromText, UpdateDictionaryJson, KratosVectorToList


def Factory(settings):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return MaskingProcess(settings["parameters"])

class MaskingProcess(PreprocessingProcess):

    def __init__(self, settings):
        super().__init__(settings)
        """ The default constructor of the class

        Keyword arguments:
        self -- It signifies an instance of a class.
        model -- the container of the different model parts.
        settings -- Kratos parameters containing process settings.
        """
        if not self.load_from_log:
            self.variable_ids = settings["variable_ids"].GetVector()
        self.objective = settings["objective"].GetString()

    def Preprocess(self, data_in, data_out):

        # Check if masking is in the input
        if self.objective == "input": 
            input_log = ImportDictionaryFromText(self.input_log_name)
            # Mask from log
            if self.load_from_log:
                
                if "masking" in input_log:
                    masking_ids = input_log["masking"]
                else:
                    raise Exception("No masking parameters in the input log file.")
            # Mask from JSON
            else:
                masking_ids = list(map(round, KratosVectorToList(self.variable_ids)))
                input_log.update({'masking': masking_ids})
                UpdateDictionaryJson(self.input_log_name, input_log)
            
            data_in = np.delete(data_in, masking_ids, 1)

        # Check if masking is in the output
        elif self.objective == "output": 
            output_log = ImportDictionaryFromText(self.output_log_name)
            # Mask from log
            if self.load_from_log:
                
                if "masking" in output_log:
                    masking_ids = output_log["masking"]
                else:
                    raise Exception("No masking parameters in the output log file.")
            # Mask from JSON
            else:
                masking_ids = list(map(round, KratosVectorToList(self.variable_ids)))
                output_log.update({'masking': masking_ids})
                UpdateDictionaryJson(self.output_log_name, output_log)

            data_out = np.delete(data_out, masking_ids, 1)
     
        else:
            raise Exception("Masking objective not supported. Supported objectives are input and output")    
            
        return [data_in, data_out]

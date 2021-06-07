import numpy as np
import KratosMultiphysics as KM
from KratosMultiphysics.NeuralNetworkApplication.preprocessing_process import PreprocessingProcess
from KratosMultiphysics.NeuralNetworkApplication.data_loading_utilities import ImportDictionaryFromText, UpdateDictionaryJson, KratosVectorToList


def Factory(settings):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return LookbackProcess(settings["parameters"])

class LookbackProcess(PreprocessingProcess):

    def __init__(self, settings):
        super().__init__(settings)
        """ The default constructor of the class

        Keyword arguments:
        self -- It signifies an instance of a class.
        model -- the container of the different model parts.
        settings -- Kratos parameters containing process settings.
        """
        if not self.load_from_log:
            self.lookback = settings["lookback"].GetInt()

    def Preprocess(self, data_in, data_out):

        new_data_in, new_data_out = [], []

        for i in range(len(data_in)-self.lookback-1):
            value = data_out[i:(i+self.lookback),0]
            new_data_in.append(value)
            new_data_out.append(data_out[i+self.lookback]) 
        new_data_in = np.array(new_data_in)
        new_data_out = np.array(new_data_out)  
        new_data_in = np.reshape(new_data_in, (new_data_in.shape[0],new_data_in.shape[1],1))
        return [new_data_in, new_data_out]

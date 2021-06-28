import numpy as np
import KratosMultiphysics as KM
from KratosMultiphysics.NeuralNetworkApplication.preprocessing_process import PreprocessingProcess
from KratosMultiphysics.NeuralNetworkApplication.data_loading_utilities import ImportDictionaryFromText, UpdateDictionaryJson, KratosVectorToList


def Factory(settings):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return DivideDatasetProcess(settings["parameters"])

class DivideDatasetProcess(PreprocessingProcess):

    def __init__(self, settings):
        super().__init__(settings)
        """ The default constructor of the class

        Keyword arguments:
        self -- It signifies an instance of a class.
        model -- the container of the different model parts.
        settings -- Kratos parameters containing process settings.
        """
        try:
            self.proportion = settings["proportion"].GetDouble()
        except AttributeError:
            print("No proportion indicated for the division of the dataset. The chosen value is 0.0.")
        
        if (self.proportion >= 1.0) or (self.proportion <= 0.0):
            print("The division proportion must be a value between 0.0 and 1.0.")

    def Preprocess(self, data_in, data_out):
        """ This method records the lookback of timeseries and saves them. 
            Normally applied before LSTM and RNN layers. """

        # new_data_in, new_data_out = [], []

        new_data_in = data_in[:np.int64(len(data_in)*self.proportion)]
        new_data_out = data_out[:np.int64(len(data_out)*self.proportion)]
        # for i in range(np.int64(len(data_in)*self.proportion)):
        #     new_data_in.append(data_in[i])
        #     new_data_out.append(data_out[i])
        return [new_data_in, new_data_out]

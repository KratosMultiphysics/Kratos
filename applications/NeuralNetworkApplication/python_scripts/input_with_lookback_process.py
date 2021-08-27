import numpy as np
import KratosMultiphysics as KM
from KratosMultiphysics.NeuralNetworkApplication.preprocessing_process import PreprocessingProcess
from KratosMultiphysics.NeuralNetworkApplication.input_dataclasses import ListDataWithLookback, DataWithLookback


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

    def Preprocess(self, data_structure_in, data_structure_out):
        """ This method records the lookback of timeseries and saves them. 
            Normally applied before LSTM and RNN layers. """

        if not isinstance(data_structure_in, ListDataWithLookback) or not isinstance(data_structure_in, DataWithLookback):
            new_data_structure_in = ListDataWithLookback(lookback_index = self.lookback)
            new_data_structure_in.ExtendFromNeuralNetworkData(data_structure_in)
        

        for i in reversed(range(self.lookback)):
            
            data_in = data_structure_out.ExportAsArray()
            new_data_in = np.zeros_like(data_in)
            new_data_in[i+1:] = data_in[:-i-1]
            # if len(new_data_in.shape) == 1: 
            #     new_data_in = np.reshape(new_data_in, (1, 1,new_data_in.shape[0]))

            # if len(new_data_in.shape) == 2: 
            #     new_data_in = np.reshape(new_data_in, (new_data_in.shape[0],1, new_data_in.shape[1]))

            new_data_structure_in.CheckLookbackAndUpdate(new_data_in)

        return [new_data_structure_in, data_structure_out]
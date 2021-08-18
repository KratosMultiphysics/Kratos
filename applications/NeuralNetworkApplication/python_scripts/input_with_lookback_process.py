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
        """ This method records the lookback of timeseries and saves them. 
            Normally applied before LSTM and RNN layers. """

        new_data_in, new_data_out = [], []
        if data_in == None:
            data_in = np.array([0.0])
        if data_out == None:
            data_out = np.array([0.0])

        lookback_indexes = len(data_in)-self.lookback-1
        if lookback_indexes < 1:
            initial_vector_shape_input = data_in.shape
            initial_vector_shape_output = data_out.shape
            new_data_in = data_in
            new_data_out = data_out
            for i in range(-lookback_indexes):
                new_data_in = np.concatenate((new_data_in,np.zeros(initial_vector_shape_input)))
                new_data_out = np.concatenate((new_data_out,np.zeros(initial_vector_shape_output)))

        for i in range(len(data_in)-self.lookback-1):
            value = data_out[i:(i+self.lookback)]
            new_data_in.append(np.concatenate((data_in[i], value.T[0])))
            new_data_out.append(data_out[i+self.lookback]) 
        new_data_in = np.array(new_data_in)
        new_data_out = np.array(new_data_out)

        if len(new_data_in.shape) == 1: 
            new_data_in = np.reshape(new_data_in, (1, 1,new_data_in.shape[0]))

        if len(new_data_in.shape) == 2: 
            new_data_in = np.reshape(new_data_in, (new_data_in.shape[0],1, new_data_in.shape[1]))
        print(new_data_in)
        return [new_data_in, new_data_out]

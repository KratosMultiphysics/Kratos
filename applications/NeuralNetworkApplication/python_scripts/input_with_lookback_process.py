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
            try:
                self.record = settings["record"].GetBool()
            except RuntimeError:
                self.record = False
            try:
                self.timesteps_as_features = settings["timesteps_as_features"].GetBool()
            except RuntimeError:
                self.timesteps_as_features = False
            try:
                self.features_as_timesteps = settings["features_as_timesteps"].GetBool()
            except RuntimeError:
                self.features_as_timesteps = False
            try:
                self.only_input = settings["only_input"].GetBool()
            except RuntimeError:
                self.only_input = False
            

    def Preprocess(self, data_structure_in, data_structure_out):
        """ This method records the lookback of timeseries and saves them. 
            Normally applied before LSTM and RNN layers. """

        if not isinstance(data_structure_in, ListDataWithLookback) or not isinstance(data_structure_in, DataWithLookback):
            new_data_structure_in = ListDataWithLookback(lookback_index = self.lookback, record_data=self.record, 
                                                        timesteps_as_features = self.timesteps_as_features, 
                                                        features_as_timesteps = self.features_as_timesteps)
            new_data_structure_in.ExtendFromNeuralNetworkData(data_structure_in)
        

        for i in reversed(range(self.lookback)):
            
            if not self.only_input:
                data_in_lookback = data_structure_out.ExportAsArray()
                new_data_in = np.zeros_like(data_in_lookback)
                new_data_in[i+1:] = data_in_lookback[:-i-1]
                new_data_structure_in.CheckLookbackAndUpdate(new_data_in)

            if self.record:
                data_in_record = data_structure_in.ExportAsArray()
                new_data_in_record = np.zeros_like(data_in_record)
                new_data_in_record[i+1:] = data_in_record[:-i-1]
                new_data_structure_in.CheckRecordAndUpdate(new_data_in_record)

        return [new_data_structure_in, data_structure_out]
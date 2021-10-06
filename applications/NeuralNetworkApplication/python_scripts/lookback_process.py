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
        if settings.Has("stack_input"):
            self.stack_input = settings["stack_input"].GetBool()
        else:
            self.stack_input = False
        try:
            self.timesteps_as_features = settings["timesteps_as_features"].GetBool()
        except RuntimeError:
            self.timesteps_as_features = False
        try:
            self.features_as_timesteps = settings["features_as_timesteps"].GetBool()
        except RuntimeError:
            self.features_as_timesteps = False
        try:
            self.flatten_lookback = settings["flatten_lookback"].GetBool()
        except RuntimeError:
            self.flatten_lookback = True

    def Preprocess(self, data_structure_in, data_structure_out):
        """ This method records the lookback of timeseries and saves them. 
            Normally applied before LSTM and RNN layers. """

        if not isinstance(data_structure_in, ListDataWithLookback) or not isinstance(data_structure_in, DataWithLookback):
            new_data_structure_in = ListDataWithLookback(lookback_index = self.lookback,
                                                        only_lookback = not self.stack_input,
                                                        timesteps_as_features = self.timesteps_as_features, 
                                                        features_as_timesteps = self.features_as_timesteps)
            new_data_structure_in.ExtendFromNeuralNetworkData(data_structure_in)
        
        data_in = data_structure_out.ExportAsArray()

        for i in reversed(range(self.lookback)):
            new_data_in = np.zeros_like(data_in)
            new_data_in[i+1:] = data_in[:-i-1]
            new_data_structure_in.CheckLookbackAndUpdate(new_data_in)

        if self.flatten_lookback:
            new_data_structure_in.FlattenTo2DLookback()
        new_data_structure_in.SetOnlyLookback(False) if self.stack_input else new_data_structure_in.SetOnlyLookback(True)


        return [new_data_structure_in, data_structure_out]

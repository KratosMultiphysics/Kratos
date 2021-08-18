import numpy as np
from dataclasses import dataclass, field

@dataclass
class NeuralNetworkData:
    data: np.ndarray = None

    def UpdateData(self, new_data):
        self.data = new_data


@dataclass
class DataWithLookback(NeuralNetworkData):
    lookback_index: int = 1
    lookback_data: np.ndarray = field(init=False)
    lookback_state: bool = False

    def ExportFull(self):
        return np.concatenate((self.input, self.lookback_data))
    
    def UpdateLookbackLast(self, new_data):
        new_array = np.zeros_like(self.lookback_data)
        if self.lookback_index > 1:
            new_array[:-2] = self.lookback_data[1:]
            new_array[-1] = new_data
        else:
            new_array = new_data
        self.lookback_data = new_array

    def UpdateLookbackAll(self, new_array):
        self.lookback_data = new_array
        self.lookback_index = new_array.shape[0]
        self.lookback_state = True

    def CheckLookbackAndUpdate(self, new_data):
        if self.lookback_state:
            self.UpdateLookbackLast(new_data)
        else:
            self.lookback_data = new_data
            if self.lookback_index > 1:
                for i in range(self.lookback_index-1):
                    self.lookback_index = np.concatenate((self.lookback_data, np.zeros_like(new_data)))
            self.lookback_state = True


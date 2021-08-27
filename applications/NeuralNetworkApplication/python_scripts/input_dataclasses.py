import numpy as np
from dataclasses import dataclass, field
from typing import List

@dataclass
class NeuralNetworkData:
    data: np.ndarray = None

    def UpdateData(self, new_data):
        self.data = new_data
    
    def ExportAsArray(self):
        return self.data

    def __len__(self):
        return len(self.data)
    
    def __getitem__(self, key):
        return self.data[key]

@dataclass
class DataWithLookback(NeuralNetworkData):
    lookback_index: int = 1
    lookback_data: np.ndarray = field(init=False)
    lookback_state: bool = False
    only_lookback: bool = False

    def ExportAsArray(self):
        if self.data.all() != None:
            if hasattr(self, "lookback_data"):
                if self.only_lookback:
                    new_array = np.squeeze(self.lookback_data)
                else:
                    new_array = np.concatenate((self.data, np.squeeze(self.lookback_data)))
                try:
                    new_array = np.reshape(new_array, (1, new_array.shape[0], new_array.shape[1]))
                except IndexError:
                    new_array = np.reshape(new_array, (1, new_array.shape[0], 1))
                return new_array
            else:
                try:
                    new_array = np.reshape(self.data, (1, self.data.shape[0], self.data.shape[1]))
                except IndexError:
                    new_array = np.reshape(self.data, (1, self.data.shape[0], 1))
                return new_array
        else:
            new_array = self.lookback_data
            try:
                new_array = np.reshape(new_array, (1, new_array.shape[0], new_array.shape[1]))
            except IndexError:
                new_array = np.reshape(new_array, (1, new_array.shape[0], 1))
            return new_array
    
    def UpdateLookbackLast(self, new_data):
        new_array = np.zeros_like(self.lookback_data)
        if self.lookback_index > 1:
            new_array[:-1] = self.lookback_data[1:]
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
                    self.lookback_data = np.vstack([np.zeros_like(new_data), self.lookback_data])

            self.lookback_state = True
    
    def ExtendFromNeuralNetworkData(self, data_structure):
        self.data = data_structure.data
    
    def SetOnlyLookback(self, switch_to = True):
        self.only_lookback = switch_to

@dataclass
class ListNeuralNetworkData:
    data_array: List[NeuralNetworkData] = field(default_factory=list)

    def AddToList(self, new_array):
        if isinstance(new_array, NeuralNetworkData):
            self.data_array.append(new_array)
        elif isinstance(new_array, np.ndarray) or isinstance(new_array, (list,float)):
            self.data_array.append(NeuralNetworkData(new_array))
        else:
            self.data_array.append(NeuralNetworkData(np.ndarray(new_array)))

    def ExportAsArray(self):
        array = []
        for i in range(len(self.data_array)):
            try:
                if len(array)>0 or not isinstance(array, list):
                    new_timestep = self.data_array[i].ExportAsArray()
                    new_timestep = np.reshape(new_timestep,  tuple((1,))+new_timestep.shape)
                    array = np.vstack([array, new_timestep]) 
                else:
                    array = self.data_array[i].ExportAsArray()
                    array = np.reshape(array,tuple((1,)+array.shape))
            except TypeError:
                array = np.reshape(array,array.shape + tuple((1,)))
                array = np.vstack([array, self.data_array[i].ExportAsArray()]) 
                
        array = np.squeeze(array)
        return array

    def ImportFromArray(self, array):
        for i in range(array.shape[0]):
            self.AddToList(array[i])

    def UpdateData(self, array):
        for i in range(len(self.data_array)):
            self.data_array[i].UpdateData(array[i])
    
    def __len__(self):
        return len(self.data_array)

    def __getitem__(self, key):
        return self.data_array[key]
    
    def KeepPart(self,key):
        self.data_array = self.data_array[key]

@dataclass
class ListDataWithLookback(ListNeuralNetworkData):
    data_array: List[DataWithLookback] = field(default_factory=list)
    lookback_index: int = 1
    only_lookback: bool = False

    def AddToList(self, new_array):
        if isinstance(new_array, DataWithLookback):
            self.data_array.append(new_array)
        elif isinstance(new_array, np.ndarray) or isinstance(new_array, list):
            self.data_array.append(DataWithLookback(new_array, self.lookback_index))
        else:
            self.data_array.append(DataWithLookback(np.ndarray(new_array), self.lookback_index))
    
    def ExportLookbackAsArray(self):
        array = []
        for i in range(len(self.data_array)):
            try:
                if len(array)>0 or not isinstance(array, list):
                    array = np.vstack([array, self.data_array[i].ExportLookbackAsArray()]) 
                else:
                    array = self.data_array[i].ExportLookbackAsArray()
            except TypeError:
                array = np.vstack([array, self.data_array[i].ExportLookbackAsArray()]) 
        return array

    def UpdateLookbackLast(self, array):
        for i in range(array.shape[0]):
            self.data_array[i].UpdateLookbackLast(array[i])

    def UpdateLookbackAll(self, array):
       for i in range(array.shape[0]):
            self.data_array[i].UpdateLookbackAll(array[i])

    def CheckLookbackAndUpdate(self, array):
        for i in range(array.shape[0]):
            self.data_array[i].CheckLookbackAndUpdate(array[i])
    
    def ExtendFromNeuralNetworkData(self, data_structure):
        for i in range(len(data_structure)):
            new_entry = DataWithLookback(lookback_index=self.lookback_index)
            new_entry.ExtendFromNeuralNetworkData(data_structure[i])
            self.AddToList(new_entry)
    
    def SetOnlyLookback(self, switch_to = True):
        for i in range(len(self.data_array)):
            self.data_array[i].SetOnlyLookback(switch_to)
    
    





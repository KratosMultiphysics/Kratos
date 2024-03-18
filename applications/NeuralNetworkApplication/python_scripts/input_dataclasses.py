import numpy as np
from dataclasses import dataclass, field
from typing import List


@dataclass
class NeuralNetworkData:
    data: np.ndarray = None

    def UpdateData(self, new_data):
        if isinstance(new_data, NeuralNetworkData):
            self.data = new_data.data
        else:
            self.data = new_data

    def ExportAsArray(self):
        try:
            """for converting pytorch predictions to numpy"""
            self.data = self.data.detach().numpy()
            return self.data
        except:
            return self.data

    def ExportElementAsArray(self, element_id):
        return self.ExportAsArray()

    def ExportDataOnly(self):
        return self.ExportAsArray()

    def __len__(self):
        if self.data is None:
            return 0
        else:
            try:
                return len(self.data)
            except TypeError:
                return self.data.size

    def __getitem__(self, key):
        return self.data[key]

    def GetDataShape(self):
        return self.data.shape


@dataclass
class DataWithLookback(NeuralNetworkData):
    lookback_index: int = 1
    lookback_data: np.ndarray = field(init=False)
    lookback_state: bool = False
    only_lookback: bool = False
    timesteps_as_features: bool = False
    features_as_timesteps: bool = False
    record_data: bool = False
    reorder_partitions: int = 1

    def ExportAsArray(self):
        if self.data.all() != None:
            if hasattr(self, "lookback_data"):
                if self.only_lookback:
                    new_array = np.squeeze(self.lookback_data)
                else:
                    if self.data.size == 1:
                        new_array = np.concatenate(
                            (
                                np.reshape(np.squeeze(self.data), (1,)),
                                np.squeeze(self.lookback_data),
                            )
                        )
                    else:
                        try:
                            new_array = np.hstack((self.data, self.lookback_data))
                        except ValueError:
                            raise Exception(
                                "Dimension mismatch between the lookback and the data. Check if the dimensions are compatible."
                            )
                try:
                    new_array = np.reshape(
                        new_array, (1, new_array.shape[0], new_array.shape[1])
                    )
                except IndexError:
                    # if self.lookback_index == 1:
                    #     new_array = np.reshape(new_array, (1, 1, new_array.shape[0]))
                    # else:
                    new_array = np.reshape(new_array, (1, new_array.shape[0], 1))
                if self.reorder_partitions > 1:
                    new_array = self.Reorder(new_array, self.reorder_partitions)
                return new_array
            else:
                try:
                    new_array = np.reshape(
                        self.data, (1, self.data.shape[0], self.data.shape[1])
                    )
                except IndexError:
                    new_array = np.reshape(self.data, (1, self.data.shape[0], 1))
                if self.reorder_partitions > 1:
                    new_array = self.Reorder(new_array, self.reorder_partitions)
                return new_array
        else:
            new_array = self.lookback_data
            try:
                new_array = np.reshape(
                    new_array, (1, new_array.shape[0], new_array.shape[1])
                )
            except IndexError:
                new_array = np.reshape(new_array, (1, new_array.shape[0], 1))
            if self.reorder_partitions > 1:
                new_array = self.Reorder(new_array, self.reorder_partitions)
            return new_array

    def Reorder(self, array, number_of_partitions):
        partitions = np.dsplit(array, number_of_partitions)
        new_array = array
        for i in range(len(partitions)):
            new_array[..., i::number_of_partitions] = partitions[i]

        return new_array

    def ExportDataOnly(self):
        return self.data

    def UpdateLookbackLast(self, new_data):
        new_array = np.zeros_like(self.lookback_data)
        if self.lookback_index > 1:
            # if len(new_data.shape)>1:
            #     new_array[:-new_data.shape[0]] = self.lookback_data[new_data.shape[0]:]
            #     new_array[-new_data.shape[0]:] = new_data
            # else:
            new_array[:-1] = self.lookback_data[1:]
            new_array[-1] = new_data
        else:
            new_array = new_data
        self.lookback_data = new_array

    def UpdateRecordLast(self, new_data):
        new_array = np.zeros_like(self.data)
        if self.lookback_index > 1:
            new_array[:-1] = self.data[1:]
            new_array[-1] = new_data
        else:
            new_array = new_data
        self.data = new_array

    def UpdateLookbackAll(self, new_array):
        self.lookback_data = new_array
        try:
            if len(new_array.shape) > 1:
                self.lookback_index = new_array.shape[0]
        except AttributeError:
            pass
        self.lookback_state = True

    def CheckLookbackAndUpdate(self, new_data):
        if self.lookback_state:
            self.UpdateLookbackLast(new_data)
        else:
            if self.lookback_index > 1:
                self.lookback_data = np.stack(
                    list(np.zeros_like(new_data) for _ in range(self.lookback_index))
                )
                self.lookback_data[-1] = new_data
            else:
                self.lookback_data = new_data

            self.lookback_state = True

    def FlattenTo2DLookback(self):
        if len(self.lookback_data.shape) > 2:
            self.lookback_data = np.reshape(
                self.lookback_data,
                (
                    self.lookback_index,
                    int(self.lookback_data.size / self.lookback_index),
                ),
            )
        else:
            pass

    def CheckRecordAndUpdate(self, new_data):
        if self.record_data:
            self.UpdateRecordLast(new_data)
        else:
            if self.lookback_index > 1:
                self.data = np.stack(
                    list(np.zeros_like(new_data) for _ in range(self.lookback_index))
                )
                self.data[-1] = new_data
            else:
                self.data = new_data

            self.record_data = True

    # TODO Input with a set value as entry

    def ExtendFromNeuralNetworkData(self, data_structure):
        self.data = data_structure.data

    def SetOnlyLookback(self, switch_to=True):
        self.only_lookback = switch_to

    def GetLookbackIndex(self):
        return self.lookback_index

    def GetLookbackShape(self):
        if self.lookback_index == 1:
            return tuple((1, len(self.lookback_data)))
        else:
            try:
                return self.lookback_data.shape
            except AttributeError:
                return tuple((self.lookback_index, 0))

    def ExportElementAsArray(self, id):
        array = self.ExportAsArray()

        try:
            if self.only_lookback and (
                self.timesteps_as_features or self.features_as_timesteps
            ):
                feature_length = self.GetLookbackShape()[1] * self.lookback_index
            elif self.only_lookback:
                feature_length = self.GetLookbackShape()[1]
            elif self.timesteps_as_features or self.features_as_timesteps:
                feature_length = (
                    self.GetDataShape()[1]
                    + self.GetLookbackShape()[1] * self.lookback_index
                )
            else:
                feature_length = self.GetDataShape()[1] + self.GetLookbackShape()[1]
        except IndexError:
            if self.timesteps_as_features or self.features_as_timesteps:
                feature_length = (
                    self.GetDataShape()[0]
                    + self.GetLookbackShape()[1] * self.lookback_index
                )
            else:
                feature_length = self.GetDataShape()[0] + self.GetLookbackShape()[1]

        if self.timesteps_as_features:
            array = np.reshape(array, (1, 1, feature_length))
        elif self.features_as_timesteps:
            array = np.reshape(array, (1, feature_length, 1))
        else:
            array = np.reshape(array, (1, self.lookback_index, feature_length))
        return array

    def SetTimestepsAsFeatures(self, switch_to=True):
        self.timesteps_as_features = switch_to

    def SetFeaturesAsTimesteps(self, switch_to=True):
        self.features_as_timesteps = switch_to


@dataclass
class ListNeuralNetworkData:
    data_array: List[NeuralNetworkData] = field(default_factory=list)

    def AddToList(self, new_array):
        if isinstance(new_array, NeuralNetworkData):
            self.data_array.append(new_array)
        elif isinstance(new_array, np.ndarray) or isinstance(
            new_array, (list, float, np.floating)
        ):
            self.data_array.append(NeuralNetworkData(new_array))
        else:
            self.data_array.append(NeuralNetworkData(np.ndarray(new_array)))

    def ExportAsArray(self):
        array = []
        for i in range(len(self.data_array)):
            try:
                if len(array) > 0 or not isinstance(array, list):
                    new_timestep = self.data_array[i].ExportAsArray()
                    array = np.vstack([array, new_timestep])
                else:
                    array = self.data_array[i].ExportAsArray()
            except TypeError:
                array = np.reshape(array, array.shape + tuple((1,)))
                array = np.vstack([array, self.data_array[i].ExportAsArray()])

        # THIS CAN BE MORE EFFICIENT, RESHAPE TAKES A LOT OF TIME
        if not self.data_array[0].GetDataShape() == ():
            try:
                array = np.reshape(
                    array,
                    (
                        len(self.data_array),
                        self.data_array[0].GetDataShape()[0],
                        self.data_array[0].GetDataShape()[1],
                    ),
                )
            except IndexError:
                pass
        return array

    def ExportElementAsArray(self, element_id):
        return self.data_array[element_id].ExportElementAsArray(element_id)

    def ExportDataOnly(self):
        return self.ExportAsArray()

    def ImportFromArray(self, array):
        for i in range(array.shape[0]):
            self.AddToList(array[i])

    def UpdateData(self, array):
        for i in range(len(self.data_array)):
            self.data_array[i].UpdateData(array[i])

    def CopyLast(self):
        self.AddToList(self.data_array[-1])

    def __len__(self):
        return len(self.data_array)

    def __getitem__(self, key):
        return self.data_array[key]

    def KeepPart(self, key):
        self.data_array = self.data_array[key]

    def GetDataShape(self):
        return self.data_array[0].GetDataShape()


@dataclass
class ListDataWithLookback(ListNeuralNetworkData):
    data_array: List[DataWithLookback] = field(default_factory=list)
    lookback_index: int = 1
    only_lookback: bool = False
    record_data: bool = False
    timesteps_as_features: bool = False
    features_as_timesteps: bool = False

    def AddToList(self, new_array):
        if isinstance(new_array, DataWithLookback):
            self.data_array.append(new_array)
        elif isinstance(new_array, np.ndarray) or isinstance(new_array, list):
            self.data_array.append(DataWithLookback(new_array, self.lookback_index))
        else:
            self.data_array.append(
                DataWithLookback(np.ndarray(new_array), self.lookback_index)
            )

    def ExportAsArray(self):
        array = []
        for i in range(len(self.data_array)):
            try:
                if len(array) > 0 or not isinstance(array, list):
                    new_timestep = self.data_array[i].ExportAsArray()
                    array = np.vstack([array, new_timestep])
                else:
                    array = self.data_array[i].ExportAsArray()
            except TypeError:
                array = np.reshape(array, array.shape + tuple((1,)))
                array = np.vstack([array, self.data_array[i].ExportAsArray()])
        out_array = self.__Reshape(array)
        return out_array

    def ExportDataOnly(self):
        array = []
        for i in range(len(self.data_array)):
            try:
                if len(array) > 0 or not isinstance(array, list):
                    new_timestep = self.data_array[i].ExportDataOnly()
                    array = np.vstack([array, new_timestep])
                else:
                    array = self.data_array[i].ExportDataOnly()
            except TypeError:
                array = np.reshape(array, array.shape + tuple((1,)))
                array = np.vstack([array, self.data_array[i].ExportDataOnly()])

        out_array = self.__Reshape(array, data_only=True)
        return out_array

    def ExportElementAsArray(self, element_id):
        array = self.data_array[element_id].ExportElementAsArray(element_id)

        # try:
        #     if self.only_lookback and (self.timesteps_as_features or self.features_as_timesteps):
        #         feature_length = self.data_array[0].GetLookbackShape()[1] * self.lookback_index
        #     elif self.only_lookback:
        #         feature_length = self.data_array[0].GetLookbackShape()[1]
        #     elif self.timesteps_as_features or self.features_as_timesteps:
        #         feature_length = self.data_array[0].GetDataShape()[1]+self.data_array[0].GetLookbackShape()[1] * self.lookback_index
        #     else:
        #         feature_length = self.data_array[0].GetDataShape()[1]+self.data_array[0].GetLookbackShape()[1]
        # except IndexError:
        #     if self.timesteps_as_features or self.features_as_timesteps:
        #         feature_length = self.data_array[0].GetDataShape()[0]+self.data_array[0].GetLookbackShape()[1] * self.lookback_index
        #     else:
        #         feature_length = self.data_array[0].GetDataShape()[0]+self.data_array[0].GetLookbackShape()[1]

        # if self.timesteps_as_features:
        #     array = np.reshape(array, (1,1, feature_length))
        # elif self.features_as_timesteps:
        #     array = np.reshape(array, (1, feature_length,1))
        # else:
        #     array = np.reshape(array, (1,self.lookback_index, feature_length))
        return array

    def ExportLookbackAsArray(self):
        array = []
        for i in range(len(self.data_array)):
            self.data_array[i].SetOnlyLookback()
            try:
                if len(array) > 0 or not isinstance(array, list):
                    array = np.vstack([array, self.data_array[i].ExportAsArray()])
                else:
                    array = self.data_array[i].ExportAsArray()
            except TypeError:
                array = np.vstack([array, self.data_array[i].ExportAsArray()])
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

    def CheckRecordAndUpdate(self, array):
        for i in range(array.shape[0]):
            self.data_array[i].CheckRecordAndUpdate(array[i])

    def ExtendFromNeuralNetworkData(self, data_structure):
        for i in range(len(data_structure)):
            new_entry = DataWithLookback(lookback_index=self.lookback_index)
            if isinstance(data_structure, NeuralNetworkData):
                new_entry.ExtendFromNeuralNetworkData(data_structure)
            else:
                new_entry.ExtendFromNeuralNetworkData(data_structure[i])
            self.AddToList(new_entry)

    def FlattenTo2DLookback(self):
        for i in range(len(self.data_array)):
            self.data_array[i].FlattenTo2DLookback()

    def SetOnlyLookback(self, switch_to=True):
        for i in range(len(self.data_array)):
            self.data_array[i].SetOnlyLookback(switch_to)

    def SetTimestepsAsFeatures(self, switch_to=True):
        self.timesteps_as_features = switch_to
        for entry in self.data_array:
            entry.SetTimestepsAsFeatures(switch_to)

    def SetFeaturesAsTimesteps(self, switch_to=True):
        self.features_as_timesteps = switch_to
        for entry in self.data_array:
            entry.SetFeaturesAsTimesteps(switch_to)

    def SetReorder(self, number_of_partitions):
        for entry in self.data_array:
            entry.reorder_partitions = number_of_partitions

    def GetLookbackIndex(self):
        return self.lookback_index

    def GetLookbackShape(self):
        return self.data_array[0].GetLookbackShape()

    def AddDataAndUpdate(self, new_data, new_lookback):
        self.CopyLast()
        self.data_array[-1].UpdateData(new_data)
        self.data_array[-1].CheckLookbackAndUpdate(new_lookback)

    def __Reshape(self, array, data_only=False):
        # THIS CAN BE MORE EFFICIENT, RESHAPE TAKES A LOT OF TIME

        # Calculate the length of the lookback (flattened or not)
        if self.timesteps_as_features or self.features_as_timesteps:
            lookback_length = (
                self.data_array[0].GetLookbackShape()[1] * self.lookback_index
            )
            # Calculate the lenght of the data (with record of the lookback or not)
            try:
                if self.record_data:
                    data_length = (
                        self.data_array[0].GetDataShape()[1] * self.lookback_index
                    )
                else:
                    data_length = self.data_array[0].GetDataShape()[1]
            except IndexError:
                if self.record_data:
                    data_length = (
                        self.data_array[0].GetDataShape()[0] * self.lookback_index
                    )
                else:
                    data_length = self.data_array[0].GetDataShape()[0]
        else:
            if self.lookback_index == 1:
                lookback_length = self.data_array[0].GetLookbackShape()[0]
            else:
                lookback_length = self.data_array[0].GetLookbackShape()[1]
            try:
                data_length = self.data_array[0].GetDataShape()[1]
            except IndexError:
                data_length = self.data_array[0].GetDataShape()[0]

        # Calculate the feature length
        if data_only:
            feature_length = data_length
        elif self.only_lookback:
            feature_length = lookback_length
        else:
            feature_length = data_length + lookback_length
        # Reshape for output
        if not self.data_array[0].GetDataShape() == ():
            if self.timesteps_as_features:
                out_array = np.reshape(array, (len(self.data_array), 1, feature_length))
            elif self.features_as_timesteps:
                out_array = np.reshape(array, (len(self.data_array), feature_length, 1))
            elif data_only and not self.record_data:
                out_array = np.reshape(array, (len(self.data_array), feature_length))
            else:
                out_array = np.reshape(
                    array, (len(self.data_array), self.lookback_index, feature_length)
                )
        return out_array

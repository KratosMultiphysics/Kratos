"""A set of classes for performing HDF5 output.

BSD license: HDF5Application/license.txt
"""
from abc import ABCMeta, abstractmethod
import KratosMultiphysics
import KratosMultiphysics.HDF5Application as KratosHDF5


class FileFactory(metaclass=ABCMeta):

    @abstractmethod
    def Open(self, file_name): pass


class OutputObject(metaclass=ABCMeta):

    @abstractmethod
    def Execute(self, model_part, hdf5_file): pass


class HDF5SerialFileFactory(FileFactory):

    def __init__(self, file_access_mode="truncate", echo_level=0):
        self._file_access_mode = file_access_mode
        self._echo_level = echo_level
    
    def Open(self, file_name):
        params = KratosMultiphysics.Parameters("""{
            "file_name" : "",
            "file_access_mode" : "",
            "file_driver" : "sec2",
            "echo_level" : 0
        }""")
        params["file_name"].SetString(file_name)
        params["file_access_mode"].SetString(self._file_access_mode)
        params["echo_level"].SetInt(self._echo_level)
        return KratosHDF5.HDF5FileSerial(params)


class HDF5ParallelFileFactory(FileFactory):

    def __init__(self, file_access_mode="truncate", echo_level=0):
        self._file_access_mode = file_access_mode
        self._echo_level = echo_level
    
    def Open(self, file_name):
        params = KratosMultiphysics.Parameters("""{
            "file_name" : "",
            "file_access_mode" : "",
            "file_driver" : "mpio",
            "echo_level" : 0
        }""")
        params["file_name"].SetString(file_name)
        params["file_access_mode"].SetString(self._file_access_mode)
        params["echo_level"].SetInt(self._echo_level)
        return KratosHDF5.HDF5FileParallel(params)


class ModelPartOutput(OutputObject):
    """A class for writing a model part."""

    def Execute(self, model_part, hdf5_file):
        params = KratosMultiphysics.Parameters("""{
            "prefix" : "/ModelData"
        }""")
        KratosHDF5.HDF5ModelPartIO(params, hdf5_file).WriteModelPart(model_part)


class NodalResultsOutput(OutputObject):
    """A class for writing nodal results of a model part."""

    def __init__(self, list_of_variables):
        self._list_of_variables = list_of_variables

    def Execute(self, model_part, hdf5_file):
        params = KratosMultiphysics.Parameters("""{
            "prefix" : "/ResultsData",
            "partitioned": false,
            "list_of_variables": []
        }""")
        for var in self._list_of_variables:
            params["list_of_variables"].Append(var)
        KratosHDF5.HDF5NodalSolutionStepDataIO(params, hdf5_file).WriteNodalResults(model_part.Nodes, 0)


class PartitionedModelPartOutput(OutputObject):
    """A class for writing a partitioned model part."""

    def Execute(self, model_part, hdf5_file):
        params = KratosMultiphysics.Parameters("""{
            "prefix" : "/ModelData"
        }""")
        KratosHDF5.HDF5PartitionedModelPartIO(params, hdf5_file).WriteModelPart(model_part)


class PartitionedNodalResultsOutput(OutputObject):
    """A class for writing nodal results of a partitioned model part."""

    def __init__(self, list_of_variables):
        self._list_of_variables = list_of_variables

    def Execute(self, model_part, hdf5_file):
        params = KratosMultiphysics.Parameters("""{
            "prefix" : "/ResultsData",
            "partitioned": true,
            "list_of_variables": []
        }""")
        for var in self._list_of_variables:
            params["list_of_variables"].Append(var)
        KratosHDF5.HDF5NodalSolutionStepDataIO(params, hdf5_file).WriteNodalResults(model_part.Nodes, 0)


class TemporalOutputProcess(KratosMultiphysics.Process):
    """A process for writing temporal simulation results."""

    def __init__(self, model_part, hdf5_file_factory, output_time_frequency, output_step_frequency=0):
        KratosMultiphysics.Process.__init__(self)
        self._model_part = model_part
        self._hdf5_file_factory = hdf5_file_factory
        self._output_time_frequency = output_time_frequency
        self._output_step_frequency = output_step_frequency
        self._list_of_outputs = []
        self._time_tag_precision = 4
        self._output_time = 0.0
        self._output_step = 0

    def AddOutputObject(self, output):
        self._list_of_outputs.append(output)

    def ExecuteBeforeSolutionLoop(self):
        self._output_time = 0.0
        self._output_step = 0

    def ExecuteFinalizeSolutionStep(self):
        delta_time = self._model_part.ProcessInfo[KratosMultiphysics.DELTA_TIME]
        self._output_time += delta_time
        self._output_step += 1
        if self._output_time >= self._output_time_frequency or self._output_step == self._output_step_frequency:
            hdf5_file = self._hdf5_file_factory.Open(self._get_current_file_name())
            for output in self._list_of_outputs:
                output.Execute(self._model_part, hdf5_file)
            self._output_time = 0.0
            self._output_step = 0

    def _get_current_file_name(self):
        fmt = "{:." + str(self._time_tag_precision) + "f}"
        time_tag = "-" + fmt.format(self._model_part.ProcessInfo[KratosMultiphysics.TIME])
        return self._model_part.Name + time_tag + ".h5"

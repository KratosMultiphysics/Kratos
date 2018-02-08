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

    def __init__(self, settings):
        default_settings = KratosMultiphysics.Parameters("""
            {
                "file_access_mode" : "truncate",
                "file_driver" : "sec2",
                "echo_level" : 0
            }
            """)
        self.settings = settings.Clone()
        self.settings.ValidateAndAssignDefaults(default_settings)
        self.settings.AddEmptyValue("file_name")
    
    def Open(self, file_name):
        self.settings["file_name"].SetString(file_name)
        return KratosHDF5.HDF5FileSerial(settings)


class HDF5ParallelFileFactory(FileFactory):

    def __init__(self, settings):
        default_settings = KratosMultiphysics.Parameters("""
            {
                "file_access_mode" : "truncate",
                "file_driver" : "mpio",
                "echo_level" : 0
            }
            """)
        self.settings = settings.Clone()
        self.settings.ValidateAndAssignDefaults(default_settings)
        self.settings.AddEmptyValue("file_name")
    
    def Open(self, file_name):
        self.settings["file_name"].SetString(file_name)
        return KratosHDF5.HDF5FileParallel(self.settings)


class ModelPartOutput(OutputObject):
    """A class for writing a model part."""

    def __init__(self, settings):
        default_settings = KratosMultiphysics.Parameters("""
            {
                "prefix" : "/ModelData"
            }
            """)
        self.settings = settings.Clone()
        self.settings.ValidateAndAssignDefaults(default_settings)

    def Execute(self, model_part, hdf5_file):
        KratosHDF5.HDF5ModelPartIO(self.settings, hdf5_file).WriteModelPart(model_part)


class NodalResultsOutput(OutputObject):
    """A class for writing nodal results of a model part."""

    def __init__(self, settings):
        default_settings = KratosMultiphysics.Parameters("""
            {
                "prefix" : "/ResultsData",
                "list_of_variables": []
            }
            """)
        self.settings = settings.Clone()
        self.settings.ValidateAndAssignDefaults(default_settings)
        self.settings.AddEmptyValue("partitioned")
        self.settings["partitioned"].SetBool(False)

    def Execute(self, model_part, hdf5_file):
        KratosHDF5.HDF5NodalSolutionStepDataIO(self.settings, hdf5_file).WriteNodalResults(model_part.Nodes, 0)


class PartitionedModelPartOutput(OutputObject):
    """A class for writing a partitioned model part."""

    def __init__(self, settings):
        default_settings = KratosMultiphysics.Parameters("""
            {
                "prefix" : "/ModelData"
            }
            """)
        self.settings = settings.Clone()
        self.settings.ValidateAndAssignDefaults(default_settings)

    def Execute(self, model_part, hdf5_file):
        KratosHDF5.HDF5PartitionedModelPartIO(self.settings, hdf5_file).WriteModelPart(model_part)


class PartitionedNodalResultsOutput(OutputObject):
    """A class for writing nodal results of a partitioned model part."""

    def __init__(self, settings):
        default_settings = KratosMultiphysics.Parameters("""
            {
                "prefix" : "/ResultsData",
                "list_of_variables": []
            }
            """)
        self.settings = settings.Clone()
        self.settings.ValidateAndAssignDefaults(default_settings)
        self.settings.AddEmptyValue("partitioned")
        self.settings["partitioned"].SetBool(True)

    def Execute(self, model_part, hdf5_file):
        KratosHDF5.HDF5NodalSolutionStepDataIO(self.settings, hdf5_file).WriteNodalResults(model_part.Nodes, 0)


class TemporalOutputProcess(KratosMultiphysics.Process):
    """A process for writing temporal simulation results."""

    def __init__(self, model_part, hdf5_file_factory, settings):
        KratosMultiphysics.Process.__init__(self)
        default_settings = KratosMultiphysics.Parameters("""
            {
                "output_time_frequency": 1.0,
                "output_step_frequency": 0,
                "time_tag_precision": 4
            }
            """)
        settings.ValidateAndAssignDefaults(default_settings)
        self._model_part = model_part
        self._hdf5_file_factory = hdf5_file_factory
        self._output_time_frequency = settings["output_time_frequency"].GetDouble()
        self._output_step_frequency = settings["output_step_frequency"].GetInt()
        self._time_tag_precision = settings["time_tag_precision"].GetInt()
        self._list_of_outputs = []
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

class StaticOutputProcess(KratosMultiphysics.Process):
    """A process for writing static simulation results."""

    def __init__(self, model_part, hdf5_file_factory):
        KratosMultiphysics.Process.__init__(self)
        self._model_part = model_part
        self._hdf5_file_factory = hdf5_file_factory
        self._list_of_outputs = []

    def AddOutputObject(self, output):
        self._list_of_outputs.append(output)

    def Execute(self):
        hdf5_file = self._hdf5_file_factory.Open(self._model_part.Name + ".h5")
        for output in self._list_of_outputs:
            output.Execute(self._model_part, hdf5_file)

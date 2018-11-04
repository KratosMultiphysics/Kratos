"""A set of classes for performing HDF5 input/output.

BSD license: HDF5Application/license.txt
"""
from abc import ABCMeta, abstractmethod
import KratosMultiphysics
import KratosMultiphysics.HDF5Application as KratosHDF5
import hdf5_defaults


class FileFactory(metaclass=ABCMeta):

    @abstractmethod
    def Open(self, file_name): pass


class IOObject(metaclass=ABCMeta):

    @abstractmethod
    def Execute(self, model_part, hdf5_file): pass


class HDF5SerialFileFactory(FileFactory):

    def __init__(self, settings):
        default_settings = KratosMultiphysics.Parameters("""
            {
                "file_access_mode" : "exclusive",
                "file_driver" : "sec2",
                "echo_level" : 0
            }
            """)
        self.settings = settings.Clone()
        self.settings.ValidateAndAssignDefaults(default_settings)
        self.settings.AddEmptyValue("file_name")

    def Open(self, file_name):
        self.settings["file_name"].SetString(file_name)
        return KratosHDF5.HDF5FileSerial(self.settings)


class HDF5ParallelFileFactory(FileFactory):

    def __init__(self, settings):
        default_settings = KratosMultiphysics.Parameters("""
            {
                "file_access_mode" : "exclusive",
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


class ModelPartOutput(IOObject):
    """Provides the interface for writing a model part to a file."""

    def __init__(self, settings):
        default_settings = KratosMultiphysics.Parameters(hdf5_defaults.model_part_output_default_settings)

        self.settings = settings.Clone()
        self.settings.ValidateAndAssignDefaults(default_settings)

    def Execute(self, model_part, hdf5_file):
        KratosHDF5.HDF5ModelPartIO(hdf5_file, self.settings["prefix"].GetString()).WriteModelPart(model_part)

class ElementDataValueOutput(IOObject):
    """Provides the interface for writing element results to a file."""

    def __init__(self, settings):
        default_settings = KratosMultiphysics.Parameters(hdf5_defaults.temporal_default_settings)
        self.settings = settings.Clone()
        self.settings.ValidateAndAssignDefaults(default_settings)

    def Execute(self, model_part, hdf5_file):
        KratosHDF5.HDF5ElementDataValueIO(self.settings, hdf5_file).WriteElementResults(model_part.Elements)


class NodalSolutionStepDataOutput(IOObject):
    """Provides the interface for writing nodal solution step results to a file."""

    def __init__(self, settings):
        default_settings = KratosMultiphysics.Parameters(hdf5_defaults.temporal_default_settings)
        self.settings = settings.Clone()
        self.settings.ValidateAndAssignDefaults(default_settings)

    def Execute(self, model_part, hdf5_file):
        KratosHDF5.HDF5NodalSolutionStepDataIO(self.settings, hdf5_file).WriteNodalResults(model_part.Nodes, 0)

class NodalDataValueOutput(IOObject):
    """Provides the interface for writing nodal data values to a file."""

    def __init__(self, settings):
        default_settings = KratosMultiphysics.Parameters(hdf5_defaults.temporal_default_settings)
        self.settings = settings.Clone()
        self.settings.ValidateAndAssignDefaults(default_settings)

    def Execute(self, model_part, hdf5_file):
        KratosHDF5.HDF5NodalDataValueIO(self.settings, hdf5_file).WriteNodalResults(model_part.Nodes)

class PrimalBossakOutput(IOObject):
    """Provides the interface for writing a transient primal solution to a file."""

    def __init__(self, settings, alpha_bossak):
        default_settings = KratosMultiphysics.Parameters(hdf5_defaults.temporal_default_settings)
        self.settings = settings.Clone()
        self.settings.ValidateAndAssignDefaults(default_settings)
        self.alpha_bossak = alpha_bossak

    def Execute(self, model_part, hdf5_file):
        primal_io = KratosHDF5.HDF5NodalSolutionStepBossakIO(self.settings, hdf5_file)
        primal_io.SetAlphaBossak(self.alpha_bossak)
        primal_io.WriteNodalResults(model_part.Nodes)


class PrimalBossakInput(IOObject):
    """Provides the interface for reading a transient primal solution from a file."""

    def __init__(self, settings):
        default_settings = KratosMultiphysics.Parameters(hdf5_defaults.temporal_default_settings)
        self.settings = settings.Clone()
        self.settings.ValidateAndAssignDefaults(default_settings)

    def Execute(self, model_part, hdf5_file):
        primal_io = KratosHDF5.HDF5NodalSolutionStepBossakIO(self.settings, hdf5_file)
        primal_io.ReadNodalResults(model_part.Nodes, model_part.GetCommunicator())

class NodalSolutionStepDataInput(IOObject):
    """Provides the interface for reading a transient nodal solution step data from a file."""

    def __init__(self, settings):
        default_settings = KratosMultiphysics.Parameters(hdf5_defaults.temporal_default_settings)
        self.settings = settings.Clone()
        self.settings.ValidateAndAssignDefaults(default_settings)

    def Execute(self, model_part, hdf5_file):
        KratosHDF5.HDF5NodalSolutionStepDataIO(self.settings, hdf5_file).ReadNodalResults(model_part.Nodes, model_part.GetCommunicator(), 0)

class NodalDataValueInput(IOObject):
    """Provides the interface for reading a nodal data values from a file."""

    def __init__(self, settings):
        default_settings = KratosMultiphysics.Parameters(hdf5_defaults.temporal_default_settings)
        self.settings = settings.Clone()
        self.settings.ValidateAndAssignDefaults(default_settings)

    def Execute(self, model_part, hdf5_file):
        primal_io = KratosHDF5.HDF5NodalDataValueIO(self.settings, hdf5_file)
        primal_io.ReadNodalResults(model_part.Nodes, model_part.GetCommunicator())

class ElementDataValueInput(IOObject):
    """Provides the interface for reading element data values from a file."""

    def __init__(self, settings):
        default_settings = KratosMultiphysics.Parameters(hdf5_defaults.temporal_default_settings)
        self.settings = settings.Clone()
        self.settings.ValidateAndAssignDefaults(default_settings)

    def Execute(self, model_part, hdf5_file):
        KratosHDF5.HDF5ElementDataValueIO(self.settings, hdf5_file).ReadElementResults(model_part.Elements)

class PartitionedModelPartOutput(IOObject):
    """Provides the interface for writing a partitioned model part to a file."""

    def __init__(self, settings):
        default_settings = KratosMultiphysics.Parameters(hdf5_defaults.model_part_output_default_settings)
        self.settings = settings.Clone()
        self.settings.ValidateAndAssignDefaults(default_settings)

    def Execute(self, model_part, hdf5_file):
        KratosHDF5.HDF5PartitionedModelPartIO(hdf5_file, self.settings["prefix"].GetString()).WriteModelPart(model_part)


class StaticOutputProcess(KratosMultiphysics.Process):
    """A process for writing static simulation results."""

    def __init__(self, model_part, hdf5_file_factory, file_name):
        KratosMultiphysics.Process.__init__(self)
        self._model_part = model_part
        self._hdf5_file_factory = hdf5_file_factory
        self._list_of_outputs = []
        self._file_name = file_name

    def AddOutput(self, output):
        self._list_of_outputs.append(output)

    def Execute(self):
        hdf5_file = self._hdf5_file_factory.Open(self._file_name + ".h5")
        for output in self._list_of_outputs:
            output.Execute(self._model_part, hdf5_file)


class TemporalOutputProcess(KratosMultiphysics.Process):
    """A process for writing temporal simulation results.

    Responsible for the output step control logic. Output objects to be executed
    at regular time intervals are attached using AddOutput().
    """

    def __init__(self, model_part, hdf5_file_factory, settings, list_of_initial_outputs=[]):
        KratosMultiphysics.Process.__init__(self)
        default_settings = KratosMultiphysics.Parameters("""
            {
                "output_time_frequency": 1.0,
                "output_step_frequency": 0,
                "time_tag_precision": 4,
                "file_name": "%s"
            }
            """ % model_part.Name)
        settings.ValidateAndAssignDefaults(default_settings)
        self._model_part = model_part
        self._hdf5_file_factory = hdf5_file_factory
        self._time_step_file_name = settings["file_name"].GetString()
        self._initial_output = StaticOutputProcess(model_part, hdf5_file_factory, self._time_step_file_name)
        for output in list_of_initial_outputs:
            self._initial_output.AddOutput(output)
        self._output_time_frequency = settings["output_time_frequency"].GetDouble()
        self._output_step_frequency = settings["output_step_frequency"].GetInt()
        self._time_tag_precision = settings["time_tag_precision"].GetInt()
        self._list_of_outputs = []
        self._output_time = 0.0
        self._output_step = 0

    def AddOutput(self, output):
        # assert isinstance(output, IOObject)
        self._list_of_outputs.append(output)

    def ExecuteBeforeSolutionLoop(self):
        self._output_time = 0.0
        self._output_step = 0
        self._initial_output.Execute()

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
        return self._time_step_file_name + time_tag + ".h5"


class TemporalInputProcess(KratosMultiphysics.Process):
    """A process for reading temporal simulation results."""

    def __init__(self, model_part, hdf5_file_factory, settings):
        KratosMultiphysics.Process.__init__(self)
        default_settings = KratosMultiphysics.Parameters("""
            {
                "time_tag_precision": 4,
                "file_name": "%s"
            }
            """ % model_part.Name)
        settings.ValidateAndAssignDefaults(default_settings)
        self._model_part = model_part
        self._hdf5_file_factory = hdf5_file_factory
        self._time_step_file_name = settings["file_name"].GetString()
        self._time_tag_precision = settings["time_tag_precision"].GetInt()
        self._list_of_inputs = []

    def AddInput(self, i):
        self._list_of_inputs.append(i)

    def ExecuteInitializeSolutionStep(self):
        hdf5_file = self._hdf5_file_factory.Open(self._get_current_file_name())
        for i in self._list_of_inputs:
            i.Execute(self._model_part, hdf5_file)

    def _get_current_file_name(self):
        fmt = "{:." + str(self._time_tag_precision) + "f}"
        time_tag = "-" + fmt.format(self._model_part.ProcessInfo[KratosMultiphysics.TIME])
        return self._time_step_file_name + time_tag + ".h5"

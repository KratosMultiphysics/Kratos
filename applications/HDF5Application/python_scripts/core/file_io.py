'''HDF5 file IO.

license: HDF5Application/license.txt
'''


import os
from pathlib import Path


import KratosMultiphysics
import KratosMultiphysics.HDF5Application as KratosHDF5
import KratosMultiphysics.kratos_utilities as kratos_utils
from .utils import ParametersWrapper


class _FileIO(object):

    def _FileSettings(self, model_part):
        settings = ParametersWrapper()
        settings['file_name'] = self.filename_getter.Get(model_part)
        settings['file_access_mode'] = self.file_access_mode
        settings['file_driver'] = self.file_driver
        settings['echo_level'] = self.echo_level
        return settings


class _HDF5SerialFileIO(_FileIO):

    def Get(self, model_part=None):
        # return KratosHDF5.HDF5File(self._FileSettings(model_part).Get())
        return KratosHDF5.HDF5File(KratosMultiphysics.ParallelEnvironment.GetDefaultDataCommunicator(), self._FileSettings(model_part).Get())


class _HDF5ParallelFileIO(_FileIO):

    def Get(self, model_part=None):
        # return KratosHDF5.HDF5File(self._FileSettings(model_part).Get())
        return KratosHDF5.HDF5File(KratosMultiphysics.ParallelEnvironment.GetDefaultDataCommunicator(), self._FileSettings(model_part).Get())


class _HDF5MockFileIO(_FileIO):

    class MockFile(object):

        def __init__(self, file_name):
            self.file_name = file_name

        def GetFileName(self): return self.file_name

    def Get(self, model_part=None):
        settings = self._FileSettings(model_part)
        file_name = settings['file_name']
        return self.MockFile(file_name)


def _SetDefaults(settings):
    settings.SetDefault('io_type', 'serial_hdf5_file_io')
    settings.SetDefault('file_name', 'kratos')
    if '<time>' in settings['file_name']:
        settings.SetDefault('time_format', '0.4f')
    settings.SetDefault('file_access_mode', 'exclusive')
    if 'parallel' in settings['io_type']:
        settings.SetDefault('file_driver', 'mpio')
    elif os.name == 'nt':
        settings.SetDefault('file_driver', 'windows')
    else:
        settings.SetDefault('file_driver', 'sec2')
    settings.SetDefault('echo_level', 0)


def _GetIO(io_type):
    if io_type == 'serial_hdf5_file_io':
        io = _HDF5SerialFileIO()
    elif io_type == 'parallel_hdf5_file_io':
        io = _HDF5ParallelFileIO()
    elif io_type == 'mock_hdf5_file_io':
        io = _HDF5MockFileIO()
    else:
        raise ValueError('"io_type" has invalid value "' + io_type + '"')
    return io


class _FilenameGetter(object):

    def __init__(self, settings):
        self.filename = settings['file_name']
        self.model_part = None

        if (self.filename.find("<rank>") != -1):
            raise Exception("Flag \"<rank>\" is not allowed to be used in HDF5 output file namings. Please remove it [ filename = \"" + self.filename + "\" ].")

        if (not self.filename.endswith(".h5")):
            self.filename += ".h5"

        self.max_files_to_keep = None
        if ("file_access_mode" in settings.keys() and "max_files_to_keep" in settings.keys()):
            if (settings["file_access_mode"] == "truncate" and settings["max_files_to_keep"] != "unlimited"):
                self.max_files_to_keep = int(settings["max_files_to_keep"])
                if (self.max_files_to_keep <= 0):
                    raise Exception("max_files_to_keep should be greater than zero.")

        if settings.Has('time_format'):
            self.time_format = settings['time_format']
        else:
            self.time_format = ''

    def Get(self, model_part):
        if (self.model_part != model_part):
            self.model_part = model_part
            self.file_name_data_collector = KratosMultiphysics.FileNameDataCollector(self.model_part, self.filename, {"<time>": self.time_format})

        new_file_name = self.file_name_data_collector.GetFileName()

        if (self.max_files_to_keep is not None):
            if (Path(new_file_name).parents[0].is_dir()):
                list_of_file_names = self.file_name_data_collector.GetSortedFileNamesList(["<time>"])
                if (len(list_of_file_names) >= self.max_files_to_keep):
                    # remove files from the second file in case the mesh is only written to the initial file
                    # then we need to always have that initial file.
                    for file_name in list_of_file_names[1:len(list_of_file_names) - self.max_files_to_keep + 2]:
                        kratos_utils.DeleteFileIfExisting(file_name)

        return new_file_name


class _FilenameGetterWithDirectoryInitialization(object):

    def __init__(self, settings, data_comm):
        self.filename_getter = _FilenameGetter(settings)
        self.data_comm = data_comm

    def Get(self, model_part=None):
        file_name = self.filename_getter.Get(model_part)
        self._InitializeDirectory(file_name)
        return file_name

    def _InitializeDirectory(self, file_name):
        dirname = Path(file_name).absolute().parent
        KratosMultiphysics.FilesystemExtensions.MPISafeCreateDirectories(str(dirname))

def Create(settings, data_comm):
    '''Return the IO object specified by the setting 'io_type'.

    Empty settings will contain default values after returning from the
    function call.
    '''
    _SetDefaults(settings)
    io = _GetIO(settings['io_type'])
    io.filename_getter = _FilenameGetterWithDirectoryInitialization(settings, data_comm)
    io.file_access_mode = settings['file_access_mode']
    io.file_driver = settings['file_driver']
    io.echo_level = settings['echo_level']
    return io

class OpenHDF5File(object):
    """@brief A context responsible for managing the lifetime of HDF5 files."""

    def __init__(self, file_parameters: KratosMultiphysics.Parameters, model_part: KratosMultiphysics.ModelPart):
        if isinstance(file_parameters, ParametersWrapper):
            file_parameters = file_parameters.Get()
        self.__parameters = file_parameters
        self.__model_part = model_part
        self.__file: KratosHDF5.HDF5File = None

    def __enter__(self) -> KratosHDF5.HDF5File:
        distributed = "parallel_hdf5_file_io" if self.__model_part.IsDistributed() else "serial_hdf5_file_io"
        if self.__parameters.Has("io_type"):
            self.__parameters["io_type"].SetString(distributed)
        else:
            self.__parameters.AddString("io_type", distributed)
        self.__file = Create(
            ParametersWrapper(self.__parameters),
            self.__model_part.GetCommunicator().GetDataCommunicator()
        ).Get(self.__model_part)
        return self.__file

    def __exit__(self, exit_type, exit_value, exit_traceback) -> None:
        self.__file.Close()
        self.__file = None

    @property
    def file(self) -> KratosHDF5.HDF5File:
        return self.__file

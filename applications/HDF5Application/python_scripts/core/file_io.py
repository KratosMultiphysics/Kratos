'''HDF5 file IO.

license: HDF5Application/license.txt
'''


import os


import KratosMultiphysics
import KratosMultiphysics.HDF5Application as KratosHDF5
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
        return KratosHDF5.HDF5FileSerial(self._FileSettings(model_part).Get())


class _HDF5ParallelFileIO(_FileIO):

    def Get(self, model_part=None):
        return KratosHDF5.HDF5FileParallel(self._FileSettings(model_part).Get())


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
        filename = settings['file_name']
        self.filename_parts = filename.split('<time>')
        if settings.Has('time_format'):
            self.time_format = settings['time_format']
        else:
            self.time_format = ''

    def Get(self, model_part=None):
        if hasattr(model_part, 'ProcessInfo'):
            time = model_part.ProcessInfo[KratosMultiphysics.TIME]
            filename = format(time, self.time_format).join(self.filename_parts)
        else:
            filename = ''.join(self.filename_parts)
        if hasattr(model_part, 'Name'):
            filename = filename.replace('<model_part_name>', model_part.Name)
        if not filename.endswith('.h5'):
            filename += '.h5'
        return filename


class _FilenameGetterWithDirectoryInitialization(object):

    def __init__(self, settings):
        self.filename_getter = _FilenameGetter(settings)

    def Get(self, model_part=None):
        file_name = self.filename_getter.Get(model_part)
        self._InitializeDirectory(file_name)
        return file_name

    @staticmethod
    def _InitializeDirectory(file_name):
        dirname = os.path.dirname(file_name)
        if dirname != '' and not os.path.exists(dirname):
            if KratosMultiphysics.DataCommunicator.GetDefault().Rank() == 0:
                os.makedirs(dirname)
            KratosMultiphysics.DataCommunicator.GetDefault().Barrier()


def Create(settings):
    '''Return the IO object specified by the setting 'io_type'.

    Empty settings will contain default values after returning from the
    function call.
    '''
    _SetDefaults(settings)
    io = _GetIO(settings['io_type'])
    io.filename_getter = _FilenameGetterWithDirectoryInitialization(settings)
    io.file_access_mode = settings['file_access_mode']
    io.file_driver = settings['file_driver']
    io.echo_level = settings['echo_level']
    return io

'''HDF5 model part operations.

license: HDF5Application/license.txt
'''
import KratosMultiphysics
import KratosMultiphysics.HDF5Application as KratosHDF5
from .. import utils as _utils
from importlib import import_module

def _Prefix(pattern, identifier, time_format=''):
    if hasattr(identifier, 'ProcessInfo'):
        time = identifier.ProcessInfo[KratosMultiphysics.TIME]
        prefix = format(time, time_format).join(pattern.split('<time>'))
    else:
        prefix = pattern
    if hasattr(identifier, 'Name'):
        prefix = prefix.replace('<identifier>', identifier.Name)
    return prefix


class ModelPartOutput(object):
    '''Writes a model part to a file.'''

    def __init__(self, settings):
        default_setter = _utils.DefaultSetter(settings)
        default_setter.AddString('prefix', '/ModelData')
        self.prefix = settings['prefix'].GetString()
        if '<time>' in self.prefix:
            default_setter.AddString('time_format', '0.4f')
            self.time_format = settings['time_format'].GetString()

    def __call__(self, model_part, hdf5_file):
        if hasattr(self, 'time_format'):
            prefix = _Prefix(self.prefix, model_part, self.time_format)
        else:
            prefix = _Prefix(self.prefix, model_part)
        KratosHDF5.HDF5ModelPartIO(
            hdf5_file, prefix).WriteModelPart(model_part)


class PartitionedModelPartOutput(object):
    '''Writes a partitioned model part to a file.'''

    def __init__(self, settings):
        default_setter = _utils.DefaultSetter(settings)
        default_setter.AddString('prefix', '/ModelData')
        self.prefix = settings['prefix'].GetString()
        if '<time>' in self.prefix:
            default_setter.AddString('time_format', '0.4f')
            self.time_format = settings['time_format'].GetString()

    def __call__(self, model_part, hdf5_file):
        if hasattr(self, 'time_format'):
            prefix = _Prefix(self.prefix, model_part, self.time_format)
        else:
            prefix = _Prefix(self.prefix, model_part)
        KratosHDF5.HDF5PartitionedModelPartIO(
            hdf5_file, prefix).WriteModelPart(model_part)


class VariableIO(object):
    '''Generates json settings for variable data IO.'''

    def __init__(self, settings):
        default_setter = _utils.DefaultSetter(settings)
        default_setter.AddString('prefix', '/ResultsData')
        default_setter.AddArray('list_of_variables', [])
        self.prefix = settings['prefix'].GetString()
        if '<time>' in self.prefix:
            default_setter.AddString('time_format', '0.4f')
            self.time_format = settings['time_format'].GetString()
        self.list_of_variables = settings['list_of_variables']

    def GetSettings(self, model_part):
        settings = KratosMultiphysics.Parameters()
        if hasattr(self, 'time_format'):
            prefix = _Prefix(self.prefix, model_part, self.time_format)
        else:
            prefix = _Prefix(self.prefix, model_part)
        settings.AddEmptyValue('prefix').SetString(prefix)
        settings.AddValue('list_of_variables', self.list_of_variables)
        return settings


class ElementDataValueOutput(VariableIO):
    '''Writes non-historical element data values to a file.'''

    def __init__(self, settings):
        super(ElementDataValueOutput, self).__init__(settings)

    def __call__(self, model_part, hdf5_file):
        KratosHDF5.HDF5ElementDataValueIO(
            self.GetSettings(model_part), hdf5_file).WriteElementResults(model_part.Elements)


class ElementDataValueInput(VariableIO):
    '''Reads non-historical element data values from a file.'''

    def __init__(self, settings):
        super(ElementDataValueInput, self).__init__(settings)

    def __call__(self, model_part, hdf5_file):
        KratosHDF5.HDF5ElementDataValueIO(
            self.GetSettings(model_part), hdf5_file).ReadElementResults(model_part.Elements)


class NodalSolutionStepDataOutput(VariableIO):
    '''Writes nodal solution step data to a file.'''

    def __init__(self, settings):
        super(NodalSolutionStepDataOutput, self).__init__(settings)

    def __call__(self, model_part, hdf5_file):
        KratosHDF5.HDF5NodalSolutionStepDataIO(
            self.GetSettings(model_part), hdf5_file).WriteNodalResults(model_part.Nodes, 0)


class NodalSolutionStepDataInput(VariableIO):
    '''Reads nodal solution step data from a file.'''

    def __init__(self, settings):
        super(NodalSolutionStepDataInput, self).__init__(settings)

    def __call__(self, model_part, hdf5_file):
        nodal_io = KratosHDF5.HDF5NodalSolutionStepDataIO(
            self.GetSettings(model_part), hdf5_file)
        nodal_io.ReadNodalResults(
            model_part.Nodes, model_part.GetCommunicator(), 0)


class NodalDataValueOutput(VariableIO):
    '''Writes non-historical nodal data values to a file.'''

    def __init__(self, settings):
        super(NodalDataValueOutput, self).__init__(settings)

    def __call__(self, model_part, hdf5_file):
        KratosHDF5.HDF5NodalDataValueIO(
            self.GetSettings(model_part), hdf5_file).WriteNodalResults(model_part.Nodes)


class NodalDataValueInput(VariableIO):
    '''Reads non-historical nodal data values from a file.'''

    def __init__(self, settings):
        super(NodalDataValueInput, self).__init__(settings)

    def __call__(self, model_part, hdf5_file):
        primal_io = KratosHDF5.HDF5NodalDataValueIO(
            self.GetSettings(model_part), hdf5_file)
        primal_io.ReadNodalResults(
            model_part.Nodes, model_part.GetCommunicator())


class PrimalBossakOutput(VariableIO):
    '''Writes nodal solution step data to a file for Bossak time schemes.

    Behaves the same as NodalSolutionStepDataOutput except for ACCELERATION,
    which is computed as
    (1 - alpha_bossak) * node.GetSolutionStepValue(ACCELERATION, 0) +
          alpha_bossak * node.GetSolutionStepValue(ACCELERATION, 1)
    and written to the file for the current time step. This is used by the
    transient adjoint solvers.
    '''

    def __init__(self, settings):
        super(PrimalBossakOutput, self).__init__(settings)
        default_setter = _utils.DefaultSetter(settings)
        default_setter.AddDouble('alpha_bossak', -0.3)
        self.alpha_bossak = settings['alpha_bossak'].GetDouble()

    def __call__(self, model_part, hdf5_file):
        primal_io = KratosHDF5.HDF5NodalSolutionStepBossakIO(
            self.GetSettings(model_part), hdf5_file)
        primal_io.SetAlphaBossak(self.alpha_bossak)
        primal_io.WriteNodalResults(model_part.Nodes)


class PrimalBossakInput(VariableIO):
    '''Reads nodal solution step data from a file.

    This is used by the transient adjoint solvers.
    '''

    def __init__(self, settings):
        super(PrimalBossakInput, self).__init__(settings)

    def __call__(self, model_part, hdf5_file):
        primal_io = KratosHDF5.HDF5NodalSolutionStepBossakIO(
            self.GetSettings(model_part), hdf5_file)
        primal_io.ReadNodalResults(
            model_part.Nodes, model_part.GetCommunicator())


class MoveMesh(object):
    '''Perform a mesh move operation on a model part.

    The primary use case is to set the mesh to the current configuration after
    reading the model part.
    '''

    def __call__(self, model_part, *args):
        KratosMultiphysics.SolvingStrategy(model_part, True).MoveMesh()


def Create(settings):
    '''Return the operation specified by the setting 'operation_type'.

    If the 'operation_type' is not found and the settings have a 'module_name',
    the module is imported and used to create the operation. If 'module_name'
    is not found, an exception is raised. Empty settings will contain default
    values after returning from the function call.
    '''
    default_setter = _utils.DefaultSetter(settings)
    default_setter.AddString('operation_type', 'model_part_output')
    operation_type = settings['operation_type'].GetString()
    if operation_type == 'model_part_output':
        return ModelPartOutput(settings)
    elif operation_type == 'partitioned_model_part_output':
        return PartitionedModelPartOutput(settings)
    elif operation_type == 'element_data_value_output':
        return ElementDataValueOutput(settings)
    elif operation_type == 'element_data_value_input':
        return ElementDataValueInput(settings)
    elif operation_type == 'nodal_solution_step_data_output':
        return NodalSolutionStepDataOutput(settings)
    elif operation_type == 'nodal_solution_step_data_input':
        return NodalSolutionStepDataInput(settings)
    elif operation_type == 'nodal_data_value_output':
        return NodalDataValueOutput(settings)
    elif operation_type == 'nodal_data_value_input':
        return NodalDataValueInput(settings)
    elif operation_type == 'primal_bossak_output':
        return PrimalBossakOutput(settings)
    elif operation_type == 'primal_bossak_input':
        return PrimalBossakInput(settings)
    elif operation_type == 'move_mesh':
        return MoveMesh()
    else:
        if settings.Has('module_name'):
            module_name = settings['module_name'].GetString()
            module = import_module('KratosMultiphysics.HDF5Application.core.' + module_name)
            return module.Create(settings)
        raise ValueError(
            '"operation_type" has invalid value "' + operation_type + '"')

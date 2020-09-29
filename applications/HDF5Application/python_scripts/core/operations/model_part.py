'''HDF5 model part operations.

license: HDF5Application/license.txt
'''


from importlib import import_module


import KratosMultiphysics
import KratosMultiphysics.HDF5Application as KratosHDF5
from ..utils import ParametersWrapper


def Prefix(pattern, model_part, time_format=''):
    if hasattr(model_part, 'ProcessInfo'):
        time = model_part.ProcessInfo[KratosMultiphysics.TIME]
        prefix = format(time, time_format).join(pattern.split('<time>'))
    else:
        prefix = pattern
    if hasattr(model_part, 'Name'):
        prefix = prefix.replace('<model_part_name>', model_part.Name)
    return prefix


class ModelPartOutput:
    '''Writes a model part to a file.'''

    def __init__(self, settings):
        settings.SetDefault('prefix', '/ModelData')
        self.prefix = settings['prefix']
        if '<time>' in self.prefix:
            settings.SetDefault('time_format', '0.4f')
            self.time_format = settings['time_format']

    def __call__(self, model_part, hdf5_file):
        if hasattr(self, 'time_format'):
            prefix = Prefix(self.prefix, model_part, self.time_format)
        else:
            prefix = Prefix(self.prefix, model_part)
        KratosHDF5.HDF5ModelPartIO(
            hdf5_file, prefix).WriteModelPart(model_part)


class PartitionedModelPartOutput:
    '''Writes a partitioned model part to a file.'''

    def __init__(self, settings):
        settings.SetDefault('prefix', '/ModelData')
        self.prefix = settings['prefix']
        if '<time>' in self.prefix:
            settings.SetDefault('time_format', '0.4f')
            self.time_format = settings['time_format']

    def __call__(self, model_part, hdf5_file):
        if hasattr(self, 'time_format'):
            prefix = Prefix(self.prefix, model_part, self.time_format)
        else:
            prefix = Prefix(self.prefix, model_part)
        KratosHDF5.HDF5PartitionedModelPartIO(
            hdf5_file, prefix).WriteModelPart(model_part)


class VariableIO:
    '''Generates json settings for variable data IO.'''

    def __init__(self, settings):
        settings.SetDefault('prefix', '/ResultsData')
        settings.SetDefault('list_of_variables', [])
        self.prefix = settings['prefix']
        if '<time>' in self.prefix:
            settings.SetDefault('time_format', '0.4f')
            self.time_format = settings['time_format']
        self.list_of_variables = settings['list_of_variables']

    def GetSettings(self, model_part):
        settings = ParametersWrapper()
        if hasattr(self, 'time_format'):
            prefix = Prefix(self.prefix, model_part, self.time_format)
        else:
            prefix = Prefix(self.prefix, model_part)
        settings['prefix'] = prefix
        settings['list_of_variables'] = self.list_of_variables
        return settings


class ElementDataValueOutput(VariableIO):
    '''Writes non-historical element data values to a file.'''

    def __init__(self, settings):
        super(ElementDataValueOutput, self).__init__(settings)

    def __call__(self, model_part, hdf5_file):
        KratosHDF5.HDF5ElementDataValueIO(
            self.GetSettings(model_part).Get(), hdf5_file).WriteElementResults(model_part.Elements)


class ElementDataValueInput(VariableIO):
    '''Reads non-historical element data values from a file.'''

    def __init__(self, settings):
        super(ElementDataValueInput, self).__init__(settings)

    def __call__(self, model_part, hdf5_file):
        KratosHDF5.HDF5ElementDataValueIO(
            self.GetSettings(model_part).Get(), hdf5_file).ReadElementResults(model_part.Elements,
                                                                              model_part.GetCommunicator())

class ElementFlagValueOutput(VariableIO):
    '''Writes non-historical element flag values to a file.'''

    def __init__(self, settings):
        super(ElementFlagValueOutput, self).__init__(settings)

    def __call__(self, model_part, hdf5_file):
        KratosHDF5.HDF5ElementFlagValueIO(
            self.GetSettings(model_part).Get(), hdf5_file).WriteElementFlags(model_part.Elements)


class ElementFlagValueInput(VariableIO):
    '''Reads non-historical element flag values from a file.'''

    def __init__(self, settings):
        super(ElementFlagValueInput, self).__init__(settings)

    def __call__(self, model_part, hdf5_file):
        KratosHDF5.HDF5ElementFlagValueIO(
            self.GetSettings(model_part).Get(), hdf5_file).ReadElementFlags(model_part.Elements,
                                                                            model_part.GetCommunicator())

class ConditionDataValueOutput(VariableIO):
    '''Writes non-historical element data values to a file.'''

    def __init__(self, settings):
        super(ConditionDataValueOutput, self).__init__(settings)

    def __call__(self, model_part, hdf5_file):
        KratosHDF5.HDF5ConditionDataValueIO(
            self.GetSettings(model_part).Get(), hdf5_file).WriteConditionResults(model_part.Conditions)


class ConditionDataValueInput(VariableIO):
    '''Reads non-historical element data values from a file.'''

    def __init__(self, settings):
        super(ConditionDataValueInput, self).__init__(settings)

    def __call__(self, model_part, hdf5_file):
        KratosHDF5.HDF5ConditionDataValueIO(
            self.GetSettings(model_part).Get(), hdf5_file).ReadConditionResults(model_part.Conditions,
                                                                                model_part.GetCommunicator())

class ConditionFlagValueOutput(VariableIO):
    '''Writes non-historical element flag values to a file.'''

    def __init__(self, settings):
        super(ConditionFlagValueOutput, self).__init__(settings)

    def __call__(self, model_part, hdf5_file):
        KratosHDF5.HDF5ConditionFlagValueIO(
            self.GetSettings(model_part).Get(), hdf5_file).WriteConditionFlags(model_part.Conditions)


class ConditionFlagValueInput(VariableIO):
    '''Reads non-historical element flag values from a file.'''

    def __init__(self, settings):
        super(ConditionFlagValueInput, self).__init__(settings)

    def __call__(self, model_part, hdf5_file):
        KratosHDF5.HDF5ConditionFlagValueIO(
            self.GetSettings(model_part).Get(), hdf5_file).ReadConditionFlags(model_part.Conditions,
                                                                              model_part.GetCommunicator())


class NodalSolutionStepDataOutput(VariableIO):
    '''Writes nodal solution step data to a file.'''

    def __init__(self, settings):
        super(NodalSolutionStepDataOutput, self).__init__(settings)

    def __call__(self, model_part, hdf5_file):
        KratosHDF5.HDF5NodalSolutionStepDataIO(
            self.GetSettings(model_part).Get(), hdf5_file).WriteNodalResults(model_part.Nodes, 0)


class NodalSolutionStepDataInput(VariableIO):
    '''Reads nodal solution step data from a file.'''

    def __init__(self, settings):
        super(NodalSolutionStepDataInput, self).__init__(settings)

    def __call__(self, model_part, hdf5_file):
        nodal_io = KratosHDF5.HDF5NodalSolutionStepDataIO(
            self.GetSettings(model_part).Get(), hdf5_file)
        nodal_io.ReadNodalResults(
            model_part.Nodes, model_part.GetCommunicator(), 0)


class NodalDataValueOutput(VariableIO):
    '''Writes non-historical nodal data values to a file.'''

    def __init__(self, settings):
        super(NodalDataValueOutput, self).__init__(settings)

    def __call__(self, model_part, hdf5_file):
        KratosHDF5.HDF5NodalDataValueIO(
            self.GetSettings(model_part).Get(), hdf5_file).WriteNodalResults(model_part.Nodes)


class NodalDataValueInput(VariableIO):
    '''Reads non-historical nodal data values from a file.'''

    def __init__(self, settings):
        super(NodalDataValueInput, self).__init__(settings)

    def __call__(self, model_part, hdf5_file):
        primal_io = KratosHDF5.HDF5NodalDataValueIO(
            self.GetSettings(model_part).Get(), hdf5_file)
        primal_io.ReadNodalResults(
            model_part.Nodes, model_part.GetCommunicator())

class NodalFlagValueOutput(VariableIO):
    '''Writes non-historical nodal flag values to a file.'''

    def __init__(self, settings):
        super(NodalFlagValueOutput, self).__init__(settings)

    def __call__(self, model_part, hdf5_file):
        KratosHDF5.HDF5NodalFlagValueIO(
            self.GetSettings(model_part).Get(), hdf5_file).WriteNodalFlags(model_part.Nodes)


class NodalFlagValueInput(VariableIO):
    '''Reads non-historical nodal flag values from a file.'''

    def __init__(self, settings):
        super(NodalFlagValueInput, self).__init__(settings)

    def __call__(self, model_part, hdf5_file):
        primal_io = KratosHDF5.HDF5NodalFlagValueIO(
            self.GetSettings(model_part).Get(), hdf5_file)
        primal_io.ReadNodalFlags(
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
        settings.SetDefault('alpha_bossak', -0.3)
        self.alpha_bossak = settings['alpha_bossak']

    def __call__(self, model_part, hdf5_file):
        primal_io = KratosHDF5.HDF5NodalSolutionStepBossakIO(
            self.GetSettings(model_part).Get(), hdf5_file)
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
            self.GetSettings(model_part).Get(), hdf5_file)
        primal_io.ReadNodalResults(
            model_part.Nodes, model_part.GetCommunicator())


class MoveMesh:
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
    settings.SetDefault('operation_type', 'model_part_output')
    operation_type = settings['operation_type']
    if operation_type == 'model_part_output':
        return ModelPartOutput(settings)
    elif operation_type == 'partitioned_model_part_output':
        return PartitionedModelPartOutput(settings)
    elif operation_type == 'element_data_value_output':
        return ElementDataValueOutput(settings)
    elif operation_type == 'element_flag_value_output':
        return ElementFlagValueOutput(settings)
    elif operation_type == 'element_data_value_input':
        return ElementDataValueInput(settings)
    elif operation_type == 'element_flag_value_input':
        return ElementFlagValueInput(settings)
    elif operation_type == 'condition_data_value_output':
        return ConditionDataValueOutput(settings)
    elif operation_type == 'condition_flag_value_output':
        return ConditionFlagValueOutput(settings)
    elif operation_type == 'condition_data_value_input':
        return ConditionDataValueInput(settings)
    elif operation_type == 'condition_flag_value_input':
        return ConditionFlagValueInput(settings)
    elif operation_type == 'nodal_solution_step_data_output':
        return NodalSolutionStepDataOutput(settings)
    elif operation_type == 'nodal_solution_step_data_input':
        return NodalSolutionStepDataInput(settings)
    elif operation_type == 'nodal_data_value_output':
        return NodalDataValueOutput(settings)
    elif operation_type == 'nodal_flag_value_output':
        return NodalFlagValueOutput(settings)
    elif operation_type == 'nodal_data_value_input':
        return NodalDataValueInput(settings)
    elif operation_type == 'nodal_flag_value_input':
        return NodalFlagValueInput(settings)
    elif operation_type == 'primal_bossak_output':
        return PrimalBossakOutput(settings)
    elif operation_type == 'primal_bossak_input':
        return PrimalBossakInput(settings)
    elif operation_type == 'move_mesh':
        return MoveMesh()
    else:
        if settings.Has('module_name'):
            module_name = settings['module_name']
            module = import_module(
                'KratosMultiphysics.HDF5Application.core.' + module_name)
            return module.Create(settings)
        raise ValueError(
            '"operation_type" has invalid value "' + operation_type + '"')

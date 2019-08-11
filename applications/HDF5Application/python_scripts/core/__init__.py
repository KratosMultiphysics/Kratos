'''Core HDF5 IO.

license: HDF5Application/license.txt
'''


__all__ = ["Factory"]


import KratosMultiphysics
from . import processes, controllers, operations, file_io, utils
from .utils import ParametersWrapper


def CreateControllerWithFileIO(settings, model):
    default_setter = utils.DefaultSetter(settings)
    default_setter.AddString(
        'model_part_name', 'PLEASE_SPECIFY_MODEL_PART_NAME')
    default_setter.AddString('process_step', 'initialize')
    default_setter.Add('controller_settings')
    default_setter.Add('io_settings')
    default_setter.AddArray('list_of_operations', [])
    if settings['list_of_operations'].size() == 0:
        settings['list_of_operations'].Append(
            KratosMultiphysics.Parameters())
    model_part = model[settings['model_part_name'].GetString()]
    return controllers.Create(model_part, file_io.Create(settings['io_settings']), settings['controller_settings'])


def AssignOperationsToController(operations_settings, controller):
    if not operations_settings.IsArray():
        raise ValueError('Expected settings as an array')
    for settings in operations_settings:
        controller.Add(operations.Create(settings))


def AssignControllerToProcess(settings, controller, process):
    process_step = settings['process_step'].GetString()
    if process_step == 'initialize':
        process.AddInitialize(controller)
    elif process_step == 'before_solution_loop':
        process.AddBeforeSolutionLoop(controller)
    elif process_step == 'initialize_solution_step':
        process.AddInitializeSolutionStep(controller)
    elif process_step == 'finalize_solution_step':
        process.AddFinalizeSolutionStep(controller)
    elif process_step == 'before_output_step':
        process.AddBeforeOutputStep(controller)
    elif process_step == 'after_output_step':
        process.AddAfterOutputStep(controller)
    elif process_step == 'finalize':
        process.AddFinalize(controller)
    else:
        raise ValueError(
            '"process_step" has invalid value "' + process_step + '"')


def Factory(settings, model):
    '''Return an HDF5 IO process specified by json settings.'''
    if not settings.IsArray():
        raise ValueError('Expected settings as an array')
    if settings.size() == 0:
        settings.Append(KratosMultiphysics.Parameters())
    process = processes.ControllerProcess()
    for current_settings in settings:
        controller = CreateControllerWithFileIO(current_settings, model)
        AssignOperationsToController(
            current_settings['list_of_operations'], controller)
        AssignControllerToProcess(current_settings, controller, process)
    return process

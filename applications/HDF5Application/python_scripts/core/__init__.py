'''Core HDF5 IO.

license: HDF5Application/license.txt
'''


__all__ = ["Factory"]


import KratosMultiphysics
from . import processes
from . import controllers
from . import operations
from . import file_io
from .utils import ParametersWrapper


def CreateControllerWithFileIO(settings, model):
    settings.SetDefault('model_part_name', 'PLEASE_SPECIFY_MODEL_PART_NAME')
    settings.SetDefault('process_step', 'initialize')
    settings.SetDefault('controller_settings')
    settings.SetDefault('io_settings')
    settings.SetDefault('list_of_operations', [])
    if len(settings['list_of_operations']) == 0:
        settings['list_of_operations'].Append(KratosMultiphysics.Parameters())
    model_part = model[settings['model_part_name']]
    data_comm = model_part.GetCommunicator().GetDataCommunicator()
    return controllers.Create(
        model_part, file_io.Create(settings['io_settings'], data_comm),
        settings['controller_settings'])


def AssignOperationsToController(settings, controller):
    if not settings.IsArray():
        raise ValueError('Expected settings as an array')
    for i in settings:
        controller.Add(operations.Create(settings[i]))


def AssignControllerToProcess(settings, controller, process):
    process_step = settings['process_step']
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
    if len(settings) == 0:
        settings.Append(KratosMultiphysics.Parameters())
    process = processes.ControllerProcess()
    for i in settings:
        controller = CreateControllerWithFileIO(settings[i], model)
        AssignOperationsToController(
            settings[i]['list_of_operations'], controller)
        AssignControllerToProcess(settings[i], controller, process)
    return process

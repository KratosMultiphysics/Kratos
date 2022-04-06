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

from enum import Enum


class ProcessTag(Enum):
    UNDEFINED = 1
    INPUT = 2
    OUTPUT = 3


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


def Factory(settings: KratosMultiphysics.Parameters, model: KratosMultiphysics.Model, tag: ProcessTag):
    '''Return an HDF5 IO process specified by json settings.'''
    if not settings.IsArray():
        raise ValueError('Expected settings as an array')
    if len(settings) == 0:
        settings.Append(KratosMultiphysics.Parameters())

    if tag == ProcessTag.INPUT:
        process_base = KratosMultiphysics.Process
    elif tag == ProcessTag.OUTPUT:
        process_base = KratosMultiphysics.OutputProcess
    elif tag == ProcessTag.UNDEFINED:
        process_base = KratosMultiphysics.Process
    process = processes.Factory(process_base)

    for i in settings:
        controller = CreateControllerWithFileIO(settings[i], model)
        AssignOperationsToController(
            settings[i]['list_of_operations'], controller)
        controller.AssignToProcess(process, settings[i])
    return process

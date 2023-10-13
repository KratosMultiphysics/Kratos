"""!@package HDF5Application

Core HDF5 IO.

license: HDF5Application/license.txt
"""


__all__ = ["Factory"]


# --- Core Imports ---
import KratosMultiphysics

# --- HDF5 Imports ---
from . import processes
from . import controllers
from . import operations
from .utils import ParametersWrapper

# --- STD Imports ---
import typing


##!@addtogroup HDF5Application
##!@{
def CreateController(parameters: KratosMultiphysics.Parameters,
                     model: KratosMultiphysics.Model,
                     operation: operations.AggregateOperation) -> controllers.Controller:
    parameters.AddMissingParameters(KratosMultiphysics.Parameters("""{
        "model_part_name" : "PLEASE_SPECIFY_MODEL_PART_NAME",
        "process_step" : "initialize",
        "controller_settings" : {},
        "io_settings" : {},
        "list_of_operations" : []
    }"""))
    model_part = model[parameters['model_part_name'].GetString()]
    return controllers.Factory(model_part, operation, parameters['controller_settings'])


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
    elif process_step == "output":
        # Processes assigned to the 'output' step must not have
        # TemporalControllers because they must write the output
        # on each request.
        if not isinstance(controller, controllers.DefaultController):
            raise TypeError("Processes assigned to 'output' must have a DefaultController that executes on each call. The specified controller instead is: {}".format(type(controller)))
        process.AddOutput(controller)
    else:
        raise ValueError(f'"process_step" has invalid value "{process_step}"')


def Factory(settings: ParametersWrapper, model: KratosMultiphysics.Model, process_base: type) -> typing.Union[KratosMultiphysics.Process, KratosMultiphysics.OutputProcess]:
    """Return an HDF5 IO process specified by json settings."""
    if not settings.IsArray():
        raise ValueError('Expected settings as an array')
    if len(settings) == 0:
        settings.Append(KratosMultiphysics.Parameters())
    process = processes.Factory(process_base)
    for i in settings:
        operation = operations.AggregateOperation(model, settings[i].Get())
        controller = CreateController(settings[i].Get(), model, operation)
        AssignControllerToProcess(settings[i], controller, process)
    return process
##!@}

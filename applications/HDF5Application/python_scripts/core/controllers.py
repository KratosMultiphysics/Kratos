"""!@package HDF5Application

HDF5 controllers.

This module contains the controllers which control the frequency of HDF5 IO
operations.

license: HDF5Application/license.txt
"""


# Kratos imports
import KratosMultiphysics
from . import operations

# STL imports
import abc


##!@addtogroup HDF5Application
##!@{
##!@name Kratos classes
##!@{
class Controller(metaclass=abc.ABCMeta):
    def __init__(self,
                 model_part: KratosMultiphysics.ModelPart,
                 operation: operations.AggregateOperation):
        self.model_part = model_part
        self.__operation = operation

    @abc.abstractmethod
    def IsExecuteStep(self) -> bool:
        """!Decide whether the controller should execute the registered operations at the current step."""
        pass

    @abc.abstractmethod
    def __call__(self) -> None:
        """!Execute assigned operations if a check is passed."""
        pass

    def ExecuteOperation(self) -> None:
        """!Execute all assigned operations, bypassing the controller's checks."""
        self.__operation.Execute()


class DefaultController(Controller):
    """Simple pass through controller."""

    def IsExecuteStep(self) -> bool:
        return True

    def __call__(self) -> None:
        self.ExecuteOperation()


class TemporalController(Controller):
    """!@brief Frequency-based controller.
    @detail Controls execution according to the 'time_frequency' and 'step_frequency'
    specified in the json settings.
    """

    def __init__(self,
                 model_part: KratosMultiphysics.ModelPart,
                 operation: operations.AggregateOperation,
                 settings: KratosMultiphysics.Parameters):
        super().__init__(model_part, operation)
        settings.AddMissingParameters(KratosMultiphysics.Parameters("""{
            "time_frequency" : 1.0,
            "step_frequency" : 1
        }"""))
        self.time_frequency = settings['time_frequency'].GetDouble()
        self.step_frequency = settings['step_frequency'].GetInt()
        self.current_time = 0.0
        self.current_step = 0

    def IsExecuteStep(self) -> bool:
        """!@brief Return true if the current step/time is a multiple of the output frequency.
        @detail Relative errors are compared against an epsilon, which is much larger than
        the machine epsilon, and include a lower bound based on
        https://github.com/chromium/chromium, cc::IsNearlyTheSame.
        """
        if self.current_step == self.step_frequency:
            return True
        if self.current_time > self.time_frequency:
            return True
        eps = 1e-6
        tol = eps * max(abs(self.current_time), abs(self.time_frequency), eps)
        if abs(self.current_time - self.time_frequency) < tol:
            return True
        return False

    def __call__(self) -> None:
        # TODO: separately keeping track of steps and time internally is not a good
        # idea. What happens if the solution process involves jumping back and forth
        # in time (restarts, checkpointing)? @matekelemen
        delta_time = self.model_part.ProcessInfo[KratosMultiphysics.DELTA_TIME]
        self.current_time += delta_time
        self.current_step += 1
        if self.IsExecuteStep():
            self.ExecuteOperation()
            self.current_time = 0.0
            self.current_step = 0
##!@}


def Factory(model_part: KratosMultiphysics.ModelPart,
            operation: operations.AggregateOperation,
            parameters: KratosMultiphysics.Parameters) -> Controller:
    """!@brief Return the controller specified by the setting 'controller_type'.
    @detail Empty settings will contain default values after returning from the
    function call.
    """
    parameters.AddMissingParameters(KratosMultiphysics.Parameters("""{
        "controller_type" : "default_controller"
    }"""))
    controller_type = parameters['controller_type'].GetString()
    if controller_type == 'default_controller':
        return DefaultController(model_part, operation)
    elif controller_type == 'temporal_controller':
        return TemporalController(model_part, operation, parameters)
    else:
        raise ValueError(f'"controller_type" has invalid value "{controller_type}"')
##!@}

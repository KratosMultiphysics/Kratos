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
from typing import Optional


##!@addtogroup HDF5Application
##!@{
##!@name Kratos classes
##!@{
class Controller(abc.ABC):
    def __init__(self,
                 model_part: KratosMultiphysics.ModelPart,
                 operation: operations.AggregateOperation):
        self._model_part: KratosMultiphysics.ModelPart = model_part
        self.__operation: operations.AggregateOperation = operation

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
        time_frequency: Optional[float] = None
        step_frequency: Optional[int] = None
        if settings.Has("time_frequency"):
            time_frequency = settings["time_frequency"].GetDouble()
        if settings.Has("step_frequency"):
            step_frequency = settings["step_frequency"].GetInt()

        # Neither step nor time frequency were defined => apply defaults
        # => output at every time step (step_frequency is 1, time_frequency is undefined)
        if time_frequency is None and step_frequency is None:
            settings.AddInt("step_frequency", 1)
            step_frequency = 1

        # Time frequency was not defined => the output won't be triggered by changes in TIME
        if time_frequency is None:
            time_frequency = float("inf")

        self.__time_frequency: float = time_frequency
        self.__step_frequency: Optional[int] = step_frequency
        self.__last_output_time = self._model_part.ProcessInfo[KratosMultiphysics.TIME]
        self.__last_output_step = self._model_part.ProcessInfo[KratosMultiphysics.STEP]

    def ExecuteOperation(self) -> None:
        super().ExecuteOperation()
        self.__last_output_time = self._model_part.ProcessInfo[KratosMultiphysics.TIME]
        self.__last_output_step = self._model_part.ProcessInfo[KratosMultiphysics.STEP]

    def IsExecuteStep(self) -> bool:
        """!@brief Return true if the current step/time is a multiple of the output frequency.
        @detail Relative errors are compared against an epsilon, which is much larger than
        the machine epsilon, and include a lower bound based on
        https://github.com/chromium/chromium, cc::IsNearlyTheSame.
        """
        # TODO: separately keeping track of steps and time internally is not a good
        # idea. What happens if the solution process involves jumping back and forth
        # in time (restarts, checkpointing)? @matekelemen
        if self.__step_frequency is not None:
            step_difference = self._model_part.ProcessInfo[KratosMultiphysics.STEP] - self.__last_output_step
            if self.__step_frequency <= step_difference:
                return True

        time_difference = self._model_part.ProcessInfo[KratosMultiphysics.TIME] - self.__last_output_time
        if self.__time_frequency <= time_difference:
            return True

        eps = 1e-6
        tol = eps * max(abs(time_difference), abs(self.__time_frequency), eps)
        if abs(time_difference - self.__time_frequency) < tol:
            return True
        return False

    def __call__(self) -> None:
        if self.IsExecuteStep():
            self.ExecuteOperation()
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

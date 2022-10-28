"""!@package HDF5Application

HDF5 controllers.

This module contains the controllers which control the frequency of HDF5 IO
operations.

license: HDF5Application/license.txt
"""


# Kratos imports
from tokenize import Single
import KratosMultiphysics
from . import file_io

# STL imports
import abc


##!@addtogroup HDF5Application
##!@{
##!@name Kratos classes
##!@{
class Controller(metaclass=abc.ABCMeta):
    def __init__(self, model_part: KratosMultiphysics.ModelPart, io: file_io._FileIO):
        self.model_part = model_part
        self.io = io
        self.operations = []

    def Add(self, operation) -> None:
        self.operations.append(operation)

    @abc.abstractmethod
    def IsExecuteStep(self) -> bool:
        """!Decide whether the controller should execute the registered operations at the current step."""
        pass

    @abc.abstractmethod
    def __call__(self) -> None:
        """!Execute assigned operations if a check is passed."""
        pass

    def ExecuteOperations(self) -> None:
        """!Execute all assigned operations, bypassing the controller's checks."""
        current_io = self.io.Get(self.model_part)
        for op in self.operations:
            op(self.model_part, current_io)


class DefaultController(Controller):
    """Simple pass through controller."""

    def IsExecuteStep(self) -> bool:
        return True

    def __call__(self) -> None:
        self.ExecuteOperations()


class TemporalController(Controller):
    """!@brief Frequency-based controller.
    @detail Controls execution according to the 'time_frequency' and 'step_frequency'
    specified in the json settings.
    """

    def __init__(self, model_part: KratosMultiphysics.ModelPart, io: file_io._FileIO, settings: KratosMultiphysics.Parameters):
        super().__init__(model_part, io)
        settings.SetDefault('time_frequency', 1.0)
        settings.SetDefault('step_frequency', 1)
        self.time_frequency = settings['time_frequency']
        self.step_frequency = settings['step_frequency']
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
            self.ExecuteOperations()
            self.current_time = 0.0
            self.current_step = 0
##!@}


def Factory(model_part: KratosMultiphysics.ModelPart, io: file_io._FileIO, settings: KratosMultiphysics.Parameters) -> Controller:
    """!@brief Return the controller specified by the setting 'controller_type'.
    @detail Empty settings will contain default values after returning from the
    function call.
    """
    settings.SetDefault('controller_type', 'default_controller')
    controller_type = settings['controller_type']
    if controller_type == 'default_controller':
        return DefaultController(model_part, io)
    elif controller_type == 'temporal_controller':
        return TemporalController(model_part, io, settings)
    else:
        raise ValueError(
            '"controller_type" has invalid value "' + controller_type + '"')
##!@}
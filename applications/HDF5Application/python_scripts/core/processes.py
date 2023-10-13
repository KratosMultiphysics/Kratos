"""!@package HDF5Application

HDF5 core processes.

This module only contains the most general HDF5 IO processes which should not
change frequently.

license: HDF5Application/license.txt
"""

__all__ = ["Factory"]

# Kratos imports
import KratosMultiphysics
from . import controllers

# STL imports
import inspect


##!@addtogroup HDF5Application
##!@{
##!@name Kratos classes
##!@{
class OrderedOperationProcess(KratosMultiphysics.Process):
    """!@brief A process for grouping operations.
    @detail This implements a whole-part structural decomposition. The members are
    operations or function objects with no arguments. They may be attached
    to any of the process steps during construction and are called in the same
    order at the corresponding step of the solution algorithm.
    """

    def __init__(self):
        KratosMultiphysics.Process.__init__(self)
        self.Clear()

    def ExecuteInitialize(self) -> None:
        for controller in self.initialize_controllers:
            controller()

    def ExecuteBeforeSolutionLoop(self) -> None:
        for controller in self.before_solution_loop_controllers:
            controller()

    def ExecuteInitializeSolutionStep(self) -> None:
        for controller in self.initialize_solution_step_controllers:
            controller()

    def ExecuteFinalizeSolutionStep(self) -> None:
        for controller in self.finalize_solution_step_controllers:
            controller()

    def ExecuteBeforeOutputStep(self) -> None:
        for controller in self.before_output_step_controllers:
            controller()

    def ExecuteAfterOutputStep(self) -> None:
        for controller in self.after_output_step_controllers:
            controller()

    def ExecuteFinalize(self) -> None:
        for controller in self.finalize_controllers:
            controller()

    def Clear(self) -> None:
        """!Remove all assigned controllers."""
        self.initialize_controllers: list[controllers.Controller] = []
        self.before_solution_loop_controllers: list[controllers.Controller] = []
        self.initialize_solution_step_controllers: list[controllers.Controller] = []
        self.finalize_solution_step_controllers: list[controllers.Controller] = []
        self.before_output_step_controllers: list[controllers.Controller] = []
        self.after_output_step_controllers: list[controllers.Controller] = []
        self.finalize_controllers: list[controllers.Controller] = []

    def AddInitialize(self, controller: controllers.Controller) -> None:
        """!Add a controller to be executed when ExecuteInitialize is called."""
        self.initialize_controllers.append(controller)

    def AddBeforeSolutionLoop(self, controller: controllers.Controller) -> None:
        """!Add a controller to be executed when ExecuteBeforeSolutionLoop is called."""
        self.before_solution_loop_controllers.append(controller)

    def AddInitializeSolutionStep(self, controller: controllers.Controller) -> None:
        """!Add a controller to be executed when ExecuteInitializeSolutionStep is called."""
        self.initialize_solution_step_controllers.append(controller)

    def AddFinalizeSolutionStep(self, controller: controllers.Controller) -> None:
        """!Add a controller to be executed when ExecuteFinalizeSolutionStep is called."""
        self.finalize_solution_step_controllers.append(controller)

    def AddBeforeOutputStep(self, controller: controllers.Controller) -> None:
        """!Add a controller to be executed when ExecuteBeforeOutputStep is called."""
        self.before_output_step_controllers.append(controller)

    def AddAfterOutputStep(self, controller: controllers.Controller) -> None:
        """!Add a controller to be executed when ExecuteAfterOutputStep is called."""
        self.after_output_step_controllers.append(controller)

    def AddFinalize(self, controller: controllers.Controller) -> None:
        """!Add a controller to be executed when ExecuteFinalize is called."""
        self.finalize_controllers.append(controller)


class OrderedOutputOperationProcess(KratosMultiphysics.OutputProcess, OrderedOperationProcess):

    def __init__(self):
        KratosMultiphysics.OutputProcess.__init__(self)
        OrderedOperationProcess.__init__(self)

    def IsOutputStep(self) -> bool:
        """!True if this step is an output step for any assigned controller."""
        return any(controller.IsExecuteStep() for controller in self.__GetLoopControllers())

    def PrintOutput(self) -> None:
        """!Bypass checks and execute all processes assigned to each registered controller."""
        for controller in self.__GetLoopControllers():
            controller.ExecuteOperation()

    # TODO: remove the double quotes from the return type hint after adopting python3.9 or newer
    def __GetLoopControllers(self) -> "list[controllers.Controller]":
        """!Get controllers registered to loop events."""
        return sum([self.initialize_solution_step_controllers,self.finalize_solution_step_controllers,self.before_output_step_controllers,self.after_output_step_controllers,self.finalize_controllers],[])
##!@}


def Factory(process_type: type) -> KratosMultiphysics.Process:
    """!@brief Return a regular or output ordered process depending on the input type.
    @param process_type: @ref Process or @ref OutputProcess
    @return an @ref OrdereOperationProcess or an @ref OrderedOutputOperationProcess depending on the input argument
    """
    if not inspect.isclass(process_type):
        raise TypeError("Expecting a type, but got {}".format(process_type))
    if issubclass(process_type, KratosMultiphysics.OutputProcess):
        return OrderedOutputOperationProcess()
    elif issubclass(process_type, KratosMultiphysics.Process):
        return OrderedOperationProcess()
    else:
        raise TypeError("Expecting OutputProcess or Process types, but got {}".format(process_type.__name__))
##!@}

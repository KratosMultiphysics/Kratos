'''HDF5 core processes.

This module only contains the most general HDF5 IO processes which should not
change frequently.

license: HDF5Application/license.txt
'''


__all__ = ["ControllerProcess"]


import KratosMultiphysics


class OrderedOperationProcess(KratosMultiphysics.Process):
    """A process for grouping operations.

    This implements a whole-part structural decomposition. The members are
    operations or function objects with no arguments. They may be attached
    to any of the process steps during construction and are called in the same
    order at the corresponding step of the solution algorithm.
    """

    def __init__(self):
        KratosMultiphysics.Process.__init__(self)
        self.initialize_sequence = []
        self.before_solution_loop_sequence = []
        self.initialize_solution_step_sequence = []
        self.finalize_solution_step_sequence = []
        self.before_output_step_sequence = []
        self.after_output_step_sequence = []
        self.finalize_sequence = []

    def ExecuteInitialize(self):
        for f in self.initialize_sequence:
            f()

    def ExecuteBeforeSolutionLoop(self):
        for f in self.before_solution_loop_sequence:
            f()

    def ExecuteInitializeSolutionStep(self):
        for f in self.initialize_solution_step_sequence:
            f()

    def ExecuteFinalizeSolutionStep(self):
        for f in self.finalize_solution_step_sequence:
            f()

    def ExecuteBeforeOutputStep(self):
        for f in self.before_output_step_sequence:
            f()

    def ExecuteAfterOutputStep(self):
        for f in self.after_output_step_sequence:
            f()

    def ExecuteFinalize(self):
        for f in self.finalize_sequence:
            f()

    def Clear(self):
        self.initialize_sequence = []
        self.before_solution_loop_sequence = []
        self.initialize_solution_step_sequence = []
        self.finalize_solution_step_sequence = []
        self.before_output_step_sequence = []
        self.after_output_step_sequence = []
        self.finalize_sequence = []

    def AddInitialize(self, func):
        self.initialize_sequence.append(func)

    def AddBeforeSolutionLoop(self, func):
        self.before_solution_loop_sequence.append(func)

    def AddInitializeSolutionStep(self, func):
        self.initialize_solution_step_sequence.append(func)

    def AddFinalizeSolutionStep(self, func):
        self.finalize_solution_step_sequence.append(func)

    def AddBeforeOutputStep(self, func):
        self.before_output_step_sequence.append(func)

    def AddAfterOutputStep(self, func):
        self.after_output_step_sequence.append(func)

    def AddFinalize(self, func):
        self.finalize_sequence.append(func)


class ControllerProcess(OrderedOperationProcess):
    """A process for grouping controllers.

    This adds a shorter more specific name for the process that groups
    controllers.
    """
    pass

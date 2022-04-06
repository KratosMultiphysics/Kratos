'''HDF5 core processes.

This module only contains the most general HDF5 IO processes which should not
change frequently.

license: HDF5Application/license.txt
'''


__all__ = ["Factory"]


import KratosMultiphysics


def Factory(base_type: type) -> KratosMultiphysics.Process:
    """"""
    # The base type can either be KratosMultiphysics.Process (for input processes)
    # or KratosMultiphysics.OutputProcess (for output processes). Anything else will
    # result in an exception.
    if base_type == KratosMultiphysics.Process:
        process_name = "InputControllerProcess"
    elif base_type == KratosMultiphysics.OutputProcess:
        process_name = "OutputControllerProcess"
    else:
        raise ValueError("Expecting KratosMultiphysics.Process or KratosMultiphysics.OutputProcess, but got {}".format(base_type))

    class ControllerProcess(base_type):
        """A process for grouping operations.

        This implements a whole-part structural decomposition. The members are
        operations or function objects with no arguments. They may be attached
        to any of the process steps during construction and are called in the same
        order at the corresponding step of the solution algorithm.
        """

        __name__ = process_name
        __qualname__ = process_name

        def __init__(self):
            super().__init__()
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

    return ControllerProcess()

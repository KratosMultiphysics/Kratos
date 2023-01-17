from inspect import getmro

import KratosMultiphysics as Kratos
from KratosMultiphysics.OptimizationApplication.optimization_routine import OptimizationRoutine

class OptimizationInfo:
    class __SolutionStepData(dict):
        """Solution step data holder dict

        This is a dictionary with override __getite__ to have
        proper error message.

        It also introduces HasSolutionStepDataKey which checks whether the
        given key is available in the solution step data.

        Args:
            dict (_type_): Solution step data holder dictionary.
        """
        def __getitem__(self, key):
            if key not in self.keys():
                raise RuntimeError(f"\"{key}\" is not found in the current step data. Followings are the available data: \n\t" + "\n\t".join(self.keys()))
            return super().__getitem__(key)

        def HasSolutionStepDataKey(self, key):
            return key in self.keys()

    def __init__(self, echo_level: int = 0):
        # objects registry
        self.__optimization_routines = {}

        # initializing the buffer variables
        self.__iteration_data = []
        self.__buffer_index = 0
        self.__echo_level = echo_level

    def AddOptimizationRoutine(self, routine: OptimizationRoutine):
        if not isinstance(routine, OptimizationRoutine):
            raise RuntimeError(f"Only objects of derrived types of OptimizationRoutine can be added. The object being added is: \n{routine}")

        class_hierrachy = getmro(routine.__class__)
        routine_index = class_hierrachy.index(OptimizationRoutine)

        base_class_name = class_hierrachy[routine_index - 1].__name__
        if not self.HasOptimizationRoutineType(base_class_name):
            self.__optimization_routines[base_class_name] = []

        if self.HasOptimizationRoutine(base_class_name, routine.GetName()):
            raise RuntimeError(f"Adding a routine with the name \"{routine.GetName()}\" while having another routine with the same name. OptimizationRoutine names for each type of object should be unique. [ Type of the routine: \"{base_class_name}\" ]")

        self.__optimization_routines[base_class_name].append(routine)
        self.__PrintInfo(1, f"Added {routine.GetName()} optimization routine.")

    def HasOptimizationRoutine(self, routine_class_type_name: str, routine_name: str) -> bool:
        if not self.HasOptimizationRoutineType(routine_class_type_name):
            raise RuntimeError(f"No objects of type \"{routine_class_type_name}\" [ requested routine_name = \"{routine_name}\" ]. Following are the available type options: \n\t" + "\n\t".join(self.__optimization_routines.keys()))

        return routine_name in [v.GetName() for v in self.__optimization_routines[routine_class_type_name]]

    def GetOptimizationRoutine(self, routine_class_type_name: str, routine_name: str) -> 'routine_class_type_name':
        if not self.HasOptimizationRoutineType(routine_class_type_name):
            raise RuntimeError(f"No objects of type \"{routine_class_type_name}\" [ requested routine_name = \"{routine_name}\" ]. Following are the available type options: \n\t" + "\n\t".join(self.__optimization_routines.keys()))

        if not self.HasOptimizationRoutine(routine_class_type_name, routine_name):
            raise RuntimeError(f"No routine with \"{routine_name}\" is available in the routines list for objects with type \"{routine_class_type_name}\". Followings are available options:\n\t" + "\n\t".join(v.GetName() for v in self.__optimization_routines[routine_class_type_name]))

        return self.__optimization_routines[routine_class_type_name][[v.GetName() for v in self.__optimization_routines[routine_class_type_name]].index(routine_name)]

    def HasOptimizationRoutineType(self, routine_class_type_name: str) -> bool:
        return routine_class_type_name in self.__optimization_routines.keys()

    def GetOptimizationRoutines(self, routine_class_type_name: str) -> 'list[OptimizationRoutine]':
        if not self.HasOptimizationRoutineType(routine_class_type_name):
            raise RuntimeError(f"No objects of type \"{routine_class_type_name}\". Following are the available type options: \n\t" + "\n\t".join(self.__optimization_routines.keys()))

        return self.__optimization_routines[routine_class_type_name]

    def SetBufferSize(self, buffer_size: int):
        if len(self.__iteration_data) != buffer_size:
            if len(self.__iteration_data) != 0:
                is_empty = True
                for v in self.__iteration_data:
                    if v != {}:
                        is_empty = False
                        break
                if not is_empty:
                    Kratos.Logger.PrintWarning(self.__class__.__name__, f"Changing buffer size with data will lose all the data in optimization info buffer. [ current_buffer_size = {len(self.__iteration_data)}, new_buffer_size = {buffer_size} ].")
            self.__iteration_data = [OptimizationInfo.__SolutionStepData() for _ in range(buffer_size)]
            self.__buffer_index = 0
            self.__PrintInfo(1, f"Optimization info buffer is set to {buffer_size}.")

    def AdvanceSolutionStep(self):
        self.__buffer_index = (self.__buffer_index + 1) % len(self.__iteration_data)

    def GetBufferSize(self) -> int:
        return len(self.__iteration_data)

    def HasSolutionStepDataKey(self, key):
        return self.GetSolutionStepData(0).HasSolutionStepDataKey(key)

    def GetSolutionStepData(self, solution_step_index: int) -> __SolutionStepData:
        if solution_step_index < 0:
            raise RuntimeError(f"solution_step_index should be positive. [ solution_Step_index = {solution_step_index} ].")
        current_step_index = (self.__buffer_index - solution_step_index) % len(self.__iteration_data)
        return self.__iteration_data[current_step_index]

    def __getitem__(self, key):
        return self.GetSolutionStepData(0)[key]

    def __setitem__(self, key, v):
        self.GetSolutionStepData(0)[key] = v

    def __PrintInfo(self, required_echo_level, msg, title = "OptimizationInfo"):
        if self.__echo_level >= required_echo_level:
            Kratos.Logger.PrintInfo(title, msg)

from inspect import getmro

import KratosMultiphysics as Kratos
from KratosMultiphysics.OptimizationApplication.routine import Routine

class OptimizationInfo:
    def __init__(self):
        # objects registry
        self.__objects = {}

        # initializing the buffer to one
        self.__iteration_data = [{}]
        self.__buffer_index = 0

    def AddRoutine(self, routine: Routine):
        if not isinstance(routine, Routine):
            raise RuntimeError(f"Only objects of derrived types of Routine can be added. The object being added is: \n{routine}")

        class_hierrachy = getmro(routine.__class__)
        routine_index = class_hierrachy.index(Routine)

        base_class = class_hierrachy[routine_index - 1]
        if base_class not in self.__objects.keys():
            self.__objects[base_class] = {}

        routines_dict: dict = self.__objects[base_class]
        if routine.GetName() in routines_dict.keys():
            raise RuntimeError(f"Adding a routine with the name \"{routine.GetName()}\" while having another routine with the same name. Routine names for each type of object should be unique. [ Type of the routine: \"{base_class.__name__}\" ]")

        routines_dict[routine.GetName()] = routine

    def GetRoutine(self, routine_class_type, routine_name: str) -> 'routine_class_type':
        if routine_class_type not in self.__objects.keys():
            raise RuntimeError(f"No objects of type \"{routine_class_type.__name__}\" [ requested routine_name = \"{routine_name}\" ].")

        if routine_name not in self.__objects[routine_class_type].keys():
            raise RuntimeError(f"No routine with \"{routine_name}\" is available in the routines list for objects with type \"{routine_class_type.__name__}\". Followings are available options:\n\t" + "\n\t".join(v.GetName() for v in self.__objects[routine_class_type].values()))

        return self.__objects[routine_class_type][routine_name]

    def GetRoutines(self, routine_class_type):
        if routine_class_type not in self.__objects.keys():
            raise RuntimeError(f"No objects of type \"{routine_class_type.__name__}\" [ requested routine_name = \"{routine_name}\" ].")

        return self.__objects[routine_class_type].values()

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
            self.__iteration_data = [{} for i in range(buffer_size)]
            self.__buffer_index = 0
            Kratos.Logger.PrintInfo(self.__class__.__name__, f"Optimization info buffer is set to {buffer_size}.")

    def Initialize(self):
        self["step"] = 1
        self.__ExecuteRoutineMethod("Initialize")

    def InitializeSolutionStep(self):
        self.__ExecuteRoutineMethod("InitializeSolutionStep")

    def FinalizeSolutionStep(self):
        self.__ExecuteRoutineMethod("FinalizeSolutionStep")
        self.__buffer_index = (self.__buffer_index + 1) % len(self.__iteration_data)
        self["step"] = self["step", 1] + 1

    def Finalize(self):
        self.__ExecuteRoutineMethod("Finalize")

    def __getitem__(self, key):
        if isinstance(key, tuple):
            if len(key) == 2:
                return self.GetData(key[1])[key[0]]
            else:
                raise RuntimeError(f"Unsupported key with length higher than 2 is provided for OptimizationInfo::__getitem__. [ key = {key} ].")
        else:
            return self.GetData(0)[key]

    def __setitem__(self, key, v):
        if isinstance(key, tuple):
            if len(key) == 2:
                self.GetData(key[1])[key[0]] = v
            else:
                raise RuntimeError(f"Unsupported key with length higher than 2 is provided for OptimizationInfo::__getitem__. [ key = {key} ].")
        else:
            self.GetData(0)[key] = v

    def Has(self, key):
        if isinstance(key, tuple):
            if len(key) == 2:
                return key[0] in self.GetData(key[1]).keys()
            else:
                raise RuntimeError(f"Unsupported key with length higher than 2 is provided for OptimizationInfo::__getitem__. [ key = {key} ].")
        else:
            return key in self.GetData(0).keys()

    def GetBufferSize(self) -> int:
        return len(self.__iteration_data)

    def GetData(self, solution_step_index: int) -> dict:
        if solution_step_index < 0:
            raise RuntimeError(f"solution_step_index should be positive. [ solution_Step_index = {solution_step_index} ].")
        current_step_index = (self.__buffer_index - solution_step_index) % len(self.__iteration_data)
        return self.__iteration_data[current_step_index]

    def __ExecuteRoutineMethod(self, method_name: str):
        for routines_dict in self.__objects.values():
            for routine in routines_dict.values():
                getattr(routine, method_name)()
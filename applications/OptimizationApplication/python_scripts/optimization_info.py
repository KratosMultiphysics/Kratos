from inspect import getmro

import KratosMultiphysics as Kratos
from KratosMultiphysics.OptimizationApplication.optimization_routine import OptimizationRoutine

class OptimizationInfo:
    def __init__(self):
        # objects registry
        self.__optimization_routines = {}

        # initializing the buffer to one
        self.__iteration_data = []
        self.__buffer_index = 0

    def AddRoutine(self, routine: OptimizationRoutine):
        if not isinstance(routine, OptimizationRoutine):
            raise RuntimeError(f"Only objects of derrived types of OptimizationRoutine can be added. The object being added is: \n{routine}")

        class_hierrachy = getmro(routine.__class__)
        routine_index = class_hierrachy.index(OptimizationRoutine)

        base_class_name = class_hierrachy[routine_index - 1].__name__
        if not self.HasRoutineType(base_class_name):
            self.__optimization_routines[base_class_name] = []

        if self.HasRoutine(base_class_name, routine.GetName()):
            raise RuntimeError(f"Adding a routine with the name \"{routine.GetName()}\" while having another routine with the same name. OptimizationRoutine names for each type of object should be unique. [ Type of the routine: \"{base_class_name}\" ]")

        self.__optimization_routines[base_class_name].append(routine)

    def HasRoutine(self, routine_class_type_name: str, routine_name: str) -> bool:
        if not self.HasRoutineType(routine_class_type_name):
            raise RuntimeError(f"No objects of type \"{routine_class_type_name}\" [ requested routine_name = \"{routine_name}\" ]. Following are the available type options: \n\t" + "\n\t".join(self.__optimization_routines.keys()))

        return routine_name in [v.GetName() for v in self.__optimization_routines[routine_class_type_name]]

    def GetRoutine(self, routine_class_type_name: str, routine_name: str) -> 'routine_class_type_name':
        if not self.HasRoutineType(routine_class_type_name):
            raise RuntimeError(f"No objects of type \"{routine_class_type_name}\" [ requested routine_name = \"{routine_name}\" ]. Following are the available type options: \n\t" + "\n\t".join(self.__optimization_routines.keys()))

        if not self.HasRoutine(routine_class_type_name, routine_name):
            raise RuntimeError(f"No routine with \"{routine_name}\" is available in the routines list for objects with type \"{routine_class_type_name}\". Followings are available options:\n\t" + "\n\t".join(v.GetName() for v in self.__optimization_routines[routine_class_type_name]))

        return self.__optimization_routines[routine_class_type_name][[v.GetName() for v in self.__optimization_routines[routine_class_type_name]].index(routine_name)]

    def HasRoutineType(self, routine_class_type_name: str) -> bool:
        return routine_class_type_name in self.__optimization_routines.keys()

    def GetRoutines(self, routine_class_type_name: str) -> list[OptimizationRoutine]:
        if not self.HasRoutineType(routine_class_type_name):
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
            self.__iteration_data = [{} for i in range(buffer_size)]
            self.__buffer_index = 0
            Kratos.Logger.PrintInfo(self.__class__.__name__, f"Optimization info buffer is set to {buffer_size}.")

    def AdvanceSolutionStep(self):
        self.__buffer_index = (self.__buffer_index + 1) % len(self.__iteration_data)

    def __getitem__(self, key):
        if isinstance(key, tuple):
            if len(key) == 2:
                data = self.GetData(key[1])
                if key[0] not in data.keys():
                    raise RuntimeError(f"\"{key[0]}\" is not found in the current step data. Followings are the available data: \n\t" + "\n\t".join(data.keys()))
                return data[key[0]]
            else:
                raise RuntimeError(f"Unsupported key with length higher than 2 is provided for OptimizationInfo::__getitem__. [ key = {key} ].")
        else:
            data = self.GetData(0)
            if key not in data.keys():
                raise RuntimeError(f"\"{key}\" is not found in the current step data. Followings are the available data: \n\t" + "\n\t".join(data.keys()))
            return data[key]

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

import KratosMultiphysics as Kratos

class OptimizationInfo:
    class __KeyErrorMessage(str):
        def __repr__(self):
            return str(self)

    def __init__(self, echo_level: int = 0):
        # objects registry
        self.__optimization_processes: 'dict[str, Kratos.Process]' = {}

        # initializing the buffer variables
        self.__iteration_data = []
        self.__buffer_index = 0
        self.__echo_level = echo_level

    def AddOptimizationProcess(self, category: any, routine_name: str, routine: Kratos.Process):
        if not isinstance(routine, Kratos.Process):
            raise RuntimeError(f"Only allowed to add objects derrived from Kratos.Process. [ Type of added routine: {routine.__class__.__name__}]")

        if not isinstance(routine, category):
            raise RuntimeError(f"The provided routine with name \"{routine_name}\" is not derrived from the provided category type \"{category.__name__}\".")

        if not self.HasOptimizationProcessType(category):
            self.__optimization_processes[category] = []

        if self.HasOptimizationProcess(category, routine_name):
            raise RuntimeError(f"Adding a routine with the name \"{routine_name}\" while having another routine with the same name. Kratos.Process names for each type of object should be unique. [ Type of the routine: \"{category.__name__}\" ]")

        self.__optimization_processes[category].append([routine_name, routine])
        self.__PrintInfo(1, f"Added {routine_name} optimization routine.")

    def HasOptimizationProcess(self, category: any, routine_name: str) -> bool:
        if not self.HasOptimizationProcessType(category):
            raise RuntimeError(f"No objects of type \"{category.__name__}\" [ requested routine_name = \"{routine_name}\" ]. Following are the available type options: \n\t" + "\n\t".join(self.__optimization_processes.keys()))

        return routine_name in [v[0] for v in self.__optimization_processes[category]]

    def GetOptimizationProcess(self, category: any, routine_name: str) -> Kratos.Process:
        if not self.HasOptimizationProcessType(category):
            raise RuntimeError(f"No objects of type \"{category.__name__}\" [ requested routine_name = \"{routine_name}\" ]. Following are the available type options: \n\t" + "\n\t".join(self.__optimization_processes.keys()))

        if not self.HasOptimizationProcess(category, routine_name):
            raise RuntimeError(f"No routine with \"{routine_name}\" is available in the routines list for objects with type \"{category.__name__}\". Followings are available options:\n\t" + "\n\t".join(v[0] for v in self.__optimization_processes[category]))

        return self.__optimization_processes[category][[v[0] for v in self.__optimization_processes[category]].index(routine_name)][1]

    def HasOptimizationProcessType(self, category: any) -> bool:
        return category in self.__optimization_processes.keys()

    def GetOptimizationProcesses(self, category: any) -> 'list[Kratos.Process]':
        if not self.HasOptimizationProcessType(category):
            raise RuntimeError(f"No objects of type \"{category.__name__}\". Following are the available type options: \n\t" + "\n\t".join(self.__optimization_processes.keys()))

        return [v[1] for v in self.__optimization_processes[category]]

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
            self.__iteration_data = [{} for _ in range(buffer_size)]
            self.__buffer_index = 0
            self.__PrintInfo(1, f"Optimization info buffer is set to {buffer_size}.")

    def AdvanceSolutionStep(self):
        self.__buffer_index = (self.__buffer_index + 1) % len(self.__iteration_data)

    def GetBufferSize(self) -> int:
        return len(self.__iteration_data)

    def GetSolutionStepData(self, solution_step_index: int) -> dict:
        if solution_step_index < 0:
            raise RuntimeError(f"solution_step_index should be positive. [ solution_Step_index = {solution_step_index} ].")
        current_step_index = (self.__buffer_index - solution_step_index) % len(self.__iteration_data)
        return self.__iteration_data[current_step_index]

    def HasValue(self, key: str, solution_step_index: int = 0) -> bool:
        if not isinstance(key, str):
            raise KeyError(f"Buffer initialization is only supported for string keys. Provided key = {key}.")

        current_opt_info = self.GetSolutionStepData(solution_step_index)

        sub_keys = key.split("/")
        for sub_key in sub_keys:
            if not sub_key in current_opt_info.keys():
                return False
            current_opt_info = current_opt_info[sub_key]

        return True

    def GetValue(self, key: str, solution_step_index: int = 0) -> any:
        if not isinstance(key, str):
            raise KeyError(f"Buffer initialization is only supported for string keys. Provided key = {key}.")

        current_opt_info = self.GetSolutionStepData(solution_step_index)

        sub_keys = key.split("/")
        for i, sub_key in enumerate(sub_keys):
            try:
                current_opt_info = current_opt_info[sub_key]
            except KeyError as k:
                cur_path = "\n\t" + "/".join(sub_keys[:i]) + "/"
                msg = f"Subkey {sub_key} not found for {key} for solution step index = {solution_step_index}. Followings are available keys:{cur_path}" + (cur_path).join(current_opt_info.keys())
                raise KeyError(OptimizationInfo.__KeyErrorMessage(msg)) from k

        return current_opt_info

    def SetValue(self, key: str, value: any, solution_step_index: int = 0, overwrite = False):
        if not isinstance(key, str):
            raise KeyError(f"Buffer initialization is only supported for string keys. Provided key = {key}.")

        current_opt_info = self.GetSolutionStepData(solution_step_index)

        sub_keys = key.split("/")
        for sub_key in sub_keys[:-1]:
            if sub_key not in current_opt_info.keys():
                current_opt_info[sub_key] = {}
            current_opt_info = current_opt_info[sub_key]

            if not isinstance(current_opt_info, dict):
                raise RuntimeError(f"Subkey \"{sub_key} of {key} is not a dictionary in solution step index = {solution_step_index}. Hence failing to add sub items.\"")

        if not overwrite and sub_keys[-1] in current_opt_info.keys():
            msg =  f"Value for {sub_keys[-1]} exists. [ Key = {key}, Current value = {current_opt_info[sub_keys[-1]]}, new value = {value}, current dict = {current_opt_info}]"
            raise RuntimeError(msg)

        current_opt_info[sub_keys[-1]] = value

    def __getitem__(self, key):
        return self.GetSolutionStepData(0)[key]

    def __setitem__(self, key, v):
        if not isinstance(key, str):
            raise KeyError(f"Buffer initialization is only supported for string keys. Provided key = {key}.")
        self.GetSolutionStepData(0)[key] = v

    def __PrintInfo(self, required_echo_level, msg, title = "OptimizationInfo"):
        if self.__echo_level >= required_echo_level:
            Kratos.Logger.PrintInfo(title, msg)

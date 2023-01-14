import KratosMultiphysics as Kratos

class OptimizationInfo:
    def __init__(self):
        self.__list_of_controllers = []
        self.__list_of_analysis = []
        self.__list_of_responses = []

        # initializing the buffer to one
        self.__iteration_data = [{}]
        self.__buffer_index = 0

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
        OptimizationInfo.__ExecuteMethod("Initialize", self.__list_of_analysis)
        OptimizationInfo.__ExecuteMethod("Initialize", self.__list_of_responses)
        OptimizationInfo.__ExecuteMethod("Initialize", self.__list_of_controllers)

    def InitializeSolutionStep(self):
        OptimizationInfo.__ExecuteMethod("InitializeSolutionStep", self.__list_of_analysis)
        OptimizationInfo.__ExecuteMethod("InitializeSolutionStep", self.__list_of_responses)
        OptimizationInfo.__ExecuteMethod("InitializeSolutionStep", self.__list_of_controllers)

    def FinalizeSolutionStep(self):
        OptimizationInfo.__ExecuteMethod("FinalizeSolutionStep", self.__list_of_analysis)
        OptimizationInfo.__ExecuteMethod("FinalizeSolutionStep", self.__list_of_responses)
        OptimizationInfo.__ExecuteMethod("FinalizeSolutionStep", self.__list_of_controllers)
        self.__buffer_index = (self.__buffer_index + 1) % len(self.__iteration_data)
        self["step"] = self["step", 1] + 1

    def Finalize(self):
        OptimizationInfo.__ExecuteMethod("Finalize", self.__list_of_analysis)
        OptimizationInfo.__ExecuteMethod("Finalize", self.__list_of_responses)
        OptimizationInfo.__ExecuteMethod("Finalize", self.__list_of_controllers)

    def AddExecutionPolicyWrapper(self, execution_policy_wrapper):
        OptimizationInfo.__AddObject(execution_policy_wrapper, self.__list_of_analysis)

    def AddResponseFunctionWrapper(self, response_function_wrapper):
        OptimizationInfo.__AddObject(response_function_wrapper, self.__list_of_responses)

    def AddControlsWrapper(self, controls_wrapper):
        OptimizationInfo.__AddObject(controls_wrapper, self.__list_of_controllers)

    def GetExecutionPolicyWrapper(self, name: str) -> "ExecutionPolicyWrapper":
        return OptimizationInfo.__GetObject(name, self.__list_of_analysis)

    def GetResponseFunctionWrapper(self, name: str) -> "ResponseFunctionBaseWrapper":
        return OptimizationInfo.__GetObject(name, self.__list_of_responses)

    def GetControlsWrapper(self, name: str) -> "ControlsWrapper":
        return OptimizationInfo.__GetObject(name, self.__list_of_controllers)

    def GetExecutionPolicyWrappers(self, name: str) -> "list[ExecutionPolicyWrappers]":
        return self.__list_of_analysis

    def GetResponseFunctionWrappers(self, name: str) -> "list[ResponseFunctionBaseWrappers]":
        return self.__list_of_responses

    def GetControlsWrappers(self, name: str) -> "list[ControlsWrapper]":
        return self.__list_of_controllers

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

    @staticmethod
    def __ExecuteMethod(method_name: str, objects_list: list):
        for v in objects_list:
            getattr(v, method_name)()

    @staticmethod
    def __AddObject(object, objects_list: list):
        for v in objects_list:
            if object.GetName() == v.GetName():
                raise RuntimeError(f"{object.GetName()} already exists. Please provide a unique name.")

        objects_list.append(object)

    @staticmethod
    def __GetObject(name: str, objects_list: list):
        for v in objects_list:
            if name == v.GetName():
                return v

        raise RuntimeError(f"{name} is not found. Followings are available options: \n\t" + "\n\t".join([item.GetName() for item in objects_list]))
import KratosMultiphysics as Kratos

class OptimizationInfo:
    def __init__(self):
        self.__list_of_controllers = []
        self.__list_of_analysis = []
        self.__list_of_responses = []

    def Initialize(self):
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
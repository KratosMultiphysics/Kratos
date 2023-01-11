import KratosMultiphysics as Kratos

class OptimizationInfo:
    def __init__(self):
        self.__dict_of_controllers = {}
        self.__dict_of_analysis = {}
        self.__dict_of_responses = {}

    def Initialize(self):
        self.__ExecuteMethod("Initialize", self.__dict_of_analysis)
        self.__ExecuteMethod("Initialize", self.__dict_of_responses)
        self.__ExecuteMethod("Initialize", self.__dict_of_controllers)

    def InitializeSolutionStep(self):
        self.__ExecuteMethod("InitializeSolutionStep", self.__dict_of_analysis)
        self.__ExecuteMethod("InitializeSolutionStep", self.__dict_of_responses)
        self.__ExecuteMethod("InitializeSolutionStep", self.__dict_of_controllers)

    def FinalizeSolutionStep(self):
        self.__ExecuteMethod("FinalizeSolutionStep", self.__dict_of_analysis)
        self.__ExecuteMethod("FinalizeSolutionStep", self.__dict_of_responses)
        self.__ExecuteMethod("FinalizeSolutionStep", self.__dict_of_controllers)

    def Finalize(self):
        self.__ExecuteMethod("Finalize", self.__dict_of_analysis)
        self.__ExecuteMethod("Finalize", self.__dict_of_responses)
        self.__ExecuteMethod("Finalize", self.__dict_of_controllers)

    def AddExecutionPolicyWrapper(self, execution_policy_wrapper):
        self.__AddObject(execution_policy_wrapper, self.__dict_of_analysis)

    def AddResponseFunctionWrapper(self, response_function_wrapper):
        self.__AddObject(response_function_wrapper, self.__dict_of_responses)

    def AddControlsWrapper(self, controls_wrapper):
        self.__AddObject(controls_wrapper, self.__dict_of_controllers)

    def GetExecutionPolicyWrapper(self, name: str) -> "ExecutionPolicyWrapper":
        return self.__GetObject(name, self.__dict_of_analysis)

    def GetResponseFunctionWrapper(self, name: str) -> "ResponseFunctionBaseWrapper":
        return self.__GetObject(name, self.__dict_of_responses)

    def GetControlsWrapper(self, name: str) -> "ControlsWrapper":
        return self.__GetObject(name, self.__dict_of_controllers)

    def __ExecuteMethod(self, method_name: str, objects_dict: dict):
        for v in objects_dict.values():
            getattr(v, method_name)()

    def __AddObject(self, object, objects_dict: dict):
        if object.GetName() in objects_dict.keys():
            raise RuntimeError(f"{object.GetName()} already exists. Please provide a unique names.")

        objects_dict[object.GetName()] = object

    def __GetObject(self, name: str, objects_dict: dict):
        if not name in objects_dict.keys():
            raise RuntimeError(f"{name} is not found. Followings are available options: \n\t" + "\n\t".join(objects_dict.keys()))

        return objects_dict[name]
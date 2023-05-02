from typing import Any
from KratosMultiphysics.OptimizationApplication.utilities.buffered_dict import BufferedDict
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem

class ComponentData:
    def __init__(self, component: Any, optimization_problem: OptimizationProblem):
        self.__problem_data = optimization_problem.GetProblemDataContainer()

        self.__component_type = optimization_problem.GetComponentType(component)
        self.__component_name = optimization_problem.GetComponentName(component)

        self.__data_name = f"{self.__component_type.__name__}/{self.__component_name}"

        # if the component data container not found, then create it.
        if not self.__problem_data.HasValue(self.__data_name):
            self.ResetData()

        self.__component = component

    def ResetData(self):
        self.__problem_data[self.__data_name] = {}

        self.__component_data: BufferedDict = self.__problem_data[self.__data_name]
        self.__component_buffered_data: BufferedDict = None

        # create the unbuffered data container
        self.__component_unbuffered_data = BufferedDict(1)
        self.__component_data["unbuffered"] = self.__component_unbuffered_data

    def SetDataBuffer(self, buffer_size: int):
        if not self.__component_data.HasValue("buffered"):
            self.__component_buffered_data = BufferedDict(buffer_size)
            self.__component_data["buffered"] = self.__component_buffered_data
        else:
            self.__component_buffered_data: BufferedDict = self.__component_data["buffered"]
            if self.__component_buffered_data.GetBufferSize() < buffer_size:
                raise RuntimeError(f"The required buffer size is not satisfied with the existing problem data container. [ response data container buffer size = {self.__response_buffered_data.GetBufferSize()}, required buffer size = {required_buffer_size}")

    def GetComponent(self) -> Any:
        return self.__component

    def GetBufferedData(self) -> BufferedDict:
        if self.__component_buffered_data is None:
            raise RuntimeError(f"Buffered data is not set by calling ComponentData::SetBuffer for component of type \"{self.__component_type}\" with component name \"{self.__component_name}\".")

        return self.__component_buffered_data

    def GetUnBufferedData(self) -> BufferedDict:
        return self.__component_unbuffered_data
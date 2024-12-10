from typing import Any
from KratosMultiphysics.OptimizationApplication.utilities.buffered_dict import BufferedDict
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem

class ComponentDataView:
    """A viewer to view component's data in optimization problem.

    Instances of this class create or view provided component's data in
    the OptimizationProblem container. It will view both buffered and
    unbuffered data.

    The component must be first added to the OptimizationProblem
    to use instances of this class to modify component data.

    """
    def __init__(self, component: Any, optimization_problem: OptimizationProblem):
        self.__problem_data = optimization_problem.GetProblemDataContainer()

        if isinstance(component, str):
            self.__component_type = object
            self.__component_name = component
        else:
            self.__component_type = optimization_problem.GetComponentType(component)
            self.__component_name = optimization_problem.GetComponentName(component)

        self.__data_name = f"{self.__component_type.__name__}/{self.__component_name}"
        self.__buffered_data_name = f"{self.__data_name}/buffered"
        self.__unbuffered_data_name = f"{self.__data_name}/unbuffered"
        self.__component = component

        # if the component data container not found, then create it.
        if not self.__problem_data.HasValue(self.__data_name):
            self.ResetData()

    def ResetData(self):
        if self.__problem_data.HasValue(self.__data_name):
            # if exists, delete the current problem data
            del self.__problem_data[self.__data_name]

        self.__problem_data[self.__data_name] = BufferedDict(1)
        # create the unbuffered data container
        self.__problem_data[self.__unbuffered_data_name] = BufferedDict(1, False)

    def SetDataBuffer(self, buffer_size: int):
        if not self.__problem_data.HasValue(self.__buffered_data_name):
            self.__problem_data[self.__buffered_data_name] = BufferedDict(buffer_size)
        else:
            buffered_data: BufferedDict = self.__problem_data[self.__buffered_data_name]
            if buffered_data.GetBufferSize() < buffer_size:
                raise RuntimeError(f"The required buffer size is not satisfied with the existing problem data container. [ component data container buffer size = {buffered_data.GetBufferSize()}, required buffer size = {buffer_size}, component = {self.__component_name}, component type = {self.__component_type}")

    def HasDataBuffer(self) -> bool:
        return self.__problem_data.HasValue(self.__buffered_data_name)

    def GetComponent(self) -> Any:
        return self.__component

    def GetBufferedData(self) -> BufferedDict:
        if not self.HasDataBuffer():
            raise RuntimeError(f"Buffered data is not set by calling ComponentData::SetDataBuffer for component of type \"{self.__component_type}\" with component name \"{self.__component_name}\".")

        return self.__problem_data[self.__buffered_data_name]

    def GetUnBufferedData(self) -> BufferedDict:
        return self.__problem_data[self.__unbuffered_data_name]

    def GetDataPath(self) -> str:
        return self.__data_name

    def GetComponentName(self):
        if isinstance(self.__component, str):
            return self.__component
        else:
            return self.__component.GetName()
        
    def RemoveComponentData(self, data_name: str):
        if self.__problem_data.HasValue( data_name):
            del self.__problem_data[ data_name]
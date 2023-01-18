from enum import Enum
import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA

class ContainerData:
    class ContainerEnum(Enum):
        NODES = 1
        ELEMENTS = 2
        CONDITIONS = 3
        ELEMENT_PROPERTIES = 4
        CONDITION_PROPERTIES = 5

    def __init__(self, model_part: Kratos.ModelPart, container_type: ContainerEnum):
        self.__model_part = model_part
        self.__container_type = container_type
        self.__data = Kratos.Vector()

    def GetModelPart(self) -> Kratos.ModelPart:
        return self.__model_part

    def GetContainerTpe(self):
        return self.__container_type

    def GetData(self) -> Kratos.Vector:
        return self.__data

    def SetData(self, data: Kratos.Vector):
        self.__data = data

    def Clone(self):
        other = ContainerData(self.__model_part, self.__container_type)
        other.__data = Kratos.Vector(self.__data)
        return other

    def IsSameContainer(self, other):
        return isinstance(other, ContainerData) and self.__model_part == other.__model_part and self.__container_type == other.__container_type

    def ReadDataFromContainer(self, variable: any):
        domain_size = self.__model_part.ProcessInfo[Kratos.DOMAIN_SIZE]
        if self.__container_type == ContainerData.ContainerEnum.NODES:
            KratosOA.OptimizationUtils.GetContainerVariableToVector(self.__model_part.Nodes, variable, domain_size, self.__data)
        elif self.__container_type == ContainerData.ContainerEnum.CONDITIONS:
            KratosOA.OptimizationUtils.GetContainerVariableToVector(self.__model_part.Conditions, variable, domain_size, self.__data)
        elif self.__container_type == ContainerData.ContainerEnum.ELEMENTS:
            KratosOA.OptimizationUtils.GetContainerVariableToVector(self.__model_part.Elements, variable, domain_size, self.__data)
        elif self.__container_type == ContainerData.ContainerEnum.CONDITION_PROPERTIES:
            KratosOA.OptimizationUtils.GetContainerPropertiesVariableToVector(self.__model_part.Conditions, variable, self.__data)
        elif self.__container_type == ContainerData.ContainerEnum.ELEMENT_PROPERTIES:
            KratosOA.OptimizationUtils.GetContainerPropertiesVariableToVector(self.__model_part.Elements, variable, self.__data)

    def AssignDataToContainer(self, variable: any):
        domain_size = self.__model_part.ProcessInfo[Kratos.DOMAIN_SIZE]
        if self.__container_type == ContainerData.ContainerEnum.NODES:
            KratosOA.OptimizationUtils.AssignVectorToContainer(self.__model_part.Nodes, variable, domain_size, self.__data)
        elif self.__container_type == ContainerData.ContainerEnum.CONDITIONS:
            KratosOA.OptimizationUtils.AssignVectorToContainer(self.__model_part.Conditions, variable, domain_size, self.__data)
        elif self.__container_type == ContainerData.ContainerEnum.ELEMENTS:
            KratosOA.OptimizationUtils.AssignVectorToContainer(self.__model_part.Elements, variable, domain_size, self.__data)
        elif self.__container_type == ContainerData.ContainerEnum.CONDITION_PROPERTIES:
            KratosOA.OptimizationUtils.AssignVectorToContainerProperties(self.__model_part.Conditions, variable, self.__data)
        elif self.__container_type == ContainerData.ContainerEnum.ELEMENT_PROPERTIES:
            KratosOA.OptimizationUtils.AssignVectorToContainerProperties(self.__model_part.Elements, variable, self.__data)

    def NormInf(self) -> float:
        return KratosOA.OptimizationUtils.NormInf(self.__data)

    def __add__(self, other):
        if not isinstance(other, ContainerData):
            raise NotImplementedError(f"Adding an object of ContainerData to non object of ContainerData is not supprted. [ Type of right hand value = {type(other)} ].")

        if not self.IsSameContainer(other):
            raise RuntimeError(f"Trying to add {other} to {self} which are not of the same container.")

        result = ContainerData(self.__model_part, self.__container_type)
        KratosOA.OptimizationUtils.AddVectors(result.__data, self.__data, other.__data)
        return result

    def __sub__(self, other):
        if not isinstance(other, ContainerData):
            raise NotImplementedError(f"Substracting an object of ContainerData to non object of ContainerData is not supprted. [ Type of right hand value = {type(other)} ].")

        if not self.IsSameContainer(other):
            raise RuntimeError(f"Trying to substract {other} from {self} which are not of the same container.")

        result = ContainerData(self.__model_part, self.__container_type)
        KratosOA.OptimizationUtils.SubstractVectors(result.__data, self.__data, other.__data)
        return result

    def __mul__(self, other):
        if not isinstance(other, (float, int)):
            raise NotImplementedError(f"Multiplying an object of ContainerData with non float or int is not supprted. [ Type of right hand value = {type(other)} ].")

        result = ContainerData(self.__model_part, self.__container_type)
        KratosOA.OptimizationUtils.MultiplyVector(result.__data, self.__data, float(other))
        return result

    def __truediv__(self, other):
        if not isinstance(other, (float, int)):
            raise NotImplementedError(f"Dividing an object of ContainerData with non float or int is not supprted. [ Type of right hand value = {type(other)} ].")

        result = ContainerData(self.__model_part, self.__container_type)
        KratosOA.OptimizationUtils.DivideVector(result.__data, self.__data, float(other))
        return result

    def __str__(self) -> str:
        return f"ContainerData[ Model part name = {self.__model_part.FullName()}, container type: {self.__container_type.name}, size of values = {self.__data.Size()} ]"
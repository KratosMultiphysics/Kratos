from abc import ABC
from abc import abstractmethod

import KratosMultiphysics as Kratos
from KratosMultiphysics.OptimizationApplication.optimization_info import OptimizationInfo
from KratosMultiphysics.OptimizationApplication.optimization_routine import OptimizationRoutine
from KratosMultiphysics.OptimizationApplication.utilities.control_transformation_technique import ControlTransformationTechnique
from KratosMultiphysics.OptimizationApplication.utilities.response_function_implementor import ObjectiveResponseFunctionImplementor
from KratosMultiphysics.OptimizationApplication.utilities.response_function_implementor import ConstraintResponseFunctionImplementor
from KratosMultiphysics.OptimizationApplication.utilities.helper_utils import CallOnAll

class Algorithm(OptimizationRoutine, ABC):
    def __init__(self, model: Kratos.Model, parameters: Kratos.Parameters, optimization_info: OptimizationInfo):
        super().__init__(model, parameters, optimization_info)
        self.__list_of_objectives: 'list[ObjectiveResponseFunctionImplementor]' = []
        self.__list_of_constraints: 'list[ConstraintResponseFunctionImplementor]' = []
        self.__list_of_controllers: 'list[ControlTransformationTechnique]' = []
        self.__name = None

    def SetName(self, name: str):
        self.__name = name

    def GetName(self) -> str:
        if self.__name is None:
            raise RuntimeError("Algorithm name is not set.")

        return self.__name

    def AddVariables(self):
        pass

    def AddDofs(self):
        pass

    def InitializeSolutionStep(self):
        CallOnAll(self.GetObjectives(), ObjectiveResponseFunctionImplementor.ResetResponseData)
        CallOnAll(self.GetConstraints(), ConstraintResponseFunctionImplementor.ResetResponseData)

    def SetObjectives(self, list_of_objectives: 'list[ObjectiveResponseFunctionImplementor]'):
        self.__list_of_objectives = list_of_objectives

    def SetConstraints(self, list_of_constraints: 'list[ConstraintResponseFunctionImplementor]'):
        self.__list_of_constraints = list_of_constraints

    def SetControllers(self, list_of_controllers: 'list[ControlTransformationTechnique]'):
        self.__list_of_controllers = list_of_controllers

    def GetObjectives(self) -> 'list[ObjectiveResponseFunctionImplementor]':
        return self.__list_of_objectives

    def GetConstraints(self) -> 'list[ConstraintResponseFunctionImplementor]':
        return self.__list_of_constraints

    def GetControllers(self) -> 'list[ControlTransformationTechnique]':
        return self.__list_of_controllers

    @abstractmethod
    def GetMinimumBufferSize(self) -> int:
        pass

    @abstractmethod
    def Check(self):
        pass

    @abstractmethod
    def SolveSolutionStep(self):
        pass

    @abstractmethod
    def IsConverged(self) -> bool:
        pass
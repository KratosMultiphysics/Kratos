from abc import ABC
from abc import abstractmethod

import KratosMultiphysics as Kratos
from KratosMultiphysics.python_solver import PythonSolver
from KratosMultiphysics.OptimizationApplication.optimization_info import OptimizationInfo
from KratosMultiphysics.OptimizationApplication.optimization_routine import OptimizationRoutine
from KratosMultiphysics.OptimizationApplication.model_part_controllers.model_part_controller import ModelPartController
from KratosMultiphysics.OptimizationApplication.utilities.control_transformation_technique import ControlTransformationTechnique
from KratosMultiphysics.OptimizationApplication.utilities.response_function_implementor import ObjectiveResponseFunctionImplementor
from KratosMultiphysics.OptimizationApplication.utilities.response_function_implementor import ConstraintResponseFunctionImplementor
from KratosMultiphysics.OptimizationApplication.utilities.helper_utils import CallOnAll

class Algorithm(PythonSolver, OptimizationRoutine, ABC):
    def __init__(self, model: Kratos.Model, parameters: Kratos.Parameters, optimization_info: OptimizationInfo):
        PythonSolver.__init__(self, model, parameters)
        self.optimization_info = optimization_info
        self.optimization_info.SetBufferSize(self.GetMinimumBufferSize())
        self.step = 0
        self.optimization_info["step"] = self.step

        self.__list_of_controllers = [self.optimization_info.GetOptimizationRoutine(ControlTransformationTechnique, control_name) for control_name in parameters["control_names"].GetStringArray()]
        self.__list_of_objectives = [ObjectiveResponseFunctionImplementor(objective_settings, self.optimization_info) for objective_settings in parameters["objectives"]]
        self.__list_of_constraints = [ConstraintResponseFunctionImplementor(constraint_settings, self.optimization_info) for constraint_settings in parameters["constraints"]]

    def ImportModelPart(self):
        CallOnAll(self.optimization_info.GetOptimizationRoutines(ModelPartController), ModelPartController.ImportModelPart)

    def AdvanceInTime(self, _):
        self.optimization_info.AdvanceSolutionStep()
        self.step += 1
        self.optimization_info["step"] = self.step

    def InitializeSolutionStep(self):
        CallOnAll(self.GetObjectives(), ObjectiveResponseFunctionImplementor.ResetResponseData)
        CallOnAll(self.GetConstraints(), ConstraintResponseFunctionImplementor.ResetResponseData)

    def GetObjectives(self) -> 'list[ObjectiveResponseFunctionImplementor]':
        return self.__list_of_objectives

    def GetConstraints(self) -> 'list[ConstraintResponseFunctionImplementor]':
        return self.__list_of_constraints

    def GetControllers(self) -> 'list[ControlTransformationTechnique]':
        return self.__list_of_controllers

    def GetComputingModelPart(self):
        model_part_name = self.settings["model_part_name"].GetString()
        if not self.model.HasModelPart(model_part_name):
            self.model.CreateModelPart(model_part_name)

        return self.model[model_part_name]

    @abstractmethod
    def GetMinimumBufferSize(self) -> int:
        pass

    @abstractmethod
    def Check(self):
        pass

    @abstractmethod
    def SolveSolutionStep(self) -> bool:
        pass

    @abstractmethod
    def IsConverged(self) -> bool:
        pass
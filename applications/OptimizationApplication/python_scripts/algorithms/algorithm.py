from abc import ABC
from abc import abstractmethod

import KratosMultiphysics as Kratos
from KratosMultiphysics.python_solver import PythonSolver
from KratosMultiphysics.OptimizationApplication.optimization_info import OptimizationInfo
from KratosMultiphysics.OptimizationApplication.model_part_controllers.model_part_controller import ModelPartController
from KratosMultiphysics.OptimizationApplication.utilities.control_transformation_technique import ControlTransformationTechnique
from KratosMultiphysics.OptimizationApplication.utilities.response_function_implementor import ObjectiveResponseFunctionImplementor
from KratosMultiphysics.OptimizationApplication.utilities.response_function_implementor import ConstraintResponseFunctionImplementor
from KratosMultiphysics.OptimizationApplication.utilities.helper_utilities import CallOnAll

class Algorithm(PythonSolver, ABC):
    def __init__(self, model: Kratos.Model, parameters: Kratos.Parameters, optimization_info: OptimizationInfo):
        PythonSolver.__init__(self, model, parameters)

        self.optimization_info = optimization_info
        self.optimization_info.SetBufferSize(self.GetMinimumBufferSize())
        self.step = 0
        self.optimization_info["step"] = self.step

        self.__list_of_control_transformation_techniques = [self.optimization_info.GetOptimizationProcess(ControlTransformationTechnique, control_name) for control_name in parameters["control_names"].GetStringArray()]

    def ImportModelPart(self):
        CallOnAll(self.optimization_info.GetOptimizationProcesses(ModelPartController), ModelPartController.ImportModelPart)

    def AdvanceInTime(self, _):
        self.step += 1

        # set the model part for the next iteration by applying update from the previous iteration
        if self.step > 1:
            CallOnAll(self.GetControlTransformationTechniques(), ControlTransformationTechnique.ApplyControlUpdate)

        # now we march into next iteration
        self.optimization_info.AdvanceSolutionStep()
        self.optimization_info["step"] = self.step

    def Initialize(self):
        self.__list_of_objectives = [ObjectiveResponseFunctionImplementor(objective_settings, self.optimization_info) for objective_settings in self.settings["objectives"]]
        self.__list_of_constraints = [ConstraintResponseFunctionImplementor(constraint_settings, self.optimization_info) for constraint_settings in self.settings["constraints"]]

    def GetObjectives(self) -> 'list[ObjectiveResponseFunctionImplementor]':
        return self.__list_of_objectives

    def GetConstraints(self) -> 'list[ConstraintResponseFunctionImplementor]':
        return self.__list_of_constraints

    def GetControlTransformationTechniques(self) -> 'list[ControlTransformationTechnique]':
        return self.__list_of_control_transformation_techniques

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
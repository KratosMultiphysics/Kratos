# --- Core Imports ---
import KratosMultiphysics

# --- Structural Imports ---
from KratosMultiphysics.StructuralMechanicsApplication import InsertDirichletPreTensionOperation, InsertNeumannPreTensionOperation

# --- STD Imports ---
from typing import Union
import abc


class __PreTensionProcessMeta(type(KratosMultiphysics.Process), type(abc.ABC)):
     pass


class PreTensionProcessBase(KratosMultiphysics.Process, abc.ABC, metaclass = __PreTensionProcessMeta):
    """ @brief Common base class of @ref NeumannPreTensionProcess and @ref DirichletPreTensionProcess.
        @see @ref InsertNeumannPreTensionOperation
        @see @ref InsertDirichletPreTensionProcess
    """

    def __init__(
        self,
        model: KratosMultiphysics.Model,
        parameters: KratosMultiphysics.Parameters) -> None:
            super().__init__()
            parameters.AddMissingParameters(self.GetDefaultParameters())

            # Strip the "magnitude" setting because it has to be handled separately,
            # since it can either be a numeric value or a string.
            magnitude_setting: KratosMultiphysics.Parameters = parameters["magnitude"].Clone()
            parameters.RemoveValue("magnitude")

            self.__magnitude: Union[float,KratosMultiphysics.GenericFunctionUtility]
            if magnitude_setting.IsDouble():
                self.__magnitude = magnitude_setting.GetDouble()
            elif magnitude_setting.IsString():
                self.__magnitude = KratosMultiphysics.GenericFunctionUtility(magnitude_setting.GetString())
                if self.__magnitude.DependsOnSpace():
                    raise ValueError(f"\"magnitude\" cannot depend on spatial variables, only (pseudo-) time")
            else:
                raise ValueError(f"expecting either a numeric value or a string defining a function for \"magnitude\", but got {magnitude_setting}")

            self.__model_part: KratosMultiphysics.ModelPart = model.GetModelPart(parameters["model_part_name"].GetString())
            self.__operation: Union[InsertDirichletPreTensionOperation,InsertNeumannPreTensionOperation] = self._MakeOperation(model, parameters)


    def ExecuteBeforeSolutionLoop(self) -> None:
        self.__operation.Execute()


    def ExecuteInitializeSolutionStep(self) -> None:
        self.__operation.Apply(self.__GetMagnitude())


    def GetDefaultParameters(self) -> KratosMultiphysics.Parameters:
        return KratosMultiphysics.Parameters(R"""{
            "model_part_name" : "",
            "magnitude" : 0.0,
            "verbosity" : 1
        }""")


    @classmethod
    @abc.abstractmethod
    def _MakeOperation(
        cls,
        model: KratosMultiphysics.Model,
        parameters: KratosMultiphysics.Parameters) -> Union[InsertDirichletPreTensionOperation,InsertNeumannPreTensionOperation]:
            pass


    def __GetMagnitude(self) -> float:
        if isinstance(self.__magnitude, float):
            return self.__magnitude
        elif isinstance(self.__magnitude, KratosMultiphysics.GenericFunctionUtility):
            return self.__magnitude.CallFunction(
                0.0,
                0.0,
                0.0,
                self.__model_part.ProcessInfo[KratosMultiphysics.TIME],
                0.0,
                0.0,
                0.0)
        else:
            raise RuntimeError(f"unhandled magnitude type {type(self.__magnitude)}")

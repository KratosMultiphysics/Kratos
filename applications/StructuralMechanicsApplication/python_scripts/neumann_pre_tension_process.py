# --- Core Imports ---
import KratosMultiphysics

# --- Structural Imports ---
from KratosMultiphysics import StructuralMechanicsApplication as StructuralMechanics


class NeumannPretensionProcess(KratosMultiphysics.Process):

    def __init__(self,
                 model: KratosMultiphysics.Model,
                 parameters: KratosMultiphysics.Parameters) -> None:
        super().__init__()
        parameters.ValidateAndAssignDefaults(self.GetDefaultParameters())
        self.__insert_pre_tension_operation = StructuralMechanics.InsertNeumannPreTensionOperation(
            model,
            parameters)

    def ExecuteBeforeSolutionLoop(self) -> None:
        self.__insert_pre_tension_operation.Execute()

    def GetDefaultParameters(self) -> KratosMultiphysics.Parameters:
        return KratosMultiphysics.Parameters(R"""{
            "model_part_name" : "",
            "magnitude" : 0.0,
            "verbosity" : 1
        }""")


def Factory(parameters: KratosMultiphysics.Parameters,
            model: KratosMultiphysics.Model) -> KratosMultiphysics.Process:
    return NeumannPretensionProcess(
        model,
        parameters["Parameters"])

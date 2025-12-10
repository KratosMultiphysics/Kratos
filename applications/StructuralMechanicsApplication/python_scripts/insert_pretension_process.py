# --- Core Imports ---
import KratosMultiphysics

# --- Structural Imports ---
from KratosMultiphysics import StructuralMechanicsApplication as StructuralMechanics


class PretensionProcess(KratosMultiphysics.Process):

    def __init__(self,
                 model: KratosMultiphysics.Model,
                 parameters: KratosMultiphysics.Parameters) -> None:
        parameters.ValidateAndAssignDefaults(self.GetDefaultParameters())
        #self.__model_part = model.GetModelPart(parameters["model_part_name"].GetString())
        self.__insert_pretension_operation = StructuralMechanics.InsertPretensionOperation(
            model,
            parameters)

    def ExecuteBeforeSolutionLoop(self) -> None:
        self.__insert_pretension_operation.Execute()

    def GetDefaultParameters(self) -> KratosMultiphysics.Parameters:
        return KratosMultiphysics.Parameters(R"""{
            "model_part_name" : "",
            "pretension_value" : 0.0
        }""")


def Factory(parameters: KratosMultiphysics.Parameters,
            model: KratosMultiphysics.Model) -> KratosMultiphysics.Process:
    return PretensionProcess(
        model,
        parameters["Parameters"])

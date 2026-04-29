# --- Core Imports ---
import KratosMultiphysics

# --- Structural Imports ---
from KratosMultiphysics import StructuralMechanicsApplication as StructuralMechanics
from KratosMultiphysics.StructuralMechanicsApplication.pre_tension_process_base import PreTensionProcessBase, InsertNeumannPreTensionOperation


class NeumannPreTensionProcess(PreTensionProcessBase):

    def __init__(self,
                 model: KratosMultiphysics.Model,
                 parameters: KratosMultiphysics.Parameters) -> None:
        super().__init__(model, parameters)

    @classmethod
    def _MakeOperation(
        cls,
        model: KratosMultiphysics.Model,
        parameters: KratosMultiphysics.Parameters) -> InsertNeumannPreTensionOperation:
            return InsertNeumannPreTensionOperation(model, parameters)


def Factory(parameters: KratosMultiphysics.Parameters,
            model: KratosMultiphysics.Model) -> KratosMultiphysics.Process:
    return NeumannPreTensionProcess(
        model,
        parameters["Parameters"])

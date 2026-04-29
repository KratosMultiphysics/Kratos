# --- Core Imports ---
import KratosMultiphysics

# --- Structural Imports ---
from KratosMultiphysics import StructuralMechanicsApplication as StructuralMechanics
from KratosMultiphysics.StructuralMechanicsApplication.pre_tension_process_base import PreTensionProcessBase, InsertDirichletPreTensionOperation


class DirichletPretensionProcess(PreTensionProcessBase):

    def __init__(self,
                 model: KratosMultiphysics.Model,
                 parameters: KratosMultiphysics.Parameters) -> None:
        super().__init__(model, parameters)

    @classmethod
    def _MakeOperation(
        cls,
        model: KratosMultiphysics.Model,
        parameters: KratosMultiphysics.Parameters) -> InsertDirichletPreTensionOperation:
            return InsertDirichletPreTensionOperation(model, parameters)


def Factory(parameters: KratosMultiphysics.Parameters,
            model: KratosMultiphysics.Model) -> KratosMultiphysics.Process:
    return DirichletPretensionProcess(
        model,
        parameters["Parameters"])

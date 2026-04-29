# --- Core Imports ---
import KratosMultiphysics

# --- Structural Imports ---
from KratosMultiphysics import StructuralMechanicsApplication as StructuralMechanics
from KratosMultiphysics.StructuralMechanicsApplication.pre_tension_process_base import PreTensionProcessBase, InsertDirichletPreTensionOperation


class DirichletPretensionProcess(PreTensionProcessBase):
    """
        @brief Pre-tensioning defined by a surface and a prescribed displacement.
        @details Cut the mesh along the provided surface and fix the average
                 out-of-plane displacement to the provided value while forbidding relative
                 in-plane displacements.
        @see @ref pre_tensioning "Pre-Tensioning"
        @classname DirichletPreTensionProcess
        @ingroup pre_tensioning
    """

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

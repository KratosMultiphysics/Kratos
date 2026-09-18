from abc import ABC, abstractmethod
import KratosMultiphysics as Kratos

class ModelPartController(ABC):
    def Initialize(self) -> None:
        pass

    def Finalize(self) -> None:
        pass

    def SupportsRestart(self) -> bool:
        """Whether this controller's model part can be saved/loaded via restart_settings'
        model-part-level (Kratos.Serializer) checkpointing.

        Only controllers that own an independent, non-derived import of their model part (e.g.
        MdpaModelPartController) should return True. A derived/shared-node model part (e.g.
        ConnectivityPreservingModelPartController, which has no single "its model part" to begin
        with -- see its GetModelPart()) cannot safely be restart-loaded from its own serialized
        data without breaking the node-sharing invariant with the model part it derives from, or
        colliding with entities it unconditionally (re-)creates on every ImportModelPart() call --
        it should keep re-deriving itself every run instead, restart or not.

        OptimizationAnalysis._CreateRestart() loads restart-capable model parts directly and
        synchronously in __init__, before Initialize() ever calls ImportModelPart() -- so a
        restart-capable controller's ImportModelPart() only needs to check whether its model part
        is already populated (e.g. NumberOfNodes() > 0) and skip its normal import if so; no
        separate hook/flag is needed here.
        """
        return False

    @abstractmethod
    def ImportModelPart(self) -> None:
        pass

    @abstractmethod
    def GetModelPart(self) -> Kratos.ModelPart:
        pass



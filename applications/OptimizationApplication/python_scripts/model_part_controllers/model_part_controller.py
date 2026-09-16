from abc import ABC, abstractmethod
from pathlib import Path
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
        ConnectivityPreservingModelPartController) cannot safely be restart-loaded from its own
        separate file without breaking the node-sharing invariant with the model part it derives
        from -- it should keep re-deriving itself on every ImportModelPart() call instead.
        """
        return False

    def SetRestartLoadFile(self, file_path: Path, serializer_trace: Kratos.SerializerTraceType) -> None:
        """Called by OptimizationAnalysis._CreateRestart() before Initialize() invokes
        ImportModelPart(), when a restart checkpoint file was found for this controller's model
        part. Implementations must make the next ImportModelPart() call load from file_path (a
        Kratos.FileSerializer base path, i.e. without the ".rest" suffix) instead of their normal
        import. Only meaningful if SupportsRestart() is True; no-op otherwise.
        """
        pass

    @abstractmethod
    def ImportModelPart(self) -> None:
        pass

    @abstractmethod
    def GetModelPart(self) -> Kratos.ModelPart:
        pass



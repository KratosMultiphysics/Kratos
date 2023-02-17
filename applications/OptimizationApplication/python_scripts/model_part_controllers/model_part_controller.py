import KratosMultiphysics as Kratos

class ModelPartController(Kratos.Process):
    def __init__(self):
        super().__init__()

    def ImportModelPart(self):
        raise NotImplementedError("Calling base class ModelPartController::ImportModelPart. Please implement it in the derrived class.")

    def GetModelPart(self) -> Kratos.ModelPart:
        raise NotImplementedError("Calling base class ModelPartController::GetModelPart. Please implement it in the derrived class.")



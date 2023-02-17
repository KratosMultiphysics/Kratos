import KratosMultiphysics as Kratos

class ResponseFunction(Kratos.Process):
    def __init__(self):
        super().__init__()

    def Check(self):
        raise NotImplementedError("Calling base class ResponseFunction::Check method. Please implement in the derrived class.")

    def CalculateValue(self) -> float:
        raise NotImplementedError("Calling base class ResponseFunction::CalculateValue method. Please implement in the derrived class.")

    def CalculateSensitivity(self, sensitivity_variable: any, sensitivity_model_part: Kratos.ModelPart):
        raise NotImplementedError("Calling base class ResponseFunction::CalculateSensitivity method. Please implement in the derrived class.")

    def GetModelPart(self) -> Kratos.ModelPart:
        raise NotImplementedError("Calling base class ResponseFunction::GetModelPart method. Please implement in the derrived class.")

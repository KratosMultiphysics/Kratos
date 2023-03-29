import KratosMultiphysics as Kratos

class ExecutionPolicy(Kratos.Process):
    def __init__(self):
        super().__init__()

    def GetAnalysisModelPart(self):
        raise NotImplementedError("Calling base class ExecutionPolicy::GetAnalysisModelPart method. This should be implemented in the derived class.")

    def Execute(self):
        raise NotImplementedError("Calling base class ExecutionPolicy::Execute method. This should be implemented in the derived class.")
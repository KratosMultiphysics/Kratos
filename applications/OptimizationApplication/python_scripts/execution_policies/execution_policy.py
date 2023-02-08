import KratosMultiphysics as Kratos

class ExecutionPolicy(Kratos.Process):
    def __init__(self):
        super().__init__()

    def Execute(self):
        raise NotImplementedError("Calling base class ExecutionPolicy::Execute method. This should be implemented in the derrived class.")
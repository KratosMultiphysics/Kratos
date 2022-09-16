# TODO: This class is a placeholder. It should be removed once there are real operations
# in Kratos.

import KratosMultiphysics as KMP
import KratosMultiphysics.python_operation as python_operation

class LocalOperation(python_operation.PythonOperation):
    def __init__(self):
        super().__init__()

    def Execute(self):
        KMP.Logger.PrintWarning("LocalOperation", "Calling the fake LocalOperation Execute().")

def Create(settings, model):
    return LocalOperation()
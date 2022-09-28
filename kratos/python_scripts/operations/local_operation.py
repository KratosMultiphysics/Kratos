# TODO: This class is a placeholder. It should be removed once there are real operations
# in Kratos.

import KratosMultiphysics as KMP

class LocalOperation(KMP.Operation):
    def __init__(self, model, settings):
        super().__init__()

    def Execute(self):
        KMP.Logger.PrintWarning("LocalOperation", "Calling the fake LocalOperation Execute().")

    @staticmethod
    def Create(model,settings):
        return LocalOperation(model,settings)
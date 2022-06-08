import KratosMultiphysics
from KratosMultiphysics import Operation

class FooOperation(KratosMultiphysics.Operation):
    def Execute(self):
        print("Patata")

def Factory(settings, model):
    return FooOperation().Create(model,settings)

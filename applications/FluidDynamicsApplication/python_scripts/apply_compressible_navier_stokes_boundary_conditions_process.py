import KratosMultiphysics
from KratosMultiphysics.FluidDynamicsApplication import ApplyMachDependentBoundaryConditionsTest

def Factory(settings, Model):
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplyMachDependentBoundaryConditionsTest(Model, settings["Parameters"])

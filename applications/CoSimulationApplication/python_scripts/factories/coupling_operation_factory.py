from KratosMultiphysics.CoSimulationApplication.factories import base_factory

def CreateCouplingOperation(coupling_operation_settings, *args):
    """This function creates and returns the Coupling Operation used for CoSimulation"""
    return base_factory.Create(coupling_operation_settings, [*args], "KratosMultiphysics.CoSimulationApplication.coupling_operations")

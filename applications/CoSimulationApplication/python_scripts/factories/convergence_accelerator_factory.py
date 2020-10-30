from KratosMultiphysics.CoSimulationApplication.factories import base_factory

def CreateConvergenceAccelerator(convergence_accelerator_settings):
    """This function creates and returns the Convergence Accelerator used for CoSimulation"""
    return base_factory.Create(convergence_accelerator_settings, [], "KratosMultiphysics.CoSimulationApplication.convergence_accelerators")

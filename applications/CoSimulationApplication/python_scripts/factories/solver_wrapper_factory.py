from KratosMultiphysics.CoSimulationApplication.factories import base_factory

def CreateSolverWrapper(settings, solver_name):
    """This function creates and returns the Wrapper for the Solver used for CoSimulation"""
    return base_factory.Create(settings, [solver_name], "KratosMultiphysics.CoSimulationApplication")

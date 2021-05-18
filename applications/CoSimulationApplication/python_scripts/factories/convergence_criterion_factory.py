from KratosMultiphysics.CoSimulationApplication.factories import base_factory

def CreateConvergenceCriterion(convergence_criterion_settings, *args):
    """This function creates and returns the Convergence Criterion used for CoSimulation"""
    return base_factory.Create(convergence_criterion_settings, [*args], "KratosMultiphysics.CoSimulationApplication.convergence_criteria")

from KratosMultiphysics.CoSimulationApplication.factories import base_factory

def CreateConvergenceCriterion(convergence_criterion_settings):
    """This function creates and returns the Convergence Criterion used for CoSimulation"""
    return base_factory.Create(convergence_criterion_settings, [], "KratosMultiphysics.CoSimulationApplication.convergence_criteria")

def CreateConvergenceCriterionWithoutWrapper(convergence_criterion_settings, solvers_wrappers):
    """This function creates and returns the Convergence Criterion used for CoSimulation"""

    convergence_criteria = base_factory.Create(convergence_criterion_settings, [], "KratosMultiphysics.CoSimulationApplication.convergence_criteria")
    # add access functions to base convergence class since the solves can't be passed through the factory.
    convergence_criteria.SetSolvers(convergence_criterion_settings,solvers_wrappers)
    convergence_criteria.SetInterfaceData(convergence_criterion_settings, solvers_wrappers)

    return convergence_criteria
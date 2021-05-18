from KratosMultiphysics.CoSimulationApplication.factories import base_factory

def CreatePredictor(predictor_settings, solver_wrapper):
    """This function creates and returns the Predictor used for CoSimulation"""
    return base_factory.Create(predictor_settings, [solver_wrapper], "KratosMultiphysics.CoSimulationApplication.predictors")

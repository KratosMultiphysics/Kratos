from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7


def CreatePredictor(settings, solver):
    """
    This function creates and returns the predictor.
    """
    predictor_type = settings["type"]
    module_name = predictor_type.GetString()
    module_full = "KratosMultiphysics.CoSimulationApplication.custom_predictors."+module_name+"_predictor"
    predictor_module = __import__(module_full, fromlist=[module_name])
    return predictor_module.Create(settings, solver)

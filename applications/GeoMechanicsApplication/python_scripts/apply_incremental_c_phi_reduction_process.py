import KratosMultiphysics as Core
import KratosMultiphysics.GeoMechanicsApplication as Geo

# TODO: Add adaptive and other increment strategies in the future.
def Factory(settings, model):
    if not isinstance(settings, Core.Parameters):
        raise TypeError("expected input shall be a Parameters object, encapsulating a json string")

    parameters = settings["Parameters"]

    default_parameters = Core.Parameters("""
    {
        "model_part_name": "",
        "increment_strategy": "fixed",
        "factor_increment": 0.1,
        "max_trials": 50
    }
    """)

    parameters.ValidateAndAssignDefaults(default_parameters)

    return Geo.ApplyIncrementalCPhiReductionProcess(model, parameters)

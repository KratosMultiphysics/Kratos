# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.GeoMechanicsApplication as KGM


def Factory(settings, model):
    """
    This process sets multiple moving load conditions. The 'load' is to be filled in in x,y and z direction. The 'direction'
    parameter indicates the direction of the movement of the load in x,y and z direction, this parameter is either a
    positive or a negative integer; note that the load follows a given line, thus the 'direction' parameter is not a
    normalised direction vector. The 'velocity' parameter indicates the velocity of the load in the given direction,
    this parameter can be either a double or a function of time, written as a string. The 'origin' parameter indicates
    the origin point of the moving load, note that this point needs to be located within the line condition. The configuration
    term provides the offset distance offset along the moving load line condition for each moving point load
    """
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise RuntimeError("Expected input shall be a Parameters object, encapsulating a json string")

    default_settings = KratosMultiphysics.Parameters("""
            {
                "help"                    : "This process applies a moving load condition belonging to a modelpart. The load moves over line elements.",
                "model_part_name"         : "please_specify_model_part_name",
                "compute_model_part_name" : "please_specify_model_part_name",
                "variable_name"           : "POINT_LOAD",
                "load"                    : [0.0, 1.0, 0.0],
                "direction"               : [1,1,1],
                "velocity"                : 1,
                "origin"                  : [0.0,0.0,0.0],
                "configuration"           : [0.0]
            }
            """
                                                     )
    load_settings = settings["Parameters"]
    load_settings.ValidateAndAssignDefaults(default_settings)

    # Set process
    model_part = model.GetModelPart(load_settings["model_part_name"].GetString())
    return KGM.SetMultipleMovingLoadsProcess(model_part, load_settings)


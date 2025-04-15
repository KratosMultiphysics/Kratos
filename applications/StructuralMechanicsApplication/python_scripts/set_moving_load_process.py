# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.StructuralMechanicsApplication as KSM


def Factory(settings, Model):
    """
    This process sets the moving load condition. The 'load' is to be filled in in x,y and z direction. The 'direction'
    parameter indicates the direction of the movement of the load in x,y and z direction, this parameter is either a
    positive or a negative integer; note that the load follows a given line, thus the 'direction' parameter is not a
    normalised direction vector. The 'velocity' parameter indicates the velocity of the load in the given direction,
    this parameter can be either a double or a function of time, written as a string. The 'origin' parameter indicates
    the origin point of the moving load, note that this point needs to be located within the line condition.
    """
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")

    default_settings = KratosMultiphysics.Parameters("""
            {
                "help"            : "This process applies a moving load condition belonging to a modelpart. The load moves over line elements.",
                "model_part_name" : "please_specify_model_part_name",
                "variable_name"   : "POINT_LOAD",
                "load"            : [0.0, 1.0, 0.0],
                "direction"       : [1,1,1],
                "velocity"        : 1,
                "origin"          : [0.0,0.0,0.0],
                "offset"          : 0.0
            }
            """
                                                     )
    load_settings = settings["Parameters"]
    load_settings.ValidateAndAssignDefaults(default_settings)

    # Set process
    model_part = Model.GetModelPart(load_settings["model_part_name"].GetString())
    return KSM.SetMovingLoadProcess(model_part, load_settings)

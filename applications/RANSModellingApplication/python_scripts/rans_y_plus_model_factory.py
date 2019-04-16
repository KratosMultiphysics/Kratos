import KratosMultiphysics
import KratosMultiphysics.RANSModellingApplication as KratosRANS

def Factory(model_part, settings):
    default_settings = KratosMultiphysics.Parameters(r'''
    {
        "model_type"     : "logarithmic",
        "model_settings" : {}
    }''')

    settings.ValidateAndAssignDefaults(default_settings)

    if settings["model_type"].GetString() == "logarithmic":
        return KratosRANS.RansLogarithmicYPlusModelProcess(model_part, settings["model_settings"])
    else:
        raise Exception("Unknown y_plus model_type: " + settings["model_type"].GetString())
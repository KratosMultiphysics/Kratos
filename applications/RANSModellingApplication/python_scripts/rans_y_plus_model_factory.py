import KratosMultiphysics
import KratosMultiphysics.RANSModellingApplication as KratosRANS

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    if(type(Model) != KratosMultiphysics.Model):
        raise Exception("expected input shall be a Model object")

    default_settings = KratosMultiphysics.Parameters(r'''
    {
        "model_type"      : "logarithmic",
        "model_part_name" : "PLEASE_CHOOSE_MODEL_PART_NAME",
        "model_settings"  : {}
    }''')

    settings.ValidateAndAssignDefaults(default_settings)
    model_part = Model[settings["model_part_name"].GetString()]
    y_plus_model_name = settings["model_type"].GetString()

    y_plus_model_list = ["logarithmic", "turbulent_kinetic_energy_based"]

    if (not y_plus_model_name in y_plus_model_list):
        raise Exception("Unknown y_plus \"model_type\" name: \"" + y_plus_model_name + "\".\nAllowed \"model_type\" names: " + y_plus_model_list.__str__())

    if y_plus_model_name == "logarithmic":
        return KratosRANS.RansLogarithmicYPlusModelProcess(model_part, settings["model_settings"])
    elif y_plus_model_name == "turbulent_kinetic_energy_based":
        return KratosRANS.RansTKEYPlusModelProcess(model_part, settings["model_settings"])

def InitializeModelPartName(settings, Model, model_part):
    if (not settings.Has("model_part_name")):
        settings.AddEmptyValue("model_part_name")
        settings["model_part_name"].SetString(model_part.Name)
    else:
        model_part_name = settings["model_part_name"].GetString()
        if model_part_name != model_part.Name:
            KratosMultiphysics.Logger.PrintWarning("RANSYPlusModelSettings", "Using a different model part (i.e. \"" + model_part_name + "\") over the intended model part \"" + model_part.Name + "\"")

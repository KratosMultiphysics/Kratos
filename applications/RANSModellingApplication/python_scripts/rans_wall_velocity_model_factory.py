import KratosMultiphysics
import KratosMultiphysics.RANSModellingApplication as KratosRANS


def Factory(settings, Model):
    if (type(settings) != KratosMultiphysics.Parameters):
        raise Exception(
            "expected input shall be a Parameters object, encapsulating a json string"
        )
    if (type(Model) != KratosMultiphysics.Model):
        raise Exception("expected input shall be a Model object")

    default_settings = KratosMultiphysics.Parameters(r'''
    {
        "model_type"      : "averaged_wall_velocity",
        "model_part_name" : "PLEASE_CHOOSE_MODEL_PART_NAME",
        "model_settings"  : {}
    }''')

    settings.ValidateAndAssignDefaults(default_settings)
    model_part = Model[settings["model_part_name"].GetString()]
    wall_velocity_model_name = settings["model_type"].GetString()

    wall_velocity_model_list = ["averaged_wall_velocity"]

    if (not wall_velocity_model_name in wall_velocity_model_list):
        raise Exception("Unknown wall_velocity \"model_type\" name: \"" +
                        wall_velocity_model_name +
                        "\".\nAllowed \"model_type\" names: " +
                        wall_velocity_model_list.__str__())

    if wall_velocity_model_name == "averaged_wall_velocity":
        return KratosRANS.RansWallVelocityCalculationProcess(
            model_part, settings["model_settings"])


def InitializeModelPartName(settings, Model, model_part):
    if (not settings.Has("model_part_name")):
        settings.AddEmptyValue("model_part_name")
        settings["model_part_name"].SetString(model_part.Name)
    else:
        model_part_name = settings["model_part_name"].GetString()
        if model_part_name != model_part.Name:
            KratosMultiphysics.Logger.PrintWarning(
                "RANSWallVelocityCalculationProcess",
                "Using a different model part (i.e. \"" + model_part_name +
                "\") over the intended model part \"" + model_part.Name + "\"")

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
        "process_name"      : "PLEASE_CHOOSE_PROCESS_NAME",
        "model_part_name"   : "PLEASE_CHOOSE_MODEL_PART_NAME",
        "process_settings"  : {}
    }''')

    settings.ValidateAndAssignDefaults(default_settings)
    model_part = Model[settings["model_part_name"].GetString()]
    process_name = settings["process_name"].GetString()

    process_list = [
        "nut_k_wall_function", "epsilon_wall_function",
        "scalar_value_neighbour_based_averaging", "check_scalar_bounds",
        "y_plus_logarithmic", "y_plus_based_wall_distance",
        "apply_constant_scalar_value", "scalar_cell_center_averaging"
    ]

    if (not process_name in process_list):
        raise Exception("Unknown auxiliary process_name: \"" + process_name +
                        "\".\nAllowed \"process_name\" names: " +
                        process_list.__str__())

    if process_name == "nut_k_wall_function":
        return KratosRANS.RansNutKWallFunctionProcess(
            model_part, settings["process_settings"])
    elif (process_name == "epsilon_wall_function"):
        return KratosRANS.RansEpsilonWallFunctionProcess(
            model_part, settings["process_settings"])
    elif (process_name == "scalar_value_neighbour_based_averaging"):
        return KratosRANS.RansScalarNeighbourAveragingProcess(
            model_part, settings["process_settings"])
    elif (process_name == "check_scalar_bounds"):
        return KratosRANS.RansCheckScalarBoundsProcess(
            model_part, settings["process_settings"])
    elif (process_name == "y_plus_logarithmic"):
        return KratosRANS.RansLogarithmicYPlusModelProcess(
            model_part, settings["process_settings"])
    elif (process_name == "y_plus_based_wall_distance"):
        return KratosRANS.RansYPlusWallDistanceCalculationProcess(
            model_part, settings["process_settings"])
    elif (process_name == "apply_constant_scalar_value"):
        return KratosMultiphysics.ApplyConstantScalarValueProcess(
            model_part, settings["process_settings"])
    elif (process_name == "scalar_cell_center_averaging"):
        return KratosRANS.RansScalarCellCenterAveragingProcess(
            model_part, settings["process_settings"])

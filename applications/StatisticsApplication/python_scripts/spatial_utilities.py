import KratosMultiphysics as Kratos
import KratosMultiphysics.StatisticsApplication as Statistics

def CheckDefaultVariableSettings(settings):
    default_parameters = Kratos.Parameters("""
    {
        "variable_container" : "nodal_historical",
        "variable_names"     : []
    }  """)

    settings.ValidateAndAssignDefaults(default_parameters)


def CheckVariables(settings, model_part):
    for variable_name in settings["variable_names"].GetStringArray():
        if (not Kratos.KratosGlobals.HasVariable(variable_name)):
            raise Exception(variable_name.GetString() + " not found.")

    if (settings["variable_container"].GetString() == "nodal_historical"):
        for variable_name in settings["variable_names"].GetStringArray():
            variable = Kratos.KratosGlobals.GetVariable(variable_name)
            if (not model_part.HasNodalSolutionStepVariable(variable)):
                raise Exception("Solution step variable " + variable.Name +
                                " not found in " + model_part.Name)


def GetMethodsContainer(settings):
    variable_container_type = settings["variable_container"].GetString()
    if (variable_container_type == "nodal_historical"):
        return Statistics.SpatialMethods.Historical
    elif (variable_container_type == "nodal_non_historical"):
        return Statistics.SpatialMethods.NonHistorical.Nodes
    elif (variable_container_type == "element_non_historical"):
        return Statistics.SpatialMethods.NonHistorical.Elements
    elif (variable_container_type == "condition_non_historical"):
        return Statistics.SpatialMethods.NonHistorical.Conditions
    else:
        msg = "Unknown variable container type [ variable_container_type = " + variable_container_type + "\n"
        msg += "   Allowed variable container types are:\n"
        msg += "       nodal_historical\n"
        msg += "       nodal_non_historical\n"
        msg += "       element_non_historical\n"
        msg += "       condition_non_historical\n"
        raise Exception(msg)

class SpatialOutput:
    def __init__(self, column_headers):
        Kratos.Logger.PrintInfo("StatisticsApplication", "Spatial statistics output:")
        msg = ""
        for header in column_headers:
            msg += "     " + header
        Kratos.Logger.PrintInfo("StatisticsApplication", msg)

    def PrintOutput(self, value_array):
        msg = ""
        for value in value_array:
            msg += "    " + value
        Kratos.Logger.PrintInfo("StatisticsApplication", msg)
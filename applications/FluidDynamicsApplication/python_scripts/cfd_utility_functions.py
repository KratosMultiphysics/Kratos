# Importing the Kratos Library
import KratosMultiphysics as Kratos
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD


def CheckVariableType(variable_name, required_variable_type):
    variable_type = Kratos.KratosGlobals.GetVariableType(variable_name)
    if (variable_type != required_variable_type):
        msg = "Type of " + variable_name + " is not supported. Please provide a " + required_variable_type + " [ " + variable_name + " = " + variable_type + " ]."
        raise Exception(msg)


def GetCFDUtilityFunction(params):
    """This is the factory for CFDUtilityFunctions

    This method returns the correct lambda function from CFDUtilityFunctions class

    Args:
        params (Kratos.Parameters): Parameters required for function execution

    Returns:
        [lambda model_part]: Returns correct CFDUtilityFunction lambda
    """
    if (params.Has("function_name")):
        function_name = params["function_name"].GetString()
    else:
        raise Exception("Please provide a function_name.")

    function_list = [func for func in dir(CFDUtilityFunctions)]

    if (function_name not in function_list):
        msg = "Unknown function name [ \"function_name\" = \"" + function_name + "\" ].\n"
        msg += "Supported function names are:"
        for func in function_list:
            msg += "\n\t" + func
        raise Exception(msg)

    Kratos.Logger.PrintInfo("CFDUtilityFunctions",
                            "Created " + function_name + " function.")
    return getattr(CFDUtilityFunctions, function_name)(params)


class CFDUtilityFunctions:
    """A class to hold CFDUtilities C++ methods

    This class holds wrappers for CFDUtilities developed in C++. These C++ methods
    should be standalone, and they should be able to have single point execution without
    depending on any other information rather than model_part which is being passed as an
    input argument at execution point.

    This wrapper methods are used to provide additional information for CFDUtilities methods
    based on input from Kratos.Parameters.

    The wrappers defined here needs to have following signature:
        Input args : (Kratos.Parameters)
        Output args: Returns a lambda function. The lambda function should have Kratos.ModelPart as
                     an input, and not output should be there.

    """
    @staticmethod
    def DistributeConditionVariableToNodes(params):
        default_settings = Kratos.Parameters("""
        {
            "function_name": "DistributeConditionVariableToNodes",
            "variable_name": "PLEASE_PROVIDE_A_VARIABLE_NAME"
        }""")

        params.ValidateAndAssignDefaults(default_settings)
        variable_name = params["variable_name"].GetString()
        if (not Kratos.KratosGlobals.HasVariable(variable_name)):
            raise Exception("Variable " + variable_name + " not found.")
        variable_type = Kratos.KratosGlobals.GetVariableType(variable_name)
        variable = Kratos.KratosGlobals.GetVariable(variable_name)
        if (variable_type in ["Double", "Array"]):
            return lambda model_part: KratosCFD.CFDUtilities.DistributeConditionVariableToNodes(
                model_part, variable)
        else:
            raise Exception("Unsupported variable type " + variable_type +
                            " in " + variable_name + " variable.")

    @staticmethod
    def CalculateYPlusAndUTauForConditionsBasedOnReaction(params):
        default_settings = Kratos.Parameters("""
        {
            "function_name"                    : "CalculateYPlusAndUTauForConditionsBasedOnReaction",
            "kinematic_viscosity_variable_name": "VISCOSITY",
            "reaction_variable_name"           : "REACTION"
        }""")

        params.ValidateAndAssignDefaults(default_settings)

        nu_variable_name = params[
            "kinematic_viscosity_variable_name"].GetString()
        CheckVariableType(nu_variable_name, "Double")
        nu_variable = Kratos.KratosGlobals.GetVariable(nu_variable_name)

        reaction_variable_name = params["reaction_variable_name"].GetString()
        CheckVariableType(reaction_variable_name, "Array")
        reaction_variable = Kratos.KratosGlobals.GetVariable(
            reaction_variable_name)

        return lambda model_part: KratosCFD.CFDUtilities.CalculateYPlusAndUTauForConditionsBasedOnReaction(
            model_part, nu_variable, reaction_variable)

    @staticmethod
    def CalculateYPlusAndUTauForConditionsBasedOnLinearLogarithmicWallFunction(
        params):
        default_settings = Kratos.Parameters("""
        {
            "function_name"                    : "CalculateYPlusAndUTauForConditionsBasedOnLinearLogarithmicWallFunction",
            "kinematic_viscosity_variable_name": "VISCOSITY",
            "von_karman"                       : 0.41,
            "wall_smoothness"                  : 5.2,
            "max_iterations"                   : 20,
            "tolerance"                        : 1e-6
        }""")

        params.ValidateAndAssignDefaults(default_settings)

        nu_variable_name = params[
            "kinematic_viscosity_variable_name"].GetString()
        CheckVariableType(nu_variable_name, "Double")
        nu_variable = Kratos.KratosGlobals.GetVariable(nu_variable_name)

        kappa = params["von_karman"].GetDouble()
        beta = params["wall_smoothness"].GetDouble()
        max_iterations = params["max_iterations"].GetInt()
        tolerance = params["tolerance"].GetDouble()

        return lambda model_part: KratosCFD.CFDUtilities.CalculateYPlusAndUTauForConditionsBasedOnLinearLogarithmicWallFunction(
            model_part, nu_variable, kappa, beta, max_iterations, tolerance)

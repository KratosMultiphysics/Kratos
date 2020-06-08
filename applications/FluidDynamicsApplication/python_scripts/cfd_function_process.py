# Importing the Kratos Library
import KratosMultiphysics as Kratos
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD


class CFDFunctions:
    @staticmethod
    def GetFunction(params):
        if (params.Has("function_name")):
            function_name = "Function" + params["function_name"].GetString()
        else:
            raise Exception("Please provide a function_name.")

        function_list = [
            func for func in dir(CFDFunctions) if func.startswith("Function")
        ]

        if (function_name not in function_list):
            msg = "Unknown function name [ \"function_name\" = \"" + function_name[
                8:] + "\" ].\n"
            msg += "Supported function names are:"
            for func in function_list:
                msg += "\n\t" + func[8:]
            raise Exception(msg)

        Kratos.Logger.PrintInfo("CFDFunctions",
                                "Created " + function_name[8:] + " function.")
        return getattr(CFDFunctions, function_name)(params)

    @staticmethod
    def FunctionDistributeConditionDataToNodes(params):
        default_settings = Kratos.Parameters("""
        {
            "function_name": "DistributeConditionDataToNodes",
            "variable_name": "PLEASE_PROVIDE_A_VARIABLE_NAME"
        }""")

        params.ValidateAndAssignDefaults(default_settings)
        variable_name = params["variable_name"].GetString()
        if (not Kratos.KratosGlobals.HasVariable(variable_name)):
            raise Exception("Variable " + variable_name + " not found.")
        variable_type = Kratos.KratosGlobals.GetVariableType(variable_name)
        variable = Kratos.KratosGlobals.GetVariable(variable_name)
        if (variable_type in ["Double", "Array"]):
            return lambda model_part: KratosCFD.DistributeConditionDataToNodes(
                model_part, variable)
        else:
            raise Exception("Unsupported variable type " + variable_type +
                            " in " + variable_name + " variable.")

    @staticmethod
    def FunctionCalculateYPlusAndUTauForConditionsBasedOnReaction(params):
        default_settings = Kratos.Parameters("""
        {
            "function_name"                    : "CalculateYPlusAndUTauForConditionsBasedOnReaction",
            "kinematic_viscosity_variable_name": "VISCOSITY",
            "reaction_variable_name"           : "REACTION"
        }""")

        params.ValidateAndAssignDefaults(default_settings)

        nu_variable_name = params[
            "kinematic_viscosity_variable_name"].GetString()
        CFDFunctions.__CheckVariableType(nu_variable_name, "Double")
        nu_variable = Kratos.KratosGlobals.GetVariable(nu_variable_name)

        reaction_variable_name = params["reaction_variable_name"].GetString()
        CFDFunctions.__CheckVariableType(reaction_variable_name, "Array")
        reaction_variable = Kratos.KratosGlobals.GetVariable(
            reaction_variable_name)

        return lambda model_part: KratosCFD.CFDUtilities.CalculateYPlusAndUTauForConditionsBasedOnReaction(
            model_part, nu_variable, reaction_variable)

    @staticmethod
    def FunctionCalculateYPlusAndUTauForConditionsBasedOnLinearLogarithmicWallFunction(
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
        CFDFunctions.__CheckVariableType(nu_variable_name, "Double")
        nu_variable = Kratos.KratosGlobals.GetVariable(nu_variable_name)

        kappa = params["von_karman"].GetDouble()
        beta = params["wall_smoothness"].GetDouble()
        max_iterations = params["max_iterations"].GetInt()
        tolerance = params["tolerance"].GetDouble()

        return lambda model_part: KratosCFD.CFDUtilities.CalculateYPlusAndUTauForConditionsBasedOnLinearLogarithmicWallFunction(
            model_part, nu_variable, kappa, beta, max_iterations, tolerance)

    @staticmethod
    def __CheckVariableType(variable_name, required_variable_type):
        variable_type = Kratos.KratosGlobals.GetVariableType(variable_name)
        if (variable_type != required_variable_type):
            msg = "Type of " + variable_name + " is not supported. Please provide a " + required_variable_type + " [ " + variable_name + " = " + variable_type + " ]."
            raise Exception(msg)


def Factory(settings, model):
    if (type(settings) != Kratos.Parameters):
        raise Exception(
            "expected input shall be a Parameters object, encapsulating a json string"
        )
    if (type(model) != Kratos.Model):
        raise Exception(
            "expected input shall be a Model object, encapsulating a json string"
        )

    return CFDFunctionProcess(model, settings["Parameters"])


class CFDFunctionProcess(Kratos.Process):
    def __init__(self, model, params):
        Kratos.Process.__init__(self)

        default_settings = Kratos.Parameters("""
            {
                "model_part_name"  : "",
                "function_settings": {},
                "execution_points" : ["FinalizeSolutionStep"]
            }
            """)

        params.ValidateAndAssignDefaults(default_settings)
        self.model_part = model.GetModelPart(
            params["model_part_name"].GetString())
        self.execution_points_list = params["execution_points"].GetStringArray(
        )

        self.function = CFDFunctions.GetFunction(params["function_settings"])

    def ExecuteInitialize(self):
        if ("Initialize" in self.execution_points_list):
            self.function(self.model_part)

    def ExecuteInitializeSolutionStep(self):
        if ("InitializeSolutionStep" in self.execution_points_list):
            self.function(self.model_part)

    def ExecuteFinalizeSolutionStep(self):
        if ("FinalizeSolutionStep" in self.execution_points_list):
            self.function(self.model_part)

    def ExecuteFinalize(self):
        if ("Finalize" in self.execution_points_list):
            self.function(self.model_part)
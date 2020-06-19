# Importing the Kratos Library
import KratosMultiphysics as Kratos
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD

from .cfd_utility_functions import GetCFDUtilityFunction


def Factory(settings, model):
    if (not isinstance(settings, Kratos.Parameters)):
        raise Exception(
            "expected input shall be a Parameters object, encapsulating a json string"
        )
    if (not isinstance(model, Kratos.Model)):
        raise Exception(
            "expected input shall be a Model object, encapsulating a json string"
        )

    return CFDUtilityFunctionProcess(model, settings["Parameters"])


class CFDUtilityFunctionProcess(Kratos.Process):
    """This process executes CFD functions implemented in c++ level

    This process is used to evaluate CFD functions implemented in c++ level
    and exposed through "CFDUtilities" submodule in FluidDynamicsApplication.

    Kratos Parameter settings:
        "model_part_name"  : The model part on which CFDUtility function will be executed
        "function_settings": Settings of appropriate function. ("function_name" is a must have under this group)
        "execution_points" : User defined execution points in simulation, where function is executed.
                                Allowed execution points:
                                    "Initialize",
                                    "InitializeSolutionStep",
                                    "FinalizeSolutionStep",
                                    "Finalize"

    """
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

        allowed_execution_points = [
            "Initialize", "InitializeSolutionStep", "FinalizeSolutionStep",
            "Finalize"
        ]
        for execution_point in self.execution_points_list:
            if (execution_point not in allowed_execution_points):
                msg = "Unknown function execution point [ \"execution_point\" = \"" + execution_point + "\" ].\n"
                msg += "Supported execution points are:"
                msg += "\n\t".join(allowed_execution_points)
                raise Exception(msg)

        self.function = GetCFDUtilityFunction(params["function_settings"])

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
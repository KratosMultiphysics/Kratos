from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# Check that applications were imported in the main script
KratosMultiphysics.CheckRegisteredApplications("ShapeOptimizationApplication")

def Factory(settings, model):
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    return FixVectorVariableProcess(model, settings["Parameters"])

class FixVectorVariableProcess(KratosMultiphysics.Process):
    """ This process fixes the selected components of a given vector variable
    without modifying the value of the variable.
    This works for vector variables only"""
    def __init__(self, model, settings):
        KratosMultiphysics.Process.__init__(self)

        default_settings = KratosMultiphysics.Parameters(
            """
            {
                "model_part_name"      : "SPECIFY_MODEL_PART_NAME",
                "variable_name"        : "SPECIFY_VARIABLE_NAME",
                "constrained"          : [true,true,true]
            }
            """
        )

        settings.RecursivelyValidateAndAssignDefaults(default_settings)
        self.settings = settings

        self.model_part = model[self.settings["model_part_name"].GetString()]

        self.variable = KratosMultiphysics.KratosGlobals.GetVariable(self.settings["variable_name"].GetString())

        if not isinstance(self.variable, KratosMultiphysics.Array1DVariable3) and not isinstance(self.variable, KratosMultiphysics.VectorVariable):
            msg = "Error in FixVectorVariableProcess. Variable type of variable : " + settings["variable_name"].GetString() + " is incorrect. Must be a vector or array3"
            raise Exception(msg)

        self.variable_components_fixed = {}

        components = ["_X", "_Y", "_Z"]
        for i, component in enumerate(components):
            variable = KratosMultiphysics.KratosGlobals.GetVariable(self.settings["variable_name"].GetString()+component)
            self.variable_components_fixed[variable] = self.settings["constrained"][i].GetBool()

        self.variable_utils = KratosMultiphysics.VariableUtils()

    def ExecuteInitialize(self):
        pass

    def ExecuteBeforeSolutionLoop(self):
        pass

    def ExecuteInitializeSolutionStep(self):
        for variable, is_fixed in self.variable_components_fixed.items():
            if is_fixed:
                self.variable_utils.ApplyFixity(variable, is_fixed, self.model_part.Nodes)

    def ExecuteFinalizeSolutionStep(self):
        pass

    def ExecuteBeforeOutputStep(self):
        pass

    def ExecuteAfterOutputStep(self):
        pass

    def ExecuteFinalize(self):
        #here we free all of the nodes in the mesh
        for variable, is_fixed in self.variable_components_fixed.items():
            if is_fixed:
                fixity_status = False
                self.variable_utils.ApplyFixity(variable, fixity_status, self.model_part.Nodes)



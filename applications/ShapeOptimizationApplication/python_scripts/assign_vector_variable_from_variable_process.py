from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# Check that applications were imported in the main script
KratosMultiphysics.CheckRegisteredApplications("ShapeOptimizationApplication")

# Import applications
import KratosMultiphysics.ShapeOptimizationApplication as ShapeOptimizationApplication


def Factory(settings, model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    return AssignVectorVariableFromVariableProcess(model, settings["Parameters"])

class AssignVectorVariableFromVariableProcess(KratosMultiphysics.Process):
    """ This process copies the source variable onto the variable and
    fixes the constrained components of the variable.
    This works for vector variables only"""
    def __init__(self, model, settings):
        KratosMultiphysics.Process.__init__(self)

        default_settings = KratosMultiphysics.Parameters(
            """
            {
                "model_part_name"      : "SPECIFY_MODEL_PART_NAME",
                "variable_name"        : "SPECIFY_VARIABLE_NAME",
                "source_variable_name" : "SPECIFY_SOURCE_VARIABLE_NAME",
                "constrained"          : [true,true,true]
            }
            """
        )

        settings.RecursivelyValidateAndAssignDefaults(default_settings)
        self.settings = settings

        self.model_part = model[self.settings["model_part_name"].GetString()]

        self.variable = KratosMultiphysics.KratosGlobals.GetVariable(self.settings["variable_name"].GetString())
        self.source_variable = KratosMultiphysics.KratosGlobals.GetVariable(self.settings["source_variable_name"].GetString())

        if(type(self.variable) != KratosMultiphysics.Array1DVariable3 and type(self.variable) != KratosMultiphysics.VectorVariable):
            msg = "Error in AssignVectorVariableFromVariableProcess. Variable type of variable : " + settings["variable_name"].GetString() + " is incorrect . Must be a vector or array3"
            raise Exception(msg)
        if type(self.variable) != type(self.source_variable):
            raise Exception("Error in AssignVectorVariableFromVariableProcess: Variable and source variable do not have the same type!")

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

        self.variable_utils.CopyVectorVar(self.source_variable, self.variable, self.model_part.Nodes)

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



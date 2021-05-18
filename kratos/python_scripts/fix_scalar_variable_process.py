# Importing the Kratos Library
import KratosMultiphysics

def Factory(settings, model):
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    return FixScalarVariableProcess(model, settings["Parameters"])

class FixScalarVariableProcess(KratosMultiphysics.Process):
    """ This process fixes a scalar variable without modifying the value of the variable.
    This works for scalar or component variables only"""
    def __init__(self, model, settings):
        KratosMultiphysics.Process.__init__(self)

        default_settings = KratosMultiphysics.Parameters(
            """
            {
                "model_part_name" : "SPECIFY_MODEL_PART_NAME",
                "variable_name"   : "SPECIFY_VARIABLE_NAME",
                "interval"        : [0.0, 1e30],
                "constrained"     : true
            }
            """
        )

        #assign this here since it will change the "interval" prior to validation
        self.interval = KratosMultiphysics.IntervalUtility(settings)

        settings.ValidateAndAssignDefaults(default_settings)
        self.settings = settings

        self.model_part = model[self.settings["model_part_name"].GetString()]

        self.variable = KratosMultiphysics.KratosGlobals.GetVariable(self.settings["variable_name"].GetString())

        if not isinstance(self.variable, KratosMultiphysics.DoubleVariable):
            msg = "Error in FixScalarVariableProcess. Variable type of variable : " + settings["variable_name"].GetString() + " is incorrect. Must be a scalar or a component"
            raise Exception(msg)

        self.variable_utils = KratosMultiphysics.VariableUtils()

    def ExecuteInitialize(self):
        pass

    def ExecuteBeforeSolutionLoop(self):
        pass

    def ExecuteInitializeSolutionStep(self):
        current_time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]

        if(self.interval.IsInInterval(current_time)):
            is_fixed = self.settings["constrained"].GetBool()
            if is_fixed:
                self.variable_utils.ApplyFixity(self.variable, is_fixed, self.model_part.Nodes)

    def ExecuteFinalizeSolutionStep(self):
        current_time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]

        if(self.interval.IsInInterval(current_time)):
            #here we free all of the nodes in the model part
            is_fixed = self.settings["constrained"].GetBool()
            if is_fixed:
                fixity_status = False
                self.variable_utils.ApplyFixity(self.variable, fixity_status, self.model_part.Nodes)

    def ExecuteBeforeOutputStep(self):
        pass

    def ExecuteAfterOutputStep(self):
        pass

    def ExecuteFinalize(self):
        pass



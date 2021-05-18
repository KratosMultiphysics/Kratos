# Importing the Kratos Library
import KratosMultiphysics

from KratosMultiphysics import fix_scalar_variable_process

def Factory(settings, model):
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    return FixVectorVariableProcess(model, settings["Parameters"])

class FixVectorVariableProcess(KratosMultiphysics.Process):
    """ This process fixes the selected components of a given vector variable
    without modifying the value of the variable.
    Internally a FixScalarVariableProcess is created for each component."""
    def __init__(self, model, settings):
        KratosMultiphysics.Process.__init__(self)

        default_settings = KratosMultiphysics.Parameters(
            """
            {
                "model_part_name" : "SPECIFY_MODEL_PART_NAME",
                "variable_name"   : "SPECIFY_VARIABLE_NAME",
                "interval"        : [0.0, 1e30],
                "constrained"     : [true,true,true]
            }
            """
        )

        settings.ValidateAndAssignDefaults(default_settings)
        self.settings = settings

        variable = KratosMultiphysics.KratosGlobals.GetVariable(self.settings["variable_name"].GetString())

        if not isinstance(variable, KratosMultiphysics.Array1DVariable3) and not isinstance(variable, KratosMultiphysics.VectorVariable):
            msg = "Error in FixVectorVariableProcess. Variable type of variable : " + settings["variable_name"].GetString() + " is incorrect. Must be a vector or array3"
            raise Exception(msg)

        self.aux_processes = []

        # loop over components X, Y and Z
        for index, component in enumerate(["_X", "_Y", "_Z"]):
            i_params = KratosMultiphysics.Parameters("{}")
            i_params.AddValue("model_part_name", settings["model_part_name"])
            i_params.AddEmptyValue("variable_name").SetString(settings["variable_name"].GetString() + component)
            i_params.AddValue("interval", settings["interval"])
            i_params.AddValue("constrained", settings["constrained"][index])
            self.aux_processes.append( fix_scalar_variable_process.FixScalarVariableProcess(model, i_params) )

    def ExecuteInitialize(self):
        for process in self.aux_processes:
            process.ExecuteInitialize()

    def ExecuteBeforeSolutionLoop(self):
        for process in self.aux_processes:
            process.ExecuteBeforeSolutionLoop()

    def ExecuteInitializeSolutionStep(self):
        for process in self.aux_processes:
            process.ExecuteInitializeSolutionStep()

    def ExecuteFinalizeSolutionStep(self):
        for process in self.aux_processes:
            process.ExecuteFinalizeSolutionStep()

    def ExecuteBeforeOutputStep(self):
        for process in self.aux_processes:
            process.ExecuteBeforeOutputStep()

    def ExecuteAfterOutputStep(self):
        for process in self.aux_processes:
            process.ExecuteAfterOutputStep()

    def ExecuteFinalize(self):
        for process in self.aux_processes:
            process.ExecuteFinalize()



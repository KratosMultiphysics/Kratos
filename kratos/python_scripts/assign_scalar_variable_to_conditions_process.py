import KratosMultiphysics
import sys
from math import *

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return AssignScalarVariableToConditionsProcess(Model, settings["Parameters"])

## All the processes python should be derived from "Process"
class AssignScalarVariableToConditionsProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings ):
        KratosMultiphysics.Process.__init__(self)

        default_settings = KratosMultiphysics.Parameters("""
            {
                "help"            : "This process assigns a given value (scalar) to the conditions belonging a certain submodelpart",
                "mesh_id"         : 0,
                "model_part_name" : "please_specify_model_part_name",
                "variable_name"   : "SPECIFY_VARIABLE_NAME",
                "interval"        : [0.0, 1e30],
                "value"           : 0.0,
                "local_axes"      : {}
            }
            """
            )

        # assign this here since it will change the "interval" prior to validation
        self.interval = KratosMultiphysics.IntervalUtility(settings)

        # here i do a trick, since i want to allow "value" to be a string or a double value
        if(settings.Has("value")):
            if(settings["value"].IsString()):
                default_settings["value"].SetString("0.0")

        settings.ValidateAndAssignDefaults(default_settings)

        self.model_part = Model[settings["model_part_name"].GetString()]

        self.value_is_numeric = False

        # set processes
        params = KratosMultiphysics.Parameters("{}")           
        params.AddValue("model_part_name", settings["model_part_name"])
        params.AddValue("mesh_id", settings["mesh_id"])
        params.AddValue("value", settings["value"])
        params.AddValue("variable_name", settings["variable_name"])

        # set processes
        if settings["value"].IsNumber():
            self.value_is_numeric = True
            self.AssignValueProcess = KratosMultiphysics.AssignScalarVariableToConditionsProcess(self.model_part, params)
        else:
            params.AddValue("local_axes", settings["local_axes"])
            self.AssignValueProcess = KratosMultiphysics.AssignScalarFieldToConditionsProcess(self.model_part, params)
                
        # construct a variable_utils object to speedup fixing
        self.step_is_active = False

    def ExecuteInitializeSolutionStep(self):
        current_time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]

        if(self.interval.IsInInterval(current_time)):

            self.step_is_active = True

            self.AssignValueProcess.Execute()


    def ExecuteFinalizeSolutionStep(self):
        current_time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]
        self.step_is_active = False

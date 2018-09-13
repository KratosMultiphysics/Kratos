from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.SolidMechanicsApplication as KratosSolid

## This proces sets the value of a scalar variable to conditions
import assign_modulus_and_direction_to_conditions_process as BaseProcess

def Factory(custom_settings, Model):
    if( not isinstance(custom_settings,KratosMultiphysics.Parameters) ):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return AssignModulusAndDirectionToConditionsProcess(Model, custom_settings["Parameters"])

## All the processes python should be derived from "Process"
class AssignModulusAndDirectionToConditionsProcess(BaseProcess.AssignModulusAndDirectionToConditionsProcess):
    def __init__(self, Model, custom_settings ):
        KratosMultiphysics.Process.__init__(self)

        ##settings string in json format
        default_settings = KratosMultiphysics.Parameters("""
        {
             "help" : "This process assigns a torque to conditions",
             "model_part_name": "MODEL_PART_NAME",
             "variable_name": "VARIABLE_NAME",
             "modulus" : 0.0,
             "direction": [0.0, 0.0, 0.0],
             "center": [0.0, 0.0, 0.0],
             "constrained": false,
             "interval": [0.0, "End"]
        }
        """)

        #trick to allow "value" to be a string or a double value
        if(custom_settings.Has("modulus")):
            if(custom_settings["modulus"].IsString()):
                default_settings["modulus"].SetString("0.0")


        ##overwrite the default settings with user-provided parameters
        self.settings = custom_settings
        self.settings.ValidateAndAssignDefaults(default_settings)

        self.custom_settings = custom_settings

        ###assign scalar process
        params = KratosMultiphysics.Parameters("{}")
        params.AddValue("model_part_name", self.settings["model_part_name"])
        params.AddValue("variable_name",self.settings["variable_name"])
        params.AddValue("modulus", self.settings["modulus"])
        params.AddValue("direction", self.settings["direction"])
        params.AddValue("constrained", self.settings["constrained"])
        params.AddValue("interval",self.settings["interval"])

        BaseProcess.AssignModulusAndDirectionToConditionsProcess.__init__(self, Model, params)


    def ExecuteInitialize(self):

        # set model part
        self.model_part = self.model[self.settings["model_part_name"].GetString()]
        if( self.model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED] == False ):
            self.model_part.ProcessInfo.SetValue(KratosMultiphysics.INTERVAL_END_TIME, self.interval[1])

        # set processes
        params = KratosMultiphysics.Parameters("{}")
        params.AddValue("model_part_name", self.settings["model_part_name"])
        params.AddValue("variable_name", self.settings["variable_name"])
        params.AddValue("modulus", self.custom_settings["modulus"])
        params.AddValue("direction", self.custom_settings["direction"])
        params.AddValue("center", self.custom_settings["center"])

        self.CreateAssignmentProcess(params)

        if( self.IsInsideInterval() and self.interval_string == "initial" ):
            self.AssignValueProcess.Execute()


    #
    def CreateAssignmentProcess(self, params):

        if( self.value_is_numeric ):
            self.AssignValueProcess = KratosSolid.AssignTorqueAboutAnAxisToConditionsProcess(self.model_part, params)
        else:
            self.AssignValueProcess = KratosSolid.AssignTorqueFieldAboutAnAxisToConditionsProcess(self.model_part, self.compiled_function, "function",  self.value_is_spatial_function, params)

from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.SolidMechanicsApplication as KratosSolid

## This proces sets the value of a vector variable by a direction and a modulus.
## In this case, the fixicity is given set by deffault to true.
import math
import sys

import assign_scalar_to_nodes_process as BaseProcess

def Factory(custom_settings, Model):
    if( not isinstance(custom_settings,KratosMultiphysics.Parameters) ):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return AssignRigidBodyRotationToNodesProcess(Model, custom_settings["Parameters"])

## All the processes python should be derived from "Process"
class AssignRigidBodyRotationToNodesProcess(BaseProcess.AssignScalarToNodesProcess):
    def __init__(self, Model, custom_settings ):
        KratosMultiphysics.Process.__init__(self)

        ##settings string in json format
        default_settings = KratosMultiphysics.Parameters("""
        {
             "help" : "This process assigns a rotation to nodes given by a rotation vector, modulus and direction",
             "model_part_name": "MODEL_PART_NAME",
             "variable_name": "VARIABLE_NAME",
             "modulus" : 0.0,
             "direction": [0.0, 0.0, 1.0],
             "center": [0.0, 0.0, 0.0],
             "constrained": true,
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

        ##check if variable type is a vector
        self.var = KratosMultiphysics.KratosGlobals.GetVariable(self.settings["variable_name"].GetString())
        if( not isinstance(self.var,KratosMultiphysics.Array1DVariable3) ):
            raise Exception("Variable type is incorrect. Must be a three-component vector.")


        ###assign scalar process
        params = KratosMultiphysics.Parameters("{}")
        params.AddValue("model_part_name", self.settings["model_part_name"])
        params.AddValue("variable_name",self.settings["variable_name"])
        params.AddValue("value", self.settings["modulus"])
        params.AddValue("constrained", self.settings["constrained"])
        params.AddValue("interval",self.settings["interval"])

        BaseProcess.AssignScalarToNodesProcess.__init__(self, Model, params)

    #
    def CheckVariableType(self,name):

        self.var = KratosMultiphysics.KratosGlobals.GetVariable(name)
        if( not isinstance(self.var,KratosMultiphysics.Array1DVariable3) ):
            raise Exception("Variable type is incorrect. Must be an array_1d vector")

    #
    def SetFixAndFreeProcesses(self,params):

        self.fix_derivated_variable = False
        if( self.fix_derivated_variable == False ):
            for dynamic_variable in self.LinearDynamicVariables:
                if dynamic_variable == self.variable_name:
                    self.derivated_variable_name = "DISPLACEMENT"
                    self.fix_derivated_variable = True
                    self.SetLinearTimeIntegration()
                    break


        if( self.fix_derivated_variable ):

            params["variable_name"].SetString(self.derivated_variable_name+"_X")
            fix_dof_process  =  KratosSolid.FixScalarDofProcess(self.model_part, params)
            self.FixDofsProcesses.append(fix_dof_process)
            free_dof_process = KratosSolid.FreeScalarDofProcess(self.model_part, params)
            self.FreeDofsProcesses.append(free_dof_process)

            params["variable_name"].SetString(self.derivated_variable_name+"_Y")
            fix_dof_process  =  KratosSolid.FixScalarDofProcess(self.model_part, params)
            self.FixDofsProcesses.append(fix_dof_process)
            free_dof_process = KratosSolid.FreeScalarDofProcess(self.model_part, params)
            self.FreeDofsProcesses.append(free_dof_process)

            params["variable_name"].SetString(self.derivated_variable_name+"_Z")
            fix_dof_process  =  KratosSolid.FixScalarDofProcess(self.model_part, params)
            self.FixDofsProcesses.append(fix_dof_process)
            free_dof_process = KratosSolid.FreeScalarDofProcess(self.model_part, params)
            self.FreeDofsProcesses.append(free_dof_process)

            params["variable_name"].SetString(self.settings["variable_name"].GetString())
            if( self.fix_time_integration == False ):

                params["variable_name"].SetString(self.variable_name+"_X")
                fix_dof_process  =  KratosSolid.FixScalarDofProcess(self.model_part, params)
                self.FixDofsProcesses.append(fix_dof_process)
                free_dof_process = KratosSolid.FreeScalarDofProcess(self.model_part, params)
                self.FreeDofsProcesses.append(free_dof_process)

                params["variable_name"].SetString(self.variable_name+"_Y")
                fix_dof_process  =  KratosSolid.FixScalarDofProcess(self.model_part, params)
                self.FixDofsProcesses.append(fix_dof_process)
                free_dof_process = KratosSolid.FreeScalarDofProcess(self.model_part, params)
                self.FreeDofsProcesses.append(free_dof_process)

                params["variable_name"].SetString(self.variable_name+"_Z")
                fix_dof_process  =  KratosSolid.FixScalarDofProcess(self.model_part, params)
                self.FixDofsProcesses.append(fix_dof_process)
                free_dof_process = KratosSolid.FreeScalarDofProcess(self.model_part, params)
                self.FreeDofsProcesses.append(free_dof_process)

                params["variable_name"].SetString(self.variable_name)

        else:
            if( "DISPLACEMENT" == self.variable_name ):
                self.SetLinearTimeIntegration()

            params["variable_name"].SetString(self.variable_name+"_X")
            fix_dof_process  =  KratosSolid.FixScalarDofProcess(self.model_part, params)
            self.FixDofsProcesses.append(fix_dof_process)
            free_dof_process = KratosSolid.FreeScalarDofProcess(self.model_part, params)
            self.FreeDofsProcesses.append(free_dof_process)

            params["variable_name"].SetString(self.variable_name+"_Y")
            fix_dof_process  =  KratosSolid.FixScalarDofProcess(self.model_part, params)
            self.FixDofsProcesses.append(fix_dof_process)
            free_dof_process = KratosSolid.FreeScalarDofProcess(self.model_part, params)
            self.FreeDofsProcesses.append(free_dof_process)

            params["variable_name"].SetString(self.variable_name+"_Z")
            fix_dof_process  =  KratosSolid.FixScalarDofProcess(self.model_part, params)
            self.FixDofsProcesses.append(fix_dof_process)
            free_dof_process = KratosSolid.FreeScalarDofProcess(self.model_part, params)
            self.FreeDofsProcesses.append(free_dof_process)

            params["variable_name"].SetString(self.variable_name)

    #
    def CreateAssignmentProcess(self, params):

        # the rotation or angular dynamic variables are only to set the rigid body process correctly
        # inside they assign displacement or linear dynamic variables
        if( self.fix_derivated_variable == False ):
            self.variable_name = "ROTATION"
        else:
            for dynamic_variable in self.LinearDynamicVariables:
                counter = 0
                if dynamic_variable == self.variable_name:
                    self.variable_name = self.LinearDynamicVariables[counter]
                    break
                counter = counter + 1

        params["variable_name"].SetString(self.variable_name)

        params.AddValue("direction", self.custom_settings["direction"])
        params.AddValue("center", self.custom_settings["center"])
        if( self.value_is_numeric ):
            params.AddEmptyValue("modulus").SetDouble(self.value)
            self.AssignValueProcess = KratosSolid.AssignRotationAboutAnAxisToNodesProcess(self.model_part, params)
        else:
            self.AssignValueProcess = KratosSolid.AssignRotationFieldAboutAnAxisToNodesProcess(self.model_part, self.compiled_function, "function",  self.value_is_spatial_function, params)

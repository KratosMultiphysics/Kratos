from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.SolidMechanicsApplication as KratosSolid

## This proces sets the value of a vector variable component-by-component.
## In this case, the fixicity is given set by deffault to true.
import sys

def Factory(custom_settings, Model):
    if( not isinstance(custom_settings,KratosMultiphysics.Parameters) ):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return AssignVectorComponentsToNodesProcess(Model, custom_settings["Parameters"])

## All the processes python should be derived from "Process"
class AssignVectorComponentsToNodesProcess(KratosMultiphysics.Process):
    def __init__(self, Model, custom_settings ):
        KratosMultiphysics.Process.__init__(self)

        ##settings string in json format
        default_settings = KratosMultiphysics.Parameters("""
        {
             "help": "This process assigns a vector value to a vector variable component by component",
             "model_part_name": "MODEL_PART_NAME",
             "variable_name": "VARIABLE_NAME",
             "value": [0.0, 0.0, 0.0],
             "compound_assignment": "direct",
             "constrained":true,
             "interval": [0.0, "End"],
             "local_axes" : {}
        }
        """)

        ##overwrite the default settings with user-provided parameters
        self.settings = custom_settings
        self.settings.ValidateAndAssignDefaults(default_settings)

        ##check if variable type is a vector
        self.var = KratosMultiphysics.KratosGlobals.GetVariable(self.settings["variable_name"].GetString())
        if( not isinstance(self.var,KratosMultiphysics.Array1DVariable3) ):
            raise Exception("Variable type is incorrect. Must be a three-component vector.")

        self.model         = Model
        self.variable_name = self.settings["variable_name"].GetString()

        ###check component assignation
        self.AssignValueProcesses = []

        self.constraints = []
        for i in range(0, self.settings["value"].size() ):
            if( self.settings["value"][i].IsNull() ):
                self.constraints.append(False)
            else:
                self.constraints.append(True)

        #print(" constraints ", self.constraints)

        self.BuildComponentsProcesses()

    def GetVariables(self):
        nodal_variables = [self.settings["variable_name"].GetString()]
        return nodal_variables

    def ExecuteInitialize(self):
        for process in self.AssignValueProcesses:
            process.ExecuteInitialize()


    def ExecuteInitializeSolutionStep(self):
        for process in self.AssignValueProcesses:
            process.ExecuteInitializeSolutionStep()


    def ExecuteFinalizeSolutionStep(self):
        for process in self.AssignValueProcesses:
            process.ExecuteFinalizeSolutionStep()

    #
    def BuildComponentsProcesses(self):

        counter = 0
        for imposed in self.constraints:

            if( imposed ):
                params = KratosMultiphysics.Parameters("{}")
                params.AddValue("model_part_name", self.settings["model_part_name"])

                if( counter == 0 ):
                    params.AddEmptyValue("variable_name").SetString(self.variable_name+"_X")
                if( counter == 1 ):
                    params.AddEmptyValue("variable_name").SetString(self.variable_name+"_Y")
                if( counter == 2 ):
                    params.AddEmptyValue("variable_name").SetString(self.variable_name+"_Z")

                params.AddValue("interval",self.settings["interval"])
                params.AddValue("constrained", self.settings["constrained"])
                params.AddValue("compound_assignment", self.settings["compound_assignment"])

                if( self.settings["value"][counter].IsNumber() ):
                    params.AddEmptyValue("value").SetDouble(self.settings["value"][counter].GetDouble())
                    #print(" Value ", counter," ", self.settings["value"][counter].GetDouble() )
                else:
                    params.AddEmptyValue("value").SetString(self.settings["value"][counter].GetString())

                params.AddValue("local_axes", self.settings["local_axes"])

                import assign_scalar_to_nodes_process as assign_scalar_process

                self.AssignValueProcesses.append(assign_scalar_process.AssignScalarToNodesProcess(self.model, params))

            counter +=1

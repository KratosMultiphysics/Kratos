from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.SolidMechanicsApplication as KratosSolid

## This proces sets the value of a scalar variable to conditions

import assign_scalar_to_nodes_process as BaseProcess

def Factory(custom_settings, Model):
    if( not isinstance(custom_settings,KratosMultiphysics.Parameters) ):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    return AssignScalarToConditionsProcess(Model, custom_settings["Parameters"])


class AssignScalarToConditionsProcess(BaseProcess.AssignScalarToNodesProcess):
    def __init__(self, Model, custom_settings ):
        BaseProcess.AssignScalarToNodesProcess.__init__(self, Model, custom_settings)

    def ExecuteInitialize(self):

        # set model part
        self.model_part = self.model[self.settings["model_part_name"].GetString()]

        if not self.model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED]:
            self.model_part.ProcessInfo.SetValue(KratosMultiphysics.INTERVAL_END_TIME, self.interval[1])

        # set processes
        params = KratosMultiphysics.Parameters("{}")
        params.AddValue("model_part_name", self.settings["model_part_name"])
        params.AddValue("compound_assignment", self.settings["compound_assignment"])

        self.CreateAssignmentProcess(params)

        self.SetCurrentTime()
        if( self.IsInsideInterval() and (self.interval_string == "initial" or self.interval_string == "start") ):
            self.AssignValueProcess.Execute()


    def ExecuteFinalizeSolutionStep(self):

        if not self.interval_ended:
            current_time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]
            delta_time   = self.model_part.ProcessInfo[KratosMultiphysics.DELTA_TIME]

            #arithmetic floating point tolerance
            tolerance = delta_time * 0.001

            if( (current_time + delta_time) > (self.interval[1] + tolerance) ):
                self.interval_ended = True
                if not self.finalized :
                    self.AssignValueProcess.ExecuteFinalize()
                    self.finalized = True

    def ExecuteAssignment(self):
        if self.IsInsideInterval():
            self.AssignValueProcess.Execute()

    def ExecuteUnAssignment(self):
        if self.IsInsideInterval():
            self.UnAssignValueProcess.Execute()

    def CreateAssignmentProcess(self, params):
        if( self.value_is_numeric ):
            params.AddValue("variable_name", self.settings["variable_name"])
            params.AddValue("value", self.settings["value"])
            params.AddEmptyValue("entity_type").SetString("CONDITIONS")
            self.AssignValueProcess = KratosSolid.AssignScalarToEntitiesProcess(self.model_part, params)
        else:
            #function values are assigned to a vector variable :: transformation is needed
            if( isinstance(self.var,KratosMultiphysics.DoubleVariable) ):
                variable_name = self.settings["variable_name"].GetString() + "_VECTOR"
                #print("::[--Assign_Variable--]:: "+variable_name)
                params.AddEmptyValue("variable_name")
                params["variable_name"].SetString(variable_name)
            else:
                params.AddValue("variable_name", self.settings["variable_name"])

            params.AddEmptyValue("entity_type").SetString("CONDITIONS")
            self.AssignValueProcess = KratosSolid.AssignScalarFieldToEntitiesProcess(self.model_part, self.compiled_function, "function", self.value_is_spatial_function, params)

        # in case of going to previous time step for time step reduction
        self.CreateUnAssignmentProcess(params)


    def CreateUnAssignmentProcess(self, params):
        params["compound_assignment"].SetString(self.GetInverseAssigment(self.settings["compound_assignment"].GetString()))
        if( self.value_is_numeric ):
            self.UnAssignValueProcess = KratosSolid.AssignScalarToEntitiesProcess(self.model_part, params)
        else:
            self.UnAssignValueProcess = KratosSolid.AssignScalarFieldToEntitiesProcess(self.model_part, self.compiled_function, "function", self.value_is_spatial_function, params)

from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.PfemFluidDynamicsApplication as PfemFluid

## This proces sets the value of a vector variable using a modulus and direction
## In this case, the fixicity is given set by deffault to false.
import sys
from math import *

class compiled_time_spatial_function:
    def __init__(self, compiled_function ):
        self.compiled_function = compiled_function

    def function(self,x,y,z,t):
        return eval(self.compiled_function)

def Factory(custom_settings, Model):
    if( not isinstance(custom_settings,KratosMultiphysics.Parameters) ):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return AssignModulusAndDirectionToConditionsProcess(Model, custom_settings["Parameters"])

## All the processes python should be derived from "Process"
class AssignModulusAndDirectionToConditionsProcess(KratosMultiphysics.Process):
    def __init__(self, Model, custom_settings ):
        KratosMultiphysics.Process.__init__(self)

        ##settings string in json format
        default_settings = KratosMultiphysics.Parameters("""
        {
             "help" : "This process assigns a modulus and direction value to a conditions vector variable",
             "model_part_name": "MODEL_PART_NAME",
             "variable_name": "VARIABLE_NAME",
             "modulus" : 0.0,
             "direction": [0.0, 0.0, 0.0],
             "compound_assignment": "direct",
             "constrained": false,
             "interval": [0.0, "End"],
             "local_axes" : {}
        }
        """)

        #trick to allow "value" to be a string or a double value
        if(custom_settings.Has("modulus")):
            if(custom_settings["modulus"].IsString()):
                default_settings["modulus"].SetString("0.0")


        ##overwrite the default settings with user-provided parameters
        self.settings = custom_settings
        self.settings.ValidateAndAssignDefaults(default_settings)

        ##check if variable type is a vector
        self.var = KratosMultiphysics.KratosGlobals.GetVariable(self.settings["variable_name"].GetString())
        if( (not isinstance(self.var,KratosMultiphysics.VectorVariable)) and (not isinstance(self.var,KratosMultiphysics.Array1DVariable3)) ):
            raise Exception("Variable type is incorrect. Must be a vector or an array_1d vector.")

        self.model         = Model
        self.variable_name = self.settings["variable_name"].GetString()

        ## set the interval
        self.interval_ended = False
        self.finalized = False
        self.interval  = []
        self.interval.append(self.settings["interval"][0].GetDouble())
        if( self.settings["interval"][1].IsString() ):
            if( self.settings["interval"][1].GetString() == "End" ):
                self.interval.append(sys.float_info.max)
        elif( self.settings["interval"][1].IsDouble() or  self.settings["interval"][1].IsInt() ):
            self.interval.append(self.settings["interval"][1].GetDouble())

        self.interval_string = "custom"
        if( self.interval[0] == 0.0 and self.interval[1] == 0.0 ):
            self.interval_string = "initial"
        elif( self.interval[0] < 0 ):
            self.interval_string = "start"
            self.interval[0] = 0.0
            
        ##check normalized direction
        direction   = []
        scalar_prod = 0
        for i in range(0, self.settings["direction"].size() ):
            direction.append( self.settings["direction"][i].GetDouble() )
            scalar_prod = scalar_prod + direction[i]*direction[i]

        norm = sqrt(scalar_prod)

        self.value = []
        if( norm != 0.0 ):
            for j in direction:
                self.value.append( j/norm )
        else:
            for j in direction:
                self.value.append(0.0)

        ## set the value
        self.value_is_numeric = False

        if self.settings["modulus"].IsNumber():

            self.value_is_numeric = True

            modulus = self.settings["modulus"].GetDouble()

            for i in range(0, len(self.value)):
                self.value[i] *= modulus

        else:

            self.function_expression = self.settings["modulus"].GetString()

            if (sys.version_info > (3, 0)):
                self.compiled_function = compiled_time_spatial_function(compile(self.function_expression, '', 'eval', optimize=2))
            else:
                self.compiled_function = compiled_time_spatial_function(compile(self.function_expression, '', 'eval'))

            self.value_is_spatial_function = True

            if(self.function_expression.find("x") == -1 and
               self.function_expression.find("y") == -1 and
               self.function_expression.find("z") == -1): #depends on time only
                    self.value_is_spatial_function = False


    def GetVariables(self):
        nodal_variables = [self.settings["variable_name"].GetString()]
        return nodal_variables

    def ExecuteInitialize(self):
        # set model part
        self.model_part = self.model[self.settings["model_part_name"].GetString()]
        if not self.model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED]:
            self.model_part.ProcessInfo.SetValue(KratosMultiphysics.INTERVAL_END_TIME, self.interval[1])

        # set processes
        params = KratosMultiphysics.Parameters("{}")
        params.AddValue("model_part_name", self.settings["model_part_name"])
        params.AddValue("compound_assignment", self.settings["compound_assignment"])

        params.AddEmptyValue("value")
        params.__setitem__("value", self.settings["direction"])

        self.CreateAssignmentProcess(params)

        self.SetCurrentTime()
        if( self.IsInsideInterval() and (self.interval_string == "initial" or self.interval_string == "start") ):
            self.AssignValueProcess.Execute()


    def ExecuteInitializeSolutionStep(self):
        # if self.IsRecoverStep():
        #     self.SetPreviousTime()
        #     self.ExecuteUnAssignment()

        self.SetCurrentTime()
        self.ExecuteAssignment()

    def ExecuteFinalizeSolutionStep(self):
        if not self.interval_ended:
            self.SetCurrentTime()
            #arithmetic floating point tolerance
            tolerance = self.delta_time * 0.001

            if( (self.current_time + self.delta_time) > (self.interval[1] + tolerance) ):
                self.interval_ended = True
                if not self.finalized :
                    self.AssignValueProcess.ExecuteFinalize()
                    self.finalized = True
    #
    def CreateAssignmentProcess(self, params):

        if( self.value_is_numeric ):
            params.AddValue("variable_name", self.settings["variable_name"])
            counter = 0
            for i in self.value:
                params["value"][counter].SetDouble(i)
                counter+=1
            self.AssignValueProcess = PfemFluid.AssignVectorToConditionsProcess(self.model_part, params)
        else:
            #function values are assigned to a vector variable :: transformation is needed
            if( isinstance(self.var,KratosMultiphysics.Array1DVariable3) ):
                variable_name = self.settings["variable_name"].GetString() + "_VECTOR"
                #print("::[--Assign_Variable--]:: "+variable_name)
                params.AddEmptyValue("variable_name")
                params["variable_name"].SetString(variable_name)
            else:
                params.AddValue("variable_name", self.settings["variable_name"])
            counter = 0
            for i in self.value:
                params["value"][counter].SetDouble(i)
                counter+=1
            params.AddEmptyValue("entity_type").SetString("CONDITIONS")
            self.AssignValueProcess = PfemFluid.AssignVectorFieldToEntitiesProcess(self.model_part, self.compiled_function, "function",  self.value_is_spatial_function, params)

        # in case of going to previous time step for time step reduction
        self.CreateUnAssignmentProcess(params)

    #
    def CreateUnAssignmentProcess(self, params):
        params["compound_assignment"].SetString(self.GetInverseAssigment(self.settings["compound_assignment"].GetString()))
        if( self.value_is_numeric ):
            self.UnAssignValueProcess = PfemFluid.AssignVectorToConditionsProcess(self.model_part, params)
        else:
            self.UnAssignValueProcess = PfemFluid.AssignVectorFieldToEntitiesProcess(self.model_part, self.compiled_function, "function",  self.value_is_spatial_function, params)

    #
    def ExecuteAssignment(self):
        if self.IsInsideInterval():
            self.AssignValueProcess.Execute()
    #
    def ExecuteUnAssignment(self):
        if self.IsInsideInterval():
            self.UnAssignValueProcess.Execute()
    #
    @classmethod
    def GetInverseAssigment(self,compound_assignment):
        if compound_assignment == "direct":
            return "direct"
        if compound_assignment == "addition":
            return "subtraction"
        if compound_assignment == "subtraction":
            return "addition"
        if compound_assignment == "multiplication":
            return "division"
        if compound_assignment == "division":
            return "multiplication"

    #
    def SetCurrentTime(self):
        self.delta_time = self.model_part.ProcessInfo[KratosMultiphysics.DELTA_TIME]
        self.current_time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]
        self.previous_time = self.current_time-self.delta_time
    #
    def SetPreviousTime(self):
        self.delta_time = self.model_part.ProcessInfo.GetPreviousSolutionStepInfo()[KratosMultiphysics.DELTA_TIME]
        self.current_time = self.model_part.ProcessInfo.GetPreviousSolutionStepInfo()[KratosMultiphysics.TIME]
        self.previous_time = self.current_time-self.delta_time
    #
    def IsInsideInterval(self):
        #arithmetic floating point tolerance
        tolerance = self.delta_time * 0.001

        if( self.current_time >= (self.interval[0] - tolerance) and self.current_time <= (self.interval[1] + tolerance) ):
            self.interval_ended = False
            return True
        else:
            return False

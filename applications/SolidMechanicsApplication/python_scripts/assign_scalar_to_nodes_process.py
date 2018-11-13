from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.SolidMechanicsApplication as KratosSolid


## This proces sets the value of a scalar variable
## Note that in some cases can be a vector of scalars (used in conditions with multiple nodes)
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
    return AssignScalarToNodesProcess(Model, custom_settings["Parameters"])

## All the processes python should be derived from "Process"
class AssignScalarToNodesProcess(KratosMultiphysics.Process):
    def __init__(self, Model, custom_settings ):
        KratosMultiphysics.Process.__init__(self)

        ##settings string in json format
        default_settings = KratosMultiphysics.Parameters("""
        {
             "help" : "This process assigns a scalar value to a scalar variable",
             "model_part_name": "MODEL_PART_NAME",
             "variable_name": "VARIABLE_NAME",
             "value": 0.0,
             "compound_assignment": "direct",
             "constrained": true,
             "interval": [0.0, "End"],
             "local_axes" : {}
        }
        """)


        #trick to allow "value" to be a string or a double value
        if(custom_settings.Has("value")):
            if(custom_settings["value"].IsString()):
                default_settings["value"].SetString("0.0")

        ##overwrite the default settings with user-provided parameters
        self.settings = custom_settings
        self.settings.ValidateAndAssignDefaults(default_settings)

        ##check if variable type is a scalar or a vector component
        self.CheckVariableType(self.settings["variable_name"].GetString())

        self.model         = Model
        self.variable_name = self.settings["variable_name"].GetString()

        ## set the interval
        self.finalized = False
        self.interval_started = False
        self.interval_ended   = False

        self.interval  = []
        self.interval.append(self.settings["interval"][0].GetDouble());
        if( self.settings["interval"][1].IsString() ):
            if( self.settings["interval"][1].GetString() == "End" ):
                self.interval.append(sys.float_info.max)
        elif( self.settings["interval"][1].IsDouble() or  self.settings["interval"][1].IsInt() ):
            self.interval.append(self.settings["interval"][1].GetDouble());


        self.interval_string = "custom"
        if( self.interval[0] == 0.0 and self.interval[1] == 0.0 ):
            self.interval_string = "initial"
        elif( self.interval[0] < 0 ):
            self.interval_string = "start"

        ## set the value
        self.value_is_numeric = False
        self.value_is_spatial_function = False
        self.value_is_current_value = False

        if self.settings["value"].IsNumber():
            self.value_is_numeric = True
            self.value = self.settings["value"].GetDouble()
        else:
            self.function_expression = self.settings["value"].GetString()

            if( self.function_expression == "current" ):
                self.value_is_current_value = True
            else:
                if (sys.version_info > (3, 0)):
                    self.compiled_function = compiled_time_spatial_function(compile(self.function_expression, '', 'eval', optimize=2))
                else:
                    self.compiled_function = compiled_time_spatial_function(compile(self.function_expression, '', 'eval'))

                self.value_is_spatial_function = True

                if(self.function_expression.find("x") == -1 and
                   self.function_expression.find("y") == -1 and
                   self.function_expression.find("z") == -1): #depends on time only
                    self.value_is_spatial_function = False


        self.constrained = self.settings["constrained"].GetBool()


    def ExecuteInitialAssignment(self):
        self.AssignValueProcess.Execute()
        # initial assignment no time integration
        #if( self.fix_time_integration ):
        #    for node in self.model_part.Nodes:
        #        self.TimeIntegrationMethod.Assign(node)

    def GetVariables(self):
        nodal_variables = [self.settings["variable_name"].GetString()]
        return nodal_variables

    def ExecuteInitialize(self):

        # set model part
        self.model_part = self.model[self.settings["model_part_name"].GetString()]

        if( self.model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED] == False ):
            self.model_part.ProcessInfo.SetValue(KratosMultiphysics.INTERVAL_END_TIME, self.interval[1])

        # set time integration method
        self.SetTimeIntegration()

        # set processes
        self.FixDofsProcesses     = []
        self.FreeDofsProcesses    = []

        params = KratosMultiphysics.Parameters("{}")
        params.AddValue("model_part_name", self.settings["model_part_name"])
        params.AddValue("variable_name", self.settings["variable_name"])

        self.fix_derivated_variable = False
        if( self.interval_string != "initial" and self.constrained == True ):
            self.SetFixAndFreeProcesses(params)

        params.AddValue("compound_assignment", self.settings["compound_assignment"])
        self.CreateAssignmentProcess(params)

        self.SetCurrentTime()
        if self.IsInsideInterval():
            if self.IsFixingStep():
                for process in self.FixDofsProcesses:
                    process.Execute()

            if( self.interval_string == "initial" or self.interval_string == "start" ):
                self.ExecuteInitialAssignment()

    def ExecuteInitializeSolutionStep(self):

        if self.IsRecoverStep():
            self.SetPreviousTime()
            self.ExecuteUnAssignment()

        self.SetCurrentTime()
        self.ExecuteAssignment()


    def ExecuteFinalizeSolutionStep(self):

        self.SetCurrentTime()
        if self.IsUnfixingStep():

            for process in self.FreeDofsProcesses:
                process.Execute()


    #
    def ExecuteAssignment(self):
        if self.IsInsideInterval():
            if self.IsFixingStep():
                for process in self.FixDofsProcesses:
                    process.Execute()
            self.AssignValueProcess.Execute()
            if( self.fix_time_integration ):
                for node in self.model_part.Nodes:
                    self.TimeIntegrationMethod.Assign(node)
    #
    def ExecuteUnAssignment(self):
        if self.IsInsideInterval():
            self.UnAssignValueProcess.Execute()
            if( self.fix_time_integration ):
                for node in self.model_part.Nodes:
                    self.TimeIntegrationMethod.Assign(node)

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
    def CheckVariableType(self,name):

        self.var = KratosMultiphysics.KratosGlobals.GetVariable(name)
        if( (not isinstance(self.var,KratosMultiphysics.Array1DComponentVariable)) and (not isinstance(self.var,KratosMultiphysics.DoubleVariable)) and (not isinstance(self.var,KratosMultiphysics.VectorVariable)) ):
            raise Exception("Variable type is incorrect. Must be a scalar or a component")


    #
    def SetFixAndFreeProcesses(self,params):

        if( self.fix_time_integration == True ):
            params["variable_name"].SetString(self.settings["variable_name"].GetString())
            self.primary_variable_name = self.TimeIntegrationMethod.GetPrimaryVariableName()
            if( self.primary_variable_name != self.variable_name ):
                params["variable_name"].SetString(self.primary_variable_name)

        fix_dof_process  =  KratosSolid.FixScalarDofProcess(self.model_part, params)
        self.FixDofsProcesses.append(fix_dof_process)
        free_dof_process = KratosSolid.FreeScalarDofProcess(self.model_part, params)
        self.FreeDofsProcesses.append(free_dof_process)

    #
    def SetTimeIntegration(self):
        self.fix_time_integration  = False

        self.TimeIntegrationMethod = None
        time_integration_container = KratosSolid.ComponentTimeIntegrationMethods()
        if( time_integration_container.HasProcessInfo(KratosSolid.COMPONENT_TIME_INTEGRATION_METHODS, self.model_part.ProcessInfo) ):
            time_integration_methods = time_integration_container.GetFromProcessInfo(KratosSolid.COMPONENT_TIME_INTEGRATION_METHODS, self.model_part.ProcessInfo)

            if( time_integration_methods.Has(self.variable_name) ):
                self.TimeIntegrationMethod = time_integration_methods.Get(self.variable_name).Clone()
            else:
                method_variable_name = time_integration_methods.GetMethodVariableName(self.variable_name)
                if( method_variable_name != self.variable_name ):
                    self.TimeIntegrationMethod = time_integration_methods.Get(method_variable_name).Clone()


        if( self.TimeIntegrationMethod != None ):
            self.fix_time_integration = True
            #set input variable
            input_variable = KratosMultiphysics.KratosGlobals.GetVariable(self.variable_name)
            self.TimeIntegrationMethod.SetInputVariable(input_variable)
        #else:
        #    print(self.variable_name+": No time integration ")

    #
    def CreateAssignmentProcess(self, params):

        params["variable_name"].SetString(self.settings["variable_name"].GetString())
        if( self.value_is_numeric ):
            params.AddEmptyValue("value").SetDouble(self.value)
            params.AddEmptyValue("entity_type").SetString("NODES")
            self.AssignValueProcess = KratosSolid.AssignScalarToEntitiesProcess(self.model_part, params)
        else:
            if( self.value_is_current_value ):
                self.AssignValueProcess = KratosMultiphysics.Process() #void process
            else:
                params.AddEmptyValue("entity_type").SetString("NODES")
                self.AssignValueProcess = KratosSolid.AssignScalarFieldToEntitiesProcess(self.model_part, self.compiled_function, "function",  self.value_is_spatial_function, params)


        # in case of going to previous time step for time step reduction
        self.CreateUnAssignmentProcess(params)

    #
    def CreateUnAssignmentProcess(self, params):
        params["compound_assignment"].SetString(self.GetInverseAssigment(self.settings["compound_assignment"].GetString()))
        if( self.value_is_numeric ):
            self.UnAssignValueProcess = KratosSolid.AssignScalarToEntitiesProcess(self.model_part, params)
        else:
            if( self.value_is_current_value ):
                self.UnAssignValueProcess = KratosMultiphysics.Process() #void process
            else:
                self.UnAssignValueProcess = KratosSolid.AssignScalarFieldToEntitiesProcess(self.model_part, self.compiled_function, "function",  self.value_is_spatial_function, params)

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
    def IsRecoverStep(self):
        if self.model_part.ProcessInfo.Has(KratosSolid.DELTA_TIME_CHANGED):
            if self.model_part.ProcessInfo[KratosSolid.DELTA_TIME_CHANGED] is True:
                return True
            else:
                return False
        else:
            return False

    #
    def IsInsideInterval(self):

        #arithmetic floating point tolerance
        tolerance = self.delta_time * 0.001

        if( self.current_time >= (self.interval[0] - tolerance) and self.current_time <= (self.interval[1] + tolerance) ):
            self.interval_ended = False;
            return True
        else:
            return False

    #
    def IsFixingStep(self):

        if( self.interval_started == False ):
            self.interval_started = True
            return True
        else:
            interval_time = self.model_part.ProcessInfo[KratosMultiphysics.INTERVAL_END_TIME]

            if(self.previous_time == interval_time):
                return True
            else:
                return False

    #
    def IsUnfixingStep(self):

        if( self.interval_ended == False ):

            #arithmetic floating point tolerance
            tolerance = self.delta_time * 0.001

            if( (self.current_time + self.delta_time) > (self.interval[1] + tolerance) ):
                self.interval_ended = True
                self.model_part.ProcessInfo.SetValue(KratosMultiphysics.INTERVAL_END_TIME, self.current_time)
                return True
            else:
                return False
        else:
            return False

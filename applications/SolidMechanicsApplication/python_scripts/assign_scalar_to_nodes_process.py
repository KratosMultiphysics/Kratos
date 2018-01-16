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
    if(type(custom_settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return AssignScalarToNodesProcess(Model, custom_settings["Parameters"])

## All the processes python processes should be derived from "python_process"
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

        ## dynamic variables
        self.LinearDynamicVariables  = ["ACCELERATION","VELOCITY"]
        self.AngularDynamicVariables = ["ANGULAR_ACCELERATION","ANGULAR_VELOCITY"]

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

    def GetVariables(self):
        nodal_variables = [self.settings["variable_name"].GetString()]
        return nodal_variables

    def ExecuteInitialize(self):

        # set model part
        self.model_part = self.model[self.settings["model_part_name"].GetString()]

        if( self.model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED] == False ):
            self.model_part.ProcessInfo.SetValue(KratosMultiphysics.INTERVAL_END_TIME, self.interval[1])


        self.fix_time_integration  = False
        self.time_integration_type = None
        angular_variables = self.AngularDynamicVariables[:]
        angular_variables.append("ROTATION")
        linear_variables = self.LinearDynamicVariables[:]
        linear_variables.append("DISPLACEMENT")
        variable = self.variable_name[:-2]
        try:
            index_value = angular_variables.index(variable)
            self.time_integration_type ="Angular"
        except ValueError:
            self.time_integration_type = None

        if( self.time_integration_type == None ):
            try:
                index_value = linear_variables.index(variable)
                self.time_integration_type ="Linear"
            except ValueError:
                self.time_integration_type = None

        self.TimeIntegrationMethod = None
        time_integration_method = KratosSolid.TimeIntegrationMethod()
        if( self.time_integration_type == "Angular" ):
            if( time_integration_method.HasProcessInfo(KratosSolid.ANGULAR_TIME_INTEGRATION_METHOD, self.model_part.ProcessInfo) ):
                self.TimeIntegrationMethod = time_integration_method.GetFromProcessInfo(KratosSolid.ANGULAR_TIME_INTEGRATION_METHOD, self.model_part.ProcessInfo)
            else:
                print(variable+": No time integration ")

        elif( self.time_integration_type == "Linear" ):
            if( time_integration_method.HasProcessInfo(KratosSolid.TIME_INTEGRATION_METHOD, self.model_part.ProcessInfo) ):
                self.TimeIntegrationMethod = time_integration_method.GetFromProcessInfo(KratosSolid.TIME_INTEGRATION_METHOD, self.model_part.ProcessInfo)
            else:
                print(variable+": No time integration ")
        else:
            print(variable+": No time integration ")

        # set processes
        self.FixDofsProcesses     = []
        self.FreeDofsProcesses    = []

        params = KratosMultiphysics.Parameters("{}")
        params.AddValue("model_part_name", self.settings["model_part_name"])
        params.AddValue("variable_name", self.settings["variable_name"])

        self.fix_derivated_variable = False
        if( self.interval_string != "initial" and self.constrained == True ):
            self.SetFixAndFreeProcesses(params)

        self.CreateAssignmentProcess(params)

        if self.IsInsideInterval():
            if self.IsFixingStep():
                for process in self.FixDofsProcesses:
                    process.Execute()

            if( self.interval_string == "initial" ):
                self.AssignValueProcess.Execute()

                if( self.fix_time_integration ):
                    self.SetTimeIntegration()
                    for node in self.model_part.Nodes:
                        self.TimeIntegrationMethod.Predict(node)


    def ExecuteInitializeSolutionStep(self):

        if self.IsInsideInterval():

            if self.IsFixingStep():
                for process in self.FixDofsProcesses:
                    process.Execute()

            self.AssignValueProcess.Execute()

            if( self.fix_time_integration ):
                self.SetTimeIntegration()
                for node in self.model_part.Nodes:
                    self.TimeIntegrationMethod.Predict(node)

    def ExecuteFinalizeSolutionStep(self):

        if self.IsUnfixingStep():

            for process in self.FreeDofsProcesses:
                process.Execute()


    #
    def CheckVariableType(self,name):

        self.var = KratosMultiphysics.KratosGlobals.GetVariable(name)
        if(type(self.var) != KratosMultiphysics.Array1DComponentVariable and type(self.var) != KratosMultiphysics.DoubleVariable and type(self.var) != KratosMultiphysics.VectorVariable):
            raise Exception("Variable type is incorrect. Must be a scalar or a component")


    #
    def SetFixAndFreeProcesses(self,params):

        for dynamic_variable in self.AngularDynamicVariables:
            if dynamic_variable == self.variable_name[:-2]:
                self.derivated_variable_name = "ROTATION" + self.variable_name[-2:]
                self.fix_derivated_variable = True
                self.SetAngularTimeIntegration()
                break

        if( self.fix_derivated_variable == False ):
            for dynamic_variable in self.LinearDynamicVariables:
                if dynamic_variable == self.variable_name[:-2]:
                    self.derivated_variable_name = "DISPLACEMENT" + self.variable_name[-2:]
                    self.fix_derivated_variable = True
                    self.SetLinearTimeIntegration()
                    break

        if( self.fix_derivated_variable ):
            params["variable_name"].SetString(self.derivated_variable_name)
            fix_dof_process  =  KratosSolid.FixScalarDofProcess(self.model_part, params)
            self.FixDofsProcesses.append(fix_dof_process)
            free_dof_process = KratosSolid.FreeScalarDofProcess(self.model_part, params)
            self.FreeDofsProcesses.append(free_dof_process)
            params["variable_name"].SetString(self.settings["variable_name"].GetString())

            if( self.fix_time_integration == False ):
                params["variable_name"].SetString(self.settings["variable_name"].GetString())
                fix_dof_process  =  KratosSolid.FixScalarDofProcess(self.model_part, params)
                self.FixDofsProcesses.append(fix_dof_process)
                free_dof_process = KratosSolid.FreeScalarDofProcess(self.model_part, params)
                self.FreeDofsProcesses.append(free_dof_process)
        else:
            if( "ROTATION" == self.variable_name[:-2] ):
                self.SetAngularTimeIntegration()
            elif( "DISPLACEMENT" == self.variable_name[:-2] ):
                self.SetLinearTimeIntegration()

            params["variable_name"].SetString(self.settings["variable_name"].GetString())
            fix_dof_process  =  KratosSolid.FixScalarDofProcess(self.model_part, params)
            self.FixDofsProcesses.append(fix_dof_process)
            free_dof_process = KratosSolid.FreeScalarDofProcess(self.model_part, params)
            self.FreeDofsProcesses.append(free_dof_process)

    #
    def SetTimeIntegration(self):
        if( self.time_integration_type == "Linear" ):
            self.SetLinearTimeIntegration()
        elif( self.time_integration_type == "Angular" ):
            self.SetAngularTimeIntegration()

    #
    def SetLinearTimeIntegration(self):
        if( self.TimeIntegrationMethod != None ):
            variable          = KratosMultiphysics.KratosGlobals.GetVariable("DISPLACEMENT" + self.variable_name[-2:])
            first_derivative  = KratosMultiphysics.KratosGlobals.GetVariable("VELOCITY" + self.variable_name[-2:])
            second_derivative = KratosMultiphysics.KratosGlobals.GetVariable("ACCELERATION" + self.variable_name[-2:])
            self.TimeIntegrationMethod.SetVariables(variable, first_derivative, second_derivative)
            if( self.TimeIntegrationMethod.HasStepVariable() ):
                step_variable = KratosMultiphysics.KratosGlobals.GetVariable("STEP_DISPLACEMENT" + self.variable_name[-2:])
                self.TimeIntegrationMethod.SetStepVariable(step_variable)

            self.fix_time_integration  = True
            input_variable = KratosMultiphysics.KratosGlobals.GetVariable(self.variable_name)
            self.TimeIntegrationMethod.SetInputVariable(input_variable)


    #
    def SetAngularTimeIntegration(self):
        if( self.TimeIntegrationMethod != None ):
            variable          = KratosMultiphysics.KratosGlobals.GetVariable("ROTATION" + self.variable_name[-2:])
            first_derivative  = KratosMultiphysics.KratosGlobals.GetVariable("ANGULAR_VELOCITY" + self.variable_name[-2:])
            second_derivative = KratosMultiphysics.KratosGlobals.GetVariable("ANGULAR_ACCELERATION" + self.variable_name[-2:])
            self.TimeIntegrationMethod.SetVariables(variable, first_derivative, second_derivative)
            if( self.TimeIntegrationMethod.HasStepVariable() ):
                step_variable = KratosMultiphysics.KratosGlobals.GetVariable("STEP_ROTATION" + self.variable_name[-2:])
                self.TimeIntegrationMethod.SetStepVariable(step_variable)

            self.fix_time_integration = True
            input_variable = KratosMultiphysics.KratosGlobals.GetVariable(self.variable_name)
            self.TimeIntegrationMethod.SetInputVariable(input_variable)


    #
    def CreateAssignmentProcess(self, params):

        if( self.value_is_numeric ):
            params.AddEmptyValue("value").SetDouble(self.value)
            self.AssignValueProcess = KratosSolid.AssignScalarToNodesProcess(self.model_part, params)
        else:
            if( self.value_is_current_value ):
                self.AssignValueProcess = KratosMultiphysics.Process() #void process
            else:
                self.AssignValueProcess = KratosSolid.AssignScalarFieldToNodesProcess(self.model_part, self.compiled_function, "function",  self.value_is_spatial_function, params)


    #
    def IsInsideInterval(self):

        current_time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]
        delta_time   = self.model_part.ProcessInfo[KratosMultiphysics.DELTA_TIME]

        #arithmetic floating point tolerance
        tolerance = delta_time * 0.001

        if( current_time >= (self.interval[0] - tolerance) and current_time <= (self.interval[1] + tolerance) ):
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
            previous_time = self.model_part.ProcessInfo.GetPreviousSolutionStepInfo()[KratosMultiphysics.TIME]

            if(previous_time == interval_time):
                return True
            else:
                return False

    #
    def IsUnfixingStep(self):

        if( self.interval_ended == False ):

            current_time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]
            delta_time   = self.model_part.ProcessInfo[KratosMultiphysics.DELTA_TIME]

            #arithmetic floating point tolerance
            tolerance = delta_time * 0.001

            if( (current_time + delta_time) > (self.interval[1] + tolerance) ):
                self.interval_ended = True
                self.model_part.ProcessInfo.SetValue(KratosMultiphysics.INTERVAL_END_TIME, current_time)
                return True
            else:
                return False
        else:
            return False

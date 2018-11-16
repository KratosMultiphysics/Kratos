from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.SolidMechanicsApplication as KratosSolid

## This proces sets the value of a vector variable
import sys
import assign_modulus_and_direction_to_conditions_process as BaseProcess

def Factory(custom_settings, Model):
    if( not isinstance(custom_settings,KratosMultiphysics.Parameters) ):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    return AssignVectorToConditionsProcess(Model, custom_settings["Parameters"])

class AssignVectorToConditionsProcess(BaseProcess.AssignModulusAndDirectionToConditionsProcess):
    def __init__(self, Model, custom_settings ):
        BaseProcess.AssignModulusAndDirectionToConditionsProcess.__init__(self, Model, custom_settings)

    def ExecuteInitialize(self):

        # set model part
        self.model_part = self.model[self.settings["model_part_name"].GetString()]
        if not self.model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED]:
            self.model_part.ProcessInfo.SetValue(KratosMultiphysics.INTERVAL_END_TIME, self.interval[1])

        if( self.IsInsideInterval() and (self.interval_string == "initial" or self.interval_string == "start")  ):
            self.AssignValueToConditions()

    def ExecuteInitializeSolutionStep(self):

        if self.IsInsideInterval():
            self.AssignValueToConditions()


    def ExecuteFinalizeSolutionStep(self):
        pass


    #
    def BuildComponentsProcesses(self):

        ## set the interval
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
            self.interval[0] = 0.0
            
        ## set the value
        self.value_is_numeric = False
        self.value_is_spatial_function = False

        #deprecated:
        self.function_string = self.settings["time_function"].GetString()
        if( self.function_string == "constant" ):
            self.value_is_numeric = True
            self.value = self.settings["value"].GetDouble()
        else:
            self.value_is_numeric = False
            self.value_is_spatial_function = False
            if( self.function_string == "incremental" ):
                self.function_expression = "t"
            else:
                self.function_expression = self.function_string;

            if (sys.version_info > (3, 0)):
                self.compiled_function = compile(self.function_expression, '', 'eval', optimize=2)
            else:
                self.compiled_function = compile(self.function_expression, '', 'eval')
        #deprecated:

    #
    def function(self, t):
        return eval(self.compiled_function)

    #
    def AssignValueToConditions(self):

        if( self.value_is_numeric ):
            for cond in self.model_part.Conditions:
                cond.SetValue(self.var,self.value)
        else:
            current_time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]
            if( self.value_is_spatial_function == False ):
                function_value = self.function(current_time)

                current_value = []
                for i in range(0, len(self.value) ):
                    current_value.append(self.value[i] * function_value)

                for cond in self.model_part.Conditions:
                    cond.SetValue(self.var,current_value)
            else:
                raise Exception("Spatial Functions not implemented for CONDITIONS")


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

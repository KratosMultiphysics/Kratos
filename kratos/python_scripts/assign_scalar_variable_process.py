import KratosMultiphysics
import sys
from math import *


def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return AssignScalarVariableProcess(Model, settings["Parameters"])

class AssignScalarVariableProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings ):
        KratosMultiphysics.Process.__init__(self)

        #The intervale size can be larger than 2
        #The value can be a double, a string or a vector of double or string
        #The number of values can be equal to 1, the number of time value or the number of time intervals
        default_settings = KratosMultiphysics.Parameters("""
            {
                "mesh_id"         : 0,
                "model_part_name" : "please_specify_model_part_name",
                "variable_name"   : "SPECIFY_VARIABLE_NAME",
                "interval"        : [0.0, 1e30],
                "constrained"     : true,
                "value"           : 0.0,
                "local_axes"      : {},
                "step_type"       : ""
            }
            """
            )

        #assign this here since it will change the "interval" prior to validation
        self.interval = KratosMultiphysics.IntervalUtility(settings)

        #here i do a trick, since i want to allow "value" to be a string or a double value
        if(settings.Has("value")):
            if(settings["value"].IsString()):
                default_settings["value"].SetString("0.0")

        #here i do a trick, since i want to allow "value" to be a vector or a double value
        if(settings.Has("value")):
            if(settings["value"].IsArray()):
                default_settings["value"].SetVector([0.0])

        settings.ValidateAndAssignDefaults(default_settings)

        self.variable = KratosMultiphysics.KratosGlobals.GetVariable(settings["variable_name"].GetString())
        if(type(self.variable) != KratosMultiphysics.Array1DComponentVariable and type(self.variable) != KratosMultiphysics.DoubleVariable and type(self.variable) != KratosMultiphysics.VectorVariable):
            msg = "Error in AssignScalarToNodesProcess. Variable type of variable : " + settings["variable_name"].GetString() + " is incorrect . Must be a scalar or a component"
            raise Exception(msg)

        self.model_part = Model[settings["model_part_name"].GetString()]
        self.mesh = self.model_part.GetMesh(settings["mesh_id"].GetInt())
        self.is_fixed = settings["constrained"].GetBool()

        self.values = []
        if settings["value"].IsArray():
            if (settings["value"].size()==self.interval.GetSubIntervalNb()):
                self.one_value_by_interval = True
            elif (settings["value"].size()==(self.interval.GetSubIntervalNb()+1)):
                self.one_value_by_interval = False
                self.step_type = settings["step_type"]
            else:
                msg = "Error the size of the value vector : " + settings["value"].size() + " is not consistant regarding the interval size: " + self.interval.GetSubIntervalNb()
                raise Exception(msg)
            value_size = settings["value"].size()
            for i in range(value_size):
                self.values.append(FunctionValue(settings["value"][i], settings["local_axes"], self.mesh, self.one_value_by_interval))
        else:
            self.one_value_by_interval = True
            self.interval.RemoveSubInterval()
            self.values.append(FunctionValue(settings["value"], settings["local_axes"], self.mesh, self.one_value_by_interval))

        #construct a variable_utils object to speedup fixing
        self.variable_utils = KratosMultiphysics.VariableUtils()
        self.step_is_active = False

    def ExecuteBeforeSolutionLoop(self):
        self.ExecuteInitializeSolutionStep()

    def ExecuteInitializeSolutionStep(self):
        current_time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]

        if(self.interval.IsInInterval(current_time)):

            self.step_is_active = True

            if(self.is_fixed):
                self.variable_utils.ApplyFixity(self.variable, self.is_fixed, self.mesh.Nodes)

            n = self.interval.GetSubInterval(current_time)
            if(self.one_value_by_interval):
                value = self.values[n-1]
                if not value.DependsOnSpace(): #depends on time only
                    self.variable_utils.SetScalarVar(self.variable, value.getValue(current_time), self.mesh.Nodes)
                else:
                    value.cpp_apply_function_utility.ApplyFunction(self.variable, current_time)
            else:
                if (not self.values[n-1].value_is_numeric) || (not self.values[n].value_is_numeric):
                    if (not self.values[n-1].value_is_numeric):
                        msg_val = self.values[n-1].function_string
                    elif (not self.values[n].value_is_numeric):
                        msg_val = self.values[n].function_string
                    raise Exception("Only numerical value are considered with assignment at time value:", msg_val)
                y1 = self.values[n-1].value
                y2 = self.values[n].value
                if ( abs( y2 - y1 ) < 1e-6 ):
                    value = y1
                else:
                    dy = y2 - y1;
                    t1 = self.interval.GetSubIntervalBegin(n)
                    t2 = self.interval.GetSubIntervalEnd(n)
                    dt = t2 - t1;
                    t = current_time
                    if( self.step_type.GetString() == "linear" ):
                        value = y1  + ( t - t1 ) * ( dy/dt )
                    elif( self.step_type.GetString() == "smooth" ):
                        from math import tanh as tanh  # we use hyperbolic tan to define a smooth step
                        t_av = 0.5 * (t1 + t2);
                        y_av = 0.5 * (y1 + y2);
                        c_slope = 0.15;
                        value = y_av + 0.5 * dy * tanh( (t-t_av)/(c_slope*dt) )
                    else:
                        raise Exception("Only linear and smooth are valid inputs")
                self.variable_utils.SetScalarVar(self.variable, value, self.mesh.Nodes)

    def ExecuteFinalizeSolutionStep(self):
        current_time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]

        if(self.step_is_active):
            #here we free all of the nodes in the mesh
            if(self.is_fixed):
                fixity_status  = False
                self.variable_utils.ApplyFixity(self.variable, fixity_status, self.mesh.Nodes)

        self.step_is_active = False

class FunctionValue():
    def __init__(self, value, local_axes, mesh, one_value_by_interval):
        if value.IsNumber():
            self.value_is_numeric = True
            self.value = value.GetDouble()
        else:
            self.value_is_numeric = False
            self.function_string = value.GetString()
            self.aux_function = KratosMultiphysics.PythonGenericFunctionUtility(self.function_string, local_axes)

            self.depends_on_space = self.aux_function.DependsOnSpace()
            depends_on_time = self.aux_function.DependsOnTime()
            if(depends_on_time & (not one_value_by_interval)):
                msg = "Error: time dependant value and assigned on specific time is incompatible: " + value
                raise Exception(msg)

            if(self.depends_on_space):
                self.cpp_apply_function_utility = KratosMultiphysics.ApplyFunctionToNodesUtility(mesh.Nodes, self.aux_function )

    def DependsOnSpace(self):
        if self.value_is_numeric:
            return False
        else:
            return self.depends_on_space

    def DependsOnTime(self):
        if self.value_is_numeric:
            return False
        else:
            return self.depends_on_time

    def getValue(self, time):
        if self.value_is_numeric:
            return self.value
        else:
            if self.depends_on_space:
                msg = "Error the function depends on space: " + self.function_string
                raise Exception(msg)
            return self.aux_function.CallFunction(0.0,0.0,0.0,time)

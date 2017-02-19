import KratosMultiphysics
import sys
from math import *


class aux_object_cpp_callback:
    def __init__(self, compiled_function ):
        self.compiled_function = compiled_function

    def f(self,x,y,z,t):
        return eval(self.compiled_function)

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return AssignScalarToNodesProcess(Model, settings["Parameters"])

class AssignScalarToNodesProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings ):
        KratosMultiphysics.Process.__init__(self)

        default_settings = KratosMultiphysics.Parameters("""
            {
                "mesh_id"         : 0,
                "model_part_name" : "please_specify_model_part_name",
                "variable_name"   : "SPECIFY_VARIABLE_NAME",
                "interval"        : [0.0, 1e30],
                "constrained"     : true,
                "value"           : 0.0,
                "local_axes"      : {}
            }
            """
            )

        #detect "End" as a tag and replace it by a large number
        if(settings.Has("interval")):
            if(settings["interval"][1].IsString() ):
                if(settings["interval"][1].GetString() == "End"):
                    settings["interval"][1].SetDouble(1e30) # = default_settings["interval"][1]
                else:
                    raise Exception("the second value of interval can be \"End\" or a number, interval currently:"+settings["interval"].PrettyPrintJsonString())

        #here i do a trick, since i want to allow "value" to be a string or a double value
        if(settings.Has("value")):            
            if(settings["value"].IsString()):
                default_settings["value"].SetString("0.0")



        settings.ValidateAndAssignDefaults(default_settings)

        #admissible values for local axes, are "empty" or both xorigin and axes (rotation matrix), prescribed
        #"local_axes"               :{
        #    "origin" : [0.0, 0.0, 0.0]
        #    "axes"  : [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0] ]
        #    }
        self.non_trivial_local_system = False
        if(settings["local_axes"].Has("origin")): # not empty
            self.non_trivial_local_system = True
            self.R = KratosMultiphysics.Matrix(3,3)
            self.xorigin = KratosMultiphysics.Vector(3)
            for i in range(3):
                self.xorigin[i] = settings["local_axes"]["origin"][i].GetDouble()                
                for j in range(3):
                    self.R[i,j] = settings["local_axes"]["axes"][i][j].GetDouble()

        self.variable = KratosMultiphysics.KratosGlobals.GetVariable(settings["variable_name"].GetString())
        if(type(self.variable) != KratosMultiphysics.Array1DComponentVariable and type(self.variable) != KratosMultiphysics.DoubleVariable and type(self.variable) != KratosMultiphysics.VectorVariable):
            msg = "Error in AssignScalarToNodesProcess. Variable type of variable : " + settings["variable_name"].GetString() + " is incorrect . Must be a scalar or a component"
            raise Exception(msg)

        self.model_part = Model[settings["model_part_name"].GetString()]
        self.mesh = self.model_part.GetMesh(settings["mesh_id"].GetInt())
        self.interval = KratosMultiphysics.Vector(2)
        self.interval[0] = settings["interval"][0].GetDouble()
        self.interval[1] = settings["interval"][1].GetDouble()
        self.is_fixed = settings["constrained"].GetBool()

        self.value_is_numeric = False
        self.is_time_function = False
        if settings["value"].IsNumber():
            self.value_is_numeric = True
            self.value = settings["value"].GetDouble()
        else:
            self.function_string = settings["value"].GetString()
            if (sys.version_info > (3, 0)):
                self.aux_function = aux_object_cpp_callback(compile(self.function_string, '', 'eval', optimize=2))
            else:
                self.aux_function = aux_object_cpp_callback(compile(self.function_string, '', 'eval'))

            if(self.function_string.find("x") == -1 and
               self.function_string.find("y") == -1 and
               self.function_string.find("z") == -1): #depends on time alone!
                    self.is_time_function = True
            else:
                if self.non_trivial_local_system == False: #function prescribed in the global system of coordinates
                    self.cpp_apply_function_utility = KratosMultiphysics.PythonGenericFunctionUtility(self.mesh.Nodes, self.aux_function )
                else:
                    self.cpp_apply_function_utility = KratosMultiphysics.PythonGenericFunctionUtility(self.mesh.Nodes, self.aux_function, self.R, self.xorigin )

        #construct a variable_utils object to speedup fixing
        self.variable_utils = KratosMultiphysics.VariableUtils()

    def ExecuteInitializeSolutionStep(self):
        current_time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]

        if(current_time >= self.interval[0] and  current_time<self.interval[1]):

            if(self.is_fixed):
                self.variable_utils.ApplyFixity(self.variable, self.is_fixed, self.mesh.Nodes)

            if self.value_is_numeric:
                self.variable_utils.SetScalarVar(self.variable, self.value, self.mesh.Nodes)
            else:
                if self.is_time_function:
                    self.value = self.aux_function.f(0.0,0.0,0.0,current_time)
                    self.variable_utils.SetScalarVar(self.variable, self.value, self.mesh.Nodes)
                else: #most general case - space varying function (possibly also time varying)
                    self.cpp_apply_function_utility.ApplyFunction(self.variable, current_time)
                    

    def ExecuteFinalizeSolutionStep(self):
        current_time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]

        if(current_time >= self.interval[0] and  current_time<self.interval[1]):
            #print("Freeing variable", self.variable)
            #here we free all of the nodes in the mesh
            if(self.is_fixed):
                fixity_status  = False
                self.variable_utils.ApplyFixity(self.variable, fixity_status, self.mesh.Nodes)

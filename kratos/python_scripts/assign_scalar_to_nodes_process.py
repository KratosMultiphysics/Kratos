import KratosMultiphysics
import sys
from math import *


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
        #admissible values for local axes, are "empty" or both xorigin and axes (rotation matrix), prescribed. An example follows
        #"local_axes"               :{
        #    "origin" : [0.0, 0.0, 0.0],
        #    "axes"  : [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0] ]
        #    }

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

        
        #self.non_trivial_local_system = False
        #if(settings["local_axes"].Has("origin")): # not empty
            #self.non_trivial_local_system = True
            #self.R = KratosMultiphysics.Matrix(3,3)
            #self.xorigin = KratosMultiphysics.Vector(3)
            #for i in range(3):
                #self.xorigin[i] = settings["local_axes"]["origin"][i].GetDouble()                
                #for j in range(3):
                    #self.R[i,j] = settings["local_axes"]["axes"][i][j].GetDouble()

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

        if settings["value"].IsNumber():
            self.value_is_numeric = True
            self.value = settings["value"].GetDouble()
        else:
            self.function_string = settings["value"].GetString()
            self.aux_function = KratosMultiphysics.PythonGenericFunctionUtility(self.function_string, settings["local_axes"])

            if(self.aux_function.DependsOnSpace()):
                self.cpp_apply_function_utility = KratosMultiphysics.ApplyFunctionToNodesUtility(self.mesh.Nodes, self.aux_function )
                
        #construct a variable_utils object to speedup fixing
        self.variable_utils = KratosMultiphysics.VariableUtils()
        self.step_is_active = False

    def ExecuteInitializeSolutionStep(self):
        current_time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]

        eps = max(1e-14*self.interval[0], 1e-30) #WARNING: very experimental
        if(current_time >= (self.interval[0]-eps) and  current_time< (self.interval[1]+eps) ): #WARNING: very experimental
            
            self.step_is_active = True
            
            if(self.is_fixed):
                self.variable_utils.ApplyFixity(self.variable, self.is_fixed, self.mesh.Nodes)

            if self.value_is_numeric:
                self.variable_utils.SetScalarVar(self.variable, self.value, self.mesh.Nodes)
            else:
                if self.aux_function.DependsOnSpace() == False: #depends on time only
                    self.value = self.aux_function.CallFunction(0.0,0.0,0.0,current_time)
                    self.variable_utils.SetScalarVar(self.variable, self.value, self.mesh.Nodes)
                else: #most general case - space varying function (possibly also time varying)
                    self.cpp_apply_function_utility.ApplyFunction(self.variable, current_time)
                    

    def ExecuteFinalizeSolutionStep(self):
        current_time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]

        if(self.step_is_active):
            #here we free all of the nodes in the mesh
            if(self.is_fixed):
                fixity_status  = False
                self.variable_utils.ApplyFixity(self.variable, fixity_status, self.mesh.Nodes)

        self.step_is_active = False
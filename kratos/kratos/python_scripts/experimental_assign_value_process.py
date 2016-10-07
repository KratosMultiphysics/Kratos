import KratosMultiphysics 
from math import *


class aux_object_cpp_callback:
    def __init__(self, compiled_function ):
        self.compiled_function = compiled_function

    def f(self,x,y,z,t):
        return eval(self.compiled_function)
        #return 0
        
def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return AssignValueProcess(Model, settings["Parameters"])

    
##all the processes python processes should be derived from "python_process"
class AssignValueProcess(KratosMultiphysics.Process):
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
        
        #here i do a trick, since i want to allow "value" to be a string or a double value
        if(settings.Has("value")):
            if(settings["value"].IsString()):
                default_settings["value"].SetString("0.0")
                
        
        
        settings.ValidateAndAssignDefaults(default_settings)

        #admissible values for local axes, are "empty" or 
        #"local_axes"               :{
        #    "origin" : [0.0, 0.0, 0.0]
        #    "axes"  : [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0] ] 
        #    }
        self.non_trivial_local_system = False
        if(settings["local_axes"].Has("origin")): # not empty
            self.non_trivial_local_system = True
            self.R = KratosMultiphysics.Matrix(3,3)
            self.R[0,0] = settings["local_axes"]["axes"][0][0].GetDouble()
            self.R[0,1] = settings["local_axes"]["axes"][0][1].GetDouble()
            self.R[0,2] = settings["local_axes"]["axes"][0][2].GetDouble()
            self.R[1,0] = settings["local_axes"]["axes"][1][0].GetDouble()
            self.R[1,1] = settings["local_axes"]["axes"][1][1].GetDouble()
            self.R[1,2] = settings["local_axes"]["axes"][1][2].GetDouble()
            self.R[2,0] = settings["local_axes"]["axes"][2][0].GetDouble()
            self.R[2,1] = settings["local_axes"]["axes"][2][1].GetDouble()
            self.R[2,2] = settings["local_axes"]["axes"][2][2].GetDouble()
            
            self.x0 = KratosMultiphysics.Vector(3)
            self.x0[0] = settings["local_axes"]["origin"][0].GetDouble()
            self.x0[1] = settings["local_axes"]["origin"][1].GetDouble()
            self.x0[2] = settings["local_axes"]["origin"][2].GetDouble()

        

        self.model_part = Model[settings["model_part_name"].GetString()]
        self.variable = getattr(KratosMultiphysics, settings["variable_name"].GetString())
        self.mesh = self.model_part.GetMesh(settings["mesh_id"].GetInt())
        self.interval = settings["interval"] #expected to be a vector of size 2
        self.is_fixed = settings["constrained"].GetBool()
        
        self.value_is_numeric = False
        self.is_time_function = False
        if settings["value"].IsNumber():
            self.value_is_numeric = True
            self.value = settings["value"].GetDouble()
        else:            
            self.function_string = settings["value"].GetString()
            self.aux_function = aux_object_cpp_callback(compile(self.function_string, '', 'eval', optimize=2))
            
            if(self.function_string.find("x") == -1 and 
               self.function_string.find("y") == -1 and
               self.function_string.find("z") == -1): #depends on time alone!
                    self.is_time_function = True
            else:
                self.cpp_apply_function_utility = KratosMultiphysics.PythonGenericFunctionUtility(self.mesh.Nodes, self.aux_function ) 
            
        
        #construct a variable_utils object to speedup fixing
        self.variable_utils = KratosMultiphysics.VariableUtils()
        
            
        print("finished construction of ApplyCustomFunctionProcess Process")
        
    def ExecuteInitializeSolutionStep(self):
        current_time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]
        
        if(current_time > self.interval[0].GetDouble() and  current_time<self.interval[1].GetDouble()):

            if(self.is_fixed):
                self.variable_utils.ApplyFixity(self.variable, self.is_fixed, self.mesh.Nodes)
                
            if self.value_is_numeric:
                self.variable_utils.SetScalarVar(self.variable, self.value, self.model_part.Nodes)
            else:
                if self.is_time_function:
                    self.value = self.aux_function.f(0.0,0.0,0.0,current_time)
                    self.variable_utils.SetScalarVar(self.variable, self.value, self.model_part.Nodes)
                else: #most general case - space varying function (possibly also time varying)
                    if self.non_trivial_local_system == False:
                        self.cpp_apply_function_utility.ApplyFunction(self.variable, current_time)
                    else: #TODO: OPTIMIZE!!
                        for node in self.model_part.Nodes:
                            dx = node - self.x0
                            x = self.R[0,0]*dx[0] + self.R[0,1]*dx[1] + self.R[0,2]*dx[2]
                            y = self.R[1,0]*dx[0] + self.R[1,1]*dx[1] + self.R[1,2]*dx[2]
                            z = self.R[2,0]*dx[0] + self.R[2,1]*dx[1] + self.R[2,2]*dx[2]
                            node.SetSolutionStepValue(self.variable,0,self.aux_function.f(x,y,z,current_time))
                        
            
    def ExecuteFinalizeSolutionStep(self):
        current_time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]
        
        if(current_time > self.interval[0].GetDouble() and  current_time<self.interval[1].GetDouble()):
            #here we free all of the nodes in the mesh
            if(self.is_fixed):
                fixity_status  = False
                self.variable_utils.ApplyFixity(self.variable, fixity_status, self.mesh.Nodes)

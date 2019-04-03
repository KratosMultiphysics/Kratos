import KratosMultiphysics
import sys
from math import *

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as KratosUnittest

class aux_object_cpp_callback:
    def __init__(self, compiled_function ):
        self.compiled_function = compiled_function

    def f(self,x,y,z,t):
        return eval(self.compiled_function)
        #return 0

def Factory(settings, Model):
    if( not isinstance(settings,KratosMultiphysics.Parameters) ):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return CheckScalarOnNodesProcess(Model, settings["Parameters"])


## All the processes python should be derived from "Process"
class CheckScalarOnNodesProcess(KratosMultiphysics.Process, KratosUnittest.TestCase):
    def __init__(self, Model, settings ):
        KratosMultiphysics.Process.__init__(self)

        default_settings = KratosMultiphysics.Parameters("""
            {
                "mesh_id"         : 0,
                "model_part_name" : "please_specify_model_part_name",
                "variable_name"   : "SPECIFY_VARIABLE_NAME",
                "interval"        : [0.0, 1e30],
                "value"           : 0.0,
                "tolerance_rank"  : 3,
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

        self.model    = Model
        self.settings = settings
        self.variable = getattr(KratosMultiphysics, settings["variable_name"].GetString())
        self.interval = KratosMultiphysics.Vector(2)
        self.interval[0] = settings["interval"][0].GetDouble()
        self.interval[1] = settings["interval"][1].GetDouble()

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

        # Error tolerance
        self.tol = settings["tolerance_rank"].GetInt()

        # print("Finished construction of ApplyCustomFunctionProcess Process")


    def ExecuteInitialize(self):
        self.model_part = self.model[self.settings["model_part_name"].GetString()]
        self.mesh = self.model_part.GetMesh(self.settings["mesh_id"].GetInt())

    def ExecuteInitializeSolutionStep(self):
        pass

    def ExecuteFinalizeSolutionStep(self):
        current_time = self.model_part.ProcessInfo[KratosMultiphysics.TIME] #- self.model_part.ProcessInfo[KratosMultiphysics.DELTA_TIME]

        if(current_time >= self.interval[0] and  current_time<self.interval[1]):

            if self.value_is_numeric:
                for node in self.mesh.Nodes:
                    value = node.GetSolutionStepValue(self.variable, 0)
                    self.assertAlmostEqual(self.value, value, self.tol)
            else:
                if self.is_time_function:
                    self.value = self.aux_function.f(0.0,0.0,0.0,current_time)
                    for node in self.mesh.Nodes:
                        value = node.GetSolutionStepValue(self.variable, 0)
                        self.assertAlmostEqual(self.value, value, self.tol)
                else: #most general case - space varying function (possibly also time varying)
                    if self.non_trivial_local_system == False:
                        for node in self.model_part.Nodes:
                            value = node.GetSolutionStepValue(self.variable, 0)
                            self.assertAlmostEqual(self.aux_function.f(node.X,node.Y,node.Z,current_time), value, self.tol)
                    else: #TODO: OPTIMIZE!!
                        for node in self.model_part.Nodes:
                            dx = node - self.x0
                            x = self.R[0,0]*dx[0] + self.R[0,1]*dx[1] + self.R[0,2]*dx[2]
                            y = self.R[1,0]*dx[0] + self.R[1,1]*dx[1] + self.R[1,2]*dx[2]
                            z = self.R[2,0]*dx[0] + self.R[2,1]*dx[1] + self.R[2,2]*dx[2]
                            value = node.GetSolutionStepValue(self.variable, 0)
                            self.assertAlmostEqual(self.aux_function.f(x,y,z,current_time), value, self.tol)

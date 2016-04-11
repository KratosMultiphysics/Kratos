import KratosMultiphysics 
from KratosMultiphysics.SolidMechanicsApplication import *

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplyPointLoadProcess(Model, settings["parameters"])

##all the processes python processes should be derived from "python_process"
class ApplyPointLoadProcess(KratosMultiphysics.ApplyConstantVectorValueProcess):
    def __init__(self, Model, Parameters ):
        model_part = Model[Parameters["model_part_name"].GetString()]
        Parameters.AddEmptyValue("variable_name").SetString("POINT_LOAD")
        Parameters.AddEmptyValue("is_fixed_x").SetBool(False)
        Parameters.AddEmptyValue("is_fixed_y").SetBool(False)
        Parameters.AddEmptyValue("is_fixed_z").SetBool(False)
               
        KratosMultiphysics.ApplyConstantVectorValueProcess.__init__(self,model_part, Parameters)

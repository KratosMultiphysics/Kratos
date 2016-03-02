import KratosMultiphysics 
import KratosMultiphysics.FluidDynamicsApplication 

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplyConstantVectorValue(Model, settings["Parameters"])

##all the processes python processes should be derived from "python_process"
class ApplyNoSlipProcess(KratosMultiphysics.ApplyConstantVectorValueProcess):
    def __init__(self, Model, Parameters ):
        model_part = Model[Parameters["model_part_name"].GetString()]
        Parameters["variable_name"].SetString("VELOCITY")
        Parameters["is_fixed_x"].SetBool(True)
        Parameters["is_fixed_y"].SetBool(True)
        Parameters["is_fixed_z"].SetBool(True)
        KratosMultiphysics.ApplyConstantVectorValueProcess(model_part, Parameters)
        

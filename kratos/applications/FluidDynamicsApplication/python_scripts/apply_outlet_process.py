import KratosMultiphysics 
import KratosMultiphysics.FluidDynamicsApplication as cfd

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
        
    return ApplyOutletProcess(Model, settings["Parameters"])

##all the processes python processes should be derived from "python_process"
class ApplyOutletProcess(KratosMultiphysics.ApplyConstantScalarValueProcess):
    def __init__(self, Model, Parameters ):
        model_part = Model[Parameters["model_part_name"].GetString()]
        Parameters.AddEmptyValue("variable_name").SetString("EXTERNAL_PRESSURE")
        KratosMultiphysics.ApplyConstantScalarValueProcess(model_part, Parameters)
        

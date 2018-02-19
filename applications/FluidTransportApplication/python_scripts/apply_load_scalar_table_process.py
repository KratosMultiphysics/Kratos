import KratosMultiphysics
import KratosMultiphysics.PoromechanicsApplication as KratosPoro

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplyLoadScalarTableProcess(Model, settings["Parameters"])

## All the python processes should be derived from "python_process"

class ApplyLoadScalarTableProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings ):
        KratosMultiphysics.Process.__init__(self)

        model_part = Model[settings["model_part_name"].GetString()]

        params = KratosMultiphysics.Parameters("{}")
        params.AddValue("model_part_name",settings["model_part_name"])
        params.AddValue("mesh_id",settings["mesh_id"])
        params.AddValue("variable_name",settings["variable_name"])
        params.AddValue("value",settings["value"])
        
        if settings["table"].GetInt() == 0:
            self.process = KratosMultiphysics.ApplyConstantScalarValueProcess(model_part, params)
        else:
            params.AddValue("table",settings["table"])
            self.process = KratosPoro.ApplyDoubleTableProcess(model_part, params)

    def ExecuteInitialize(self):
        
        self.process.ExecuteInitialize()

    def ExecuteInitializeSolutionStep(self):
        
        self.process.ExecuteInitializeSolutionStep()
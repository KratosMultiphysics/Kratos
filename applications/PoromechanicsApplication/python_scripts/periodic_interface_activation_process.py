import KratosMultiphysics
import KratosMultiphysics.PoromechanicsApplication as KratosPoro

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return PeriodicInterfaceActivationProcess(Model, settings["Parameters"])

## All the python processes should be derived from "python_process"

class PeriodicInterfaceActivationProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings ):
        KratosMultiphysics.Process.__init__(self)

        model_part = Model[settings["model_part_name"].GetString()]

        params = KratosMultiphysics.Parameters("{}")
        params.AddValue("model_part_name",settings["model_part_name"])
        params.AddValue("mesh_id",settings["mesh_id"])
        params.AddValue("dimension",settings["dimension"])
        params.AddValue("von_mises_limit",settings["von_mises_limit"])
        
        self.process = KratosPoro.PeriodicInterfaceProcess(model_part, params)

        if(settings["dimension"].GetInt() == 2):
            self.FindNodalNeigh = KratosMultiphysics.FindNodalNeighboursProcess(model_part,2,5)
        else:
            self.FindNodalNeigh = KratosMultiphysics.FindNodalNeighboursProcess(model_part,10,10)

    def ExecuteInitialize(self):
        
        self.FindNodalNeigh.Execute()

        self.process.ExecuteInitialize()

    def ExecuteFinalizeSolutionStep(self):

        self.process.ExecuteFinalizeSolutionStep()

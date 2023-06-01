import KratosMultiphysics
import KratosMultiphysics.GeoMechanicsApplication as KratosGeo

def Factory(settings, Model):
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise TypeError("expected input shall be a Parameters object, encapsulating a json string")
    return PeriodicInterfaceActivationProcess(Model, settings["Parameters"])

## All the python processes should be derived from "python_process"

class PeriodicInterfaceActivationProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings ):
        KratosMultiphysics.Process.__init__(self)

        model_part = Model[settings["model_part_name"].GetString()]

        params = KratosMultiphysics.Parameters("{}")
        params.AddValue("model_part_name",settings["model_part_name"])
        params.AddValue("dimension",settings["dimension"])
        params.AddValue("stress_limit",settings["stress_limit"])

        self.process = KratosGeo.PeriodicInterfaceProcess(model_part, params)

        if(settings["dimension"].GetInt() == 2):
            self.FindNodalNeigh = KratosMultiphysics.FindNodalNeighboursProcess(model_part,2,5)
        else:
            self.FindNodalNeigh = KratosMultiphysics.FindNodalNeighboursProcess(model_part,10,10)

    def ExecuteInitialize(self):

        self.FindNodalNeigh.Execute()

        self.process.ExecuteInitialize()

    def ExecuteFinalizeSolutionStep(self):

        self.process.ExecuteFinalizeSolutionStep()

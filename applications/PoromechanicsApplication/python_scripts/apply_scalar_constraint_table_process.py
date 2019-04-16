import KratosMultiphysics
import KratosMultiphysics.PoromechanicsApplication as KratosPoro

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplyScalarConstraintTableProcess(Model, settings["Parameters"])

## All the processes python should be derived from "Process"

class ApplyScalarConstraintTableProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings ):
        KratosMultiphysics.Process.__init__(self)

        model_part = Model[settings["model_part_name"].GetString()]

        params = KratosMultiphysics.Parameters("{}")
        params.AddValue("model_part_name",settings["model_part_name"])
        params.AddValue("variable_name",settings["variable_name"])
        if settings.Has("is_fixed"):
            params.AddValue("is_fixed",settings["is_fixed"])
        
        if settings.Has("hydrostatic"):
            if settings["hydrostatic"].GetBool() == False:
                params.AddValue("value",settings["value"])
                if settings["table"].GetInt() == 0:
                    self.process = KratosMultiphysics.ApplyConstantScalarValueProcess(model_part, params)
                else:
                    params.AddValue("table",settings["table"])
                    self.process = KratosPoro.ApplyDoubleTableProcess(model_part, params)
            else:
                params.AddValue("gravity_direction",settings["gravity_direction"])
                params.AddValue("reference_coordinate",settings["reference_coordinate"])
                params.AddValue("specific_weight",settings["specific_weight"])
                if settings["table"].GetInt() == 0:
                    self.process = KratosPoro.ApplyConstantHydrostaticPressureProcess(model_part, params)
                else:
                    params.AddValue("table",settings["table"])
                    self.process = KratosPoro.ApplyHydrostaticPressureTableProcess(model_part, params)
        else:
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
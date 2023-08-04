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

        self.settings = settings
        self.interval = KratosMultiphysics.IntervalUtility(self.settings)

        self.model_part = Model[settings["model_part_name"].GetString()]

        self.params = KratosMultiphysics.Parameters("{}")
        self.params.AddValue("model_part_name",settings["model_part_name"])
        self.params.AddValue("variable_name",settings["variable_name"])
        if settings.Has("is_fixed"):
            self.params.AddValue("is_fixed",settings["is_fixed"])

        if settings.Has("hydrostatic"):
            if settings["hydrostatic"].GetBool() == False:
                self.params.AddValue("value",settings["value"])
                if settings["table"].GetInt() == 0:
                    self.process = KratosMultiphysics.ApplyConstantScalarValueProcess(self.model_part, self.params)
                else:
                    self.params.AddValue("table",settings["table"])
                    self.process = KratosPoro.ApplyDoubleTableProcess(self.model_part, self.params)
            else:
                self.params.AddValue("gravity_direction",settings["gravity_direction"])
                self.params.AddValue("reference_coordinate",settings["reference_coordinate"])
                self.params.AddValue("specific_weight",settings["specific_weight"])
                if settings["table"].GetInt() == 0:
                    self.process = KratosPoro.ApplyConstantHydrostaticPressureProcess(self.model_part, self.params)
                else:
                    self.params.AddValue("table",settings["table"])
                    self.process = KratosPoro.ApplyHydrostaticPressureTableProcess(self.model_part, self.params)
        else:
            self.params.AddValue("value",settings["value"])
            if settings["table"].GetInt() == 0:
                self.process = KratosMultiphysics.ApplyConstantScalarValueProcess(self.model_part, self.params)
            else:
                self.params.AddValue("table",settings["table"])
                self.process = KratosPoro.ApplyDoubleTableProcess(self.model_part, self.params)

    def ExecuteInitialize(self):

        current_time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]

        if self.interval.IsInInterval(current_time):
            self.process.ExecuteInitialize()

    def ExecuteInitializeSolutionStep(self):

        current_time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]

        if self.interval.IsInInterval(current_time):
            self.process.ExecuteInitializeSolutionStep()

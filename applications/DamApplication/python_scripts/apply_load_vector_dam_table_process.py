from KratosMultiphysics import *
from KratosMultiphysics.DamApplication import *
from KratosMultiphysics.PoromechanicsApplication import *

def Factory(settings, Model):
    if not isinstance(settings, Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplyLoadVectorDamTableProcess(Model, settings["Parameters"])

## All the processes python should be derived from "Process"

class ApplyLoadVectorDamTableProcess(Process):
    def __init__(self, Model, settings ):
        Process.__init__(self)

        model_part = Model[settings["model_part_name"].GetString()]
        variable_name = settings["variable_name"].GetString()

        self.components_process_list = []

        self.factor = settings["modulus"].GetDouble();
        self.direction = [settings["direction"][0].GetDouble(),settings["direction"][1].GetDouble(),settings["direction"][2].GetDouble()]
        self.value = [self.direction[0]*self.factor,self.direction[1]*self.factor,self.direction[2]*self.factor]

        if abs(self.value[0])>1.0e-15:
            x_params = Parameters("{}")
            x_params.AddValue("model_part_name",settings["model_part_name"])
            x_params.AddEmptyValue("value").SetDouble(self.value[0])
            x_params.AddEmptyValue("variable_name").SetString(variable_name+"_X")
            if settings["table"].GetInt() == 0:
                self.components_process_list.append(ApplyConstantScalarValueProcess(model_part, x_params))
            else:
                x_params.AddValue("table",settings["table"])
                self.components_process_list.append(ApplyComponentTableProcess(model_part, x_params))

        if abs(self.value[1])>1.0e-15:
            y_params = Parameters("{}")
            y_params.AddValue("model_part_name",settings["model_part_name"])
            y_params.AddEmptyValue("value").SetDouble(self.value[1])
            y_params.AddEmptyValue("variable_name").SetString(variable_name+"_Y")
            if settings["table"].GetInt() == 0:
                self.components_process_list.append(ApplyConstantScalarValueProcess(model_part, y_params))
            else:
                y_params.AddValue("table",settings["table"])
                self.components_process_list.append(ApplyComponentTableProcess(model_part, y_params))

        if abs(self.value[2])>1.0e-15:
            z_params = Parameters("{}")
            z_params.AddValue("model_part_name",settings["model_part_name"])
            z_params.AddEmptyValue("value").SetDouble(self.value[2])
            z_params.AddEmptyValue("variable_name").SetString(variable_name+"_Z")
            if settings["table"].GetInt() == 0:
                self.components_process_list.append(ApplyConstantScalarValueProcess(model_part, z_params))
            else:
                z_params.AddValue("table",settings["table"])
                self.components_process_list.append(ApplyComponentTableProcess(model_part, z_params))

    def ExecuteInitialize(self):

        for component in self.components_process_list:
            component.ExecuteInitialize()

    def ExecuteInitializeSolutionStep(self):

        for component in self.components_process_list:
            component.ExecuteInitializeSolutionStep()

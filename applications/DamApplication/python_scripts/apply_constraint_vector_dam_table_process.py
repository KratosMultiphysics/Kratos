import KratosMultiphysics
import KratosMultiphysics.DamApplication as KratosDam
from KratosMultiphysics.assign_scalar_variable_process import AssignScalarVariableProcess

def Factory(settings, Model):
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplyConstraintVectorDamTableProcess(Model, settings["Parameters"])

## All the processes python should be derived from "Process"

class ApplyConstraintVectorDamTableProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings ):
        KratosMultiphysics.Process.__init__(self)

        model_part = Model[settings["model_part_name"].GetString()]
        variable_name = settings["variable_name"].GetString()

        self.components_process_list = []

        x_params = KratosMultiphysics.Parameters("{}")
        x_params.AddValue("model_part_name",settings["model_part_name"])
        x_params.AddValue("constrained",settings["constrained"][0])
        x_params.AddValue("value",settings["value"][0])
        x_params.AddEmptyValue("variable_name").SetString(variable_name+"_X")
        if settings["Value_Table"][0].GetInt() == 0:
            self.components_process_list.append(AssignScalarVariableProcess(Model, x_params))
        else:
            x_params.AddValue("table",settings["Value_Table"][0])
            self.components_process_list.append(KratosDam.ApplyComponentTableProcessDam(model_part, x_params))

        y_params = KratosMultiphysics.Parameters("{}")
        y_params.AddValue("model_part_name",settings["model_part_name"])
        y_params.AddValue("constrained",settings["constrained"][1])
        y_params.AddValue("value",settings["value"][1])
        y_params.AddEmptyValue("variable_name").SetString(variable_name+"_Y")
        if settings["Value_Table"][1].GetInt() == 0:
            self.components_process_list.append(AssignScalarVariableProcess(Model, y_params))
        else:
            y_params.AddValue("table",settings["Value_Table"][1])
            self.components_process_list.append(KratosDam.ApplyComponentTableProcessDam(model_part, y_params))

        z_params = KratosMultiphysics.Parameters("{}")
        z_params.AddValue("model_part_name",settings["model_part_name"])
        z_params.AddValue("constrained",settings["constrained"][2])
        z_params.AddValue("value",settings["value"][2])
        z_params.AddEmptyValue("variable_name").SetString(variable_name+"_Z")
        if settings["Value_Table"][2].GetInt() == 0:
            self.components_process_list.append(AssignScalarVariableProcess(Model, z_params))
        else:
            z_params.AddValue("table",settings["Value_Table"][2])
            self.components_process_list.append(KratosDam.ApplyComponentTableProcessDam(model_part, z_params))

    def ExecuteBeforeSolutionLoop(self):

        for component in self.components_process_list:
            component.ExecuteBeforeSolutionLoop()

    def ExecuteInitializeSolutionStep(self):

        for component in self.components_process_list:
            component.ExecuteInitializeSolutionStep()

import KratosMultiphysics
import KratosMultiphysics.PoromechanicsApplication as KratosPoro

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplyVectorConstraintTableProcess(Model, settings["Parameters"])

## All the processes python should be derived from "Process"

class ApplyVectorConstraintTableProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings ):
        KratosMultiphysics.Process.__init__(self)

        self.model_part = Model[settings["model_part_name"].GetString()]
        variable_name = settings["variable_name"].GetString()

        self.components_process_list = []

        if settings["active"][0].GetBool() == True:
            self.x_params = KratosMultiphysics.Parameters("{}")
            self.x_params.AddValue("model_part_name",settings["model_part_name"])
            if settings.Has("is_fixed"):
                self.x_params.AddValue("is_fixed",settings["is_fixed"][0])
            self.x_params.AddValue("value",settings["value"][0])
            self.x_params.AddEmptyValue("variable_name").SetString(variable_name+"_X")
            if settings["table"][0].GetInt() == 0:
                self.components_process_list.append(KratosMultiphysics.ApplyConstantScalarValueProcess(self.model_part, self.x_params))
            else:
                self.x_params.AddValue("table",settings["table"][0])
                self.components_process_list.append(KratosPoro.ApplyComponentTableProcess(self.model_part, self.x_params))

        if settings["active"][1].GetBool() == True:
            self.y_params = KratosMultiphysics.Parameters("{}")
            self.y_params.AddValue("model_part_name",settings["model_part_name"])
            if settings.Has("is_fixed"):
                self.y_params.AddValue("is_fixed",settings["is_fixed"][1])
            self.y_params.AddValue("value",settings["value"][1])
            self.y_params.AddEmptyValue("variable_name").SetString(variable_name+"_Y")
            if settings["table"][1].GetInt() == 0:
                self.components_process_list.append(KratosMultiphysics.ApplyConstantScalarValueProcess(self.model_part, self.y_params))
            else:
                self.y_params.AddValue("table",settings["table"][1])
                self.components_process_list.append(KratosPoro.ApplyComponentTableProcess(self.model_part, self.y_params))

        if settings["active"][2].GetBool() == True:
            self.z_params = KratosMultiphysics.Parameters("{}")
            self.z_params.AddValue("model_part_name",settings["model_part_name"])
            if settings.Has("is_fixed"):
                self.z_params.AddValue("is_fixed",settings["is_fixed"][2])
            self.z_params.AddValue("value",settings["value"][2])
            self.z_params.AddEmptyValue("variable_name").SetString(variable_name+"_Z")
            if settings["table"][2].GetInt() == 0:
                self.components_process_list.append(KratosMultiphysics.ApplyConstantScalarValueProcess(self.model_part, self.z_params))
            else:
                self.z_params.AddValue("table",settings["table"][2])
                self.components_process_list.append(KratosPoro.ApplyComponentTableProcess(self.model_part, self.z_params))

    def ExecuteInitialize(self):

        for component in self.components_process_list:
            component.ExecuteInitialize()

    def ExecuteInitializeSolutionStep(self):

        for component in self.components_process_list:
            component.ExecuteInitializeSolutionStep()
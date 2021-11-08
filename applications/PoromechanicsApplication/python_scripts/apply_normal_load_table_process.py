import KratosMultiphysics
import KratosMultiphysics.PoromechanicsApplication as KratosPoro

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplyNormalLoadTableProcess(Model, settings["Parameters"])

## All the processes python should be derived from "Process"

class ApplyNormalLoadTableProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings ):
        KratosMultiphysics.Process.__init__(self)

        self.model_part = Model[settings["model_part_name"].GetString()]

        self.components_process_list = []

        if settings["active"][0].GetBool() == True:
            self.normal_params = KratosMultiphysics.Parameters("{}")
            self.normal_params.AddValue("model_part_name",settings["model_part_name"])
            self.normal_params.AddValue("variable_name",settings["variable_name"])
            if settings["hydrostatic"].GetBool() == False:
                self.normal_params.AddValue("value",settings["value"][0])
                if settings["table"][0].GetInt() == 0:
                    self.components_process_list.append(KratosMultiphysics.ApplyConstantScalarValueProcess(self.model_part, self.normal_params))
                else:
                    self.normal_params.AddValue("table",settings["table"][0])
                    self.components_process_list.append(KratosPoro.ApplyDoubleTableProcess(self.model_part, self.normal_params))
            else:
                self.normal_params.AddValue("gravity_direction",settings["gravity_direction"])
                self.normal_params.AddValue("reference_coordinate",settings["reference_coordinate"])
                self.normal_params.AddValue("specific_weight",settings["specific_weight"])
                if settings["table"][0].GetInt() == 0:
                    self.components_process_list.append(KratosPoro.ApplyConstantHydrostaticPressureProcess(self.model_part, self.normal_params))
                else:
                    self.normal_params.AddValue("table",settings["table"][0])
                    self.components_process_list.append(KratosPoro.ApplyHydrostaticPressureTableProcess(self.model_part, self.normal_params))

        if settings["active"][1].GetBool() == True:
            tangential_params = KratosMultiphysics.Parameters("{}")
            tangential_params.AddValue("model_part_name",settings["model_part_name"])
            tangential_params.AddEmptyValue("variable_name").SetString("TANGENTIAL_CONTACT_STRESS") # Note: this is not general
            tangential_params.AddValue("value",settings["value"][1])
            if settings["table"][1].GetInt() == 0:
                self.components_process_list.append(KratosMultiphysics.ApplyConstantScalarValueProcess(self.model_part, tangential_params))
            else:
                tangential_params.AddValue("table",settings["table"][1])
                self.components_process_list.append(KratosPoro.ApplyDoubleTableProcess(self.model_part, tangential_params))

    def ExecuteInitialize(self):

        for component in self.components_process_list:
            component.ExecuteInitialize()

    def ExecuteInitializeSolutionStep(self):

        for component in self.components_process_list:
            component.ExecuteInitializeSolutionStep()
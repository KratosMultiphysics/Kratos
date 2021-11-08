import KratosMultiphysics
import KratosMultiphysics.FemToDemApplication as KratosFemDem

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplyNormalLoadTableProcess(Model, settings["Parameters"])

## All the processes python should be derived from "Process"
class ApplyNormalLoadTableProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings ):
        KratosMultiphysics.Process.__init__(self)

        model_part = Model[settings["model_part_name"].GetString()]
        
        self.components_process_list = []
        
        if settings["active"][0].GetBool() == True:
            normal_params = KratosMultiphysics.Parameters("{}")
            normal_params.AddValue("model_part_name",settings["model_part_name"])
            normal_params.AddValue("variable_name",settings["variable_name"])

            normal_params.AddValue("value",settings["value"][0])
            if settings["table"][0].GetInt() == 0:
                self.components_process_list.append(KratosMultiphysics.ApplyConstantScalarValueProcess(model_part, normal_params))
            else:
                normal_params.AddValue("table",settings["table"][0])
                self.components_process_list.append(KratosFemDem.ApplyDoubleTableProcess(model_part, normal_params))

        if settings["active"][1].GetBool() == True:
            tangential_params = KratosMultiphysics.Parameters("{}")
            tangential_params.AddValue("model_part_name",settings["model_part_name"])
            tangential_params.AddEmptyValue("variable_name").SetString("TANGENTIAL_CONTACT_STRESS") # Note: this is not general
            tangential_params.AddValue("value",settings["value"][1])
            if settings["table"][1].GetInt() == 0:
                self.components_process_list.append(KratosMultiphysics.ApplyConstantScalarValueProcess(model_part, tangential_params))
            else:
                tangential_params.AddValue("table",settings["table"][1])
                self.components_process_list.append(KratosFemDem.ApplyDoubleTableProcess(model_part, tangential_params))
                
    def ExecuteInitialize(self):
        for component in self.components_process_list:
            component.ExecuteInitialize()

    def ExecuteInitializeSolutionStep(self):
        for component in self.components_process_list:
            component.ExecuteInitializeSolutionStep()
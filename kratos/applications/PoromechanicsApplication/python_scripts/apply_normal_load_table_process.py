from KratosMultiphysics import *
from KratosMultiphysics.PoromechanicsApplication import *


def Factory(settings, Model):
    if(type(settings) != Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplyNormalLoadTableProcess(Model, settings["Parameters"])

## All the python processes should be derived from "python_process"

class ApplyNormalLoadTableProcess(Process):
    def __init__(self, Model, settings ):
        Process.__init__(self)

        model_part = Model[settings["model_part_name"].GetString()]
        
        self.components_process_list = []
        
        if settings["active"][0].GetBool() == True:
            normal_params = Parameters("{}")
            normal_params.AddValue("model_part_name",settings["model_part_name"])
            normal_params.AddValue("mesh_id",settings["mesh_id"])
            normal_params.AddValue("variable_name",settings["variable_name"])
            if settings["hydrostatic"].GetBool() == False:
                normal_params.AddValue("value",settings["value"][0])
                if settings["table"][0].GetInt() == 0:
                    self.components_process_list.append(ApplyConstantScalarValueProcess(model_part, normal_params))
                else:
                    normal_params.AddValue("table",settings["table"][0])
                    self.components_process_list.append(ApplyDoubleTableProcess(model_part, normal_params))
            else:
                normal_params.AddValue("gravity_direction",settings["gravity_direction"])
                normal_params.AddValue("reference_coordinate",settings["reference_coordinate"])
                normal_params.AddValue("specific_weight",settings["specific_weight"])
                if settings["table"][0].GetInt() == 0:
                    self.components_process_list.append(ApplyConstantHydrostaticPressureProcess(model_part, normal_params))
                else:
                    normal_params.AddValue("table",settings["table"][0])
                    self.components_process_list.append(ApplyHydrostaticPressureTableProcess(model_part, normal_params))

        if settings["active"][1].GetBool() == True:
            tangential_params = Parameters("{}")
            tangential_params.AddValue("model_part_name",settings["model_part_name"])
            tangential_params.AddValue("mesh_id",settings["mesh_id"])
            tangential_params.AddEmptyValue("variable_name").SetString("TANGENTIAL_CONTACT_STRESS") # Note: this is not general
            tangential_params.AddValue("value",settings["value"][1])
            if settings["table"][1].GetInt() == 0:
                self.components_process_list.append(ApplyConstantScalarValueProcess(model_part, tangential_params))
            else:
                tangential_params.AddValue("table",settings["table"][1])
                self.components_process_list.append(ApplyDoubleTableProcess(model_part, tangential_params))
                
    def ExecuteInitialize(self):

        for component in self.components_process_list:
            component.ExecuteInitialize()

    def ExecuteInitializeSolutionStep(self):

        for component in self.components_process_list:
            component.ExecuteInitializeSolutionStep()
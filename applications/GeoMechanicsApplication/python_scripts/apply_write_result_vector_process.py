import KratosMultiphysics
import KratosMultiphysics.GeoMechanicsApplication as KratosGeo

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplyWriteVectorProcess(Model, settings["Parameters"])

## All the python processes should be derived from "python_process"

class ApplyWriteVectorProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings ):
        KratosMultiphysics.Process.__init__(self)

        self.model_part = Model[settings["model_part_name"].GetString()]
        variable_name = settings["variable_name"].GetString()

        self.components_process_list = []

        if settings["active"][0].GetBool() == True:
            self.x_params = KratosMultiphysics.Parameters("{}")
            self.x_params.AddValue("model_part_name",settings["model_part_name"])
            self.x_params.AddValue("append_file",settings["append_file"])

            self.x_params.AddEmptyValue("variable_name").SetString(variable_name+"_X")
            self.components_process_list.append(KratosGeo.ApplyWriteScalarProcess(self.model_part, self.x_params))

        if settings["active"][1].GetBool() == True:
            self.y_params = KratosMultiphysics.Parameters("{}")
            self.y_params.AddValue("model_part_name",settings["model_part_name"])
            self.y_params.AddValue("append_file",settings["append_file"])

            self.y_params.AddEmptyValue("variable_name").SetString(variable_name+"_Y")
            self.components_process_list.append(KratosGeo.ApplyWriteScalarProcess(self.model_part, self.y_params))

        if settings["active"][2].GetBool() == True:
            self.z_params = KratosMultiphysics.Parameters("{}")
            self.z_params.AddValue("model_part_name",settings["model_part_name"])
            self.z_params.AddValue("append_file",settings["append_file"])

            self.z_params.AddEmptyValue("variable_name").SetString(variable_name+"_Z")
            self.components_process_list.append(KratosGeo.ApplyWriteScalarProcess(self.model_part, self.z_params))

    def ExecuteInitialize(self):

        for component in self.components_process_list:
            component.ExecuteInitialize()

    def ExecuteInitializeSolutionStep(self):

        for component in self.components_process_list:
            component.ExecuteInitializeSolutionStep()

    def ExecuteFinalizeSolutionStep(self):

        for component in self.components_process_list:
            component.ExecuteFinalizeSolutionStep()

    def ExecuteFinalize(self):

        for component in self.components_process_list:
            component.ExecuteFinalize()

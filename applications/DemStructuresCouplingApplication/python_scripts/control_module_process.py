import KratosMultiphysics
import KratosMultiphysics.DemStructuresCouplingApplication as DemFem

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ControlModuleProcess(Model, settings["Parameters"])

## All the processes python should be derived from "Process"

class ControlModuleProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings ):
        KratosMultiphysics.Process.__init__(self)

        model_part = Model[settings["model_part_name"].GetString()]
        variable_name = settings["variable_name"].GetString()
        reaction_variable_name = settings["reaction_variable_name"].GetString()
        target_stress_variable_name = settings["target_stress_variable_name"].GetString()
        reaction_stress_variable_name = settings["reaction_stress_variable_name"].GetString()
        loading_velocity_variable_name = settings["loading_velocity_variable_name"].GetString()

        self.components_process_list = []

        if settings["fixed"][0].GetBool() == True:
            x_params = KratosMultiphysics.Parameters("{}")
            x_params.AddValue("model_part_name",settings["model_part_name"])
            x_params.AddEmptyValue("variable_name").SetString(variable_name+"_X")
            x_params.AddEmptyValue("reaction_variable_name").SetString(reaction_variable_name+"_X")
            x_params.AddEmptyValue("target_stress_variable_name").SetString(target_stress_variable_name+"_X")
            x_params.AddEmptyValue("reaction_stress_variable_name").SetString(reaction_stress_variable_name+"_X")
            x_params.AddEmptyValue("loading_velocity_variable_name").SetString(loading_velocity_variable_name+"_X")
            x_params.AddValue("target_stress_table_id",settings["target_stress_table_id"][0])
            x_params.AddValue("initial_velocity",settings["initial_velocity"][0])
            x_params.AddValue("limit_velocity",settings["limit_velocity"][0])
            x_params.AddValue("velocity_factor",settings["velocity_factor"])
            x_params.AddValue("compression_length",settings["compression_length"])
            x_params.AddValue("young_modulus",settings["young_modulus"])
            x_params.AddValue("stress_increment_tolerance",settings["stress_increment_tolerance"])
            x_params.AddValue("update_stiffness",settings["update_stiffness"])
            x_params.AddValue("start_time",settings["start_time"])
            self.components_process_list.append(DemFem.ControlModuleProcess(model_part, x_params))

        if settings["fixed"][1].GetBool() == True:
            y_params = KratosMultiphysics.Parameters("{}")
            y_params.AddValue("model_part_name",settings["model_part_name"])
            y_params.AddEmptyValue("variable_name").SetString(variable_name+"_Y")
            y_params.AddEmptyValue("reaction_variable_name").SetString(reaction_variable_name+"_Y")
            y_params.AddEmptyValue("target_stress_variable_name").SetString(target_stress_variable_name+"_Y")
            y_params.AddEmptyValue("reaction_stress_variable_name").SetString(reaction_stress_variable_name+"_Y")
            y_params.AddEmptyValue("loading_velocity_variable_name").SetString(loading_velocity_variable_name+"_Y")
            y_params.AddValue("target_stress_table_id",settings["target_stress_table_id"][1])
            y_params.AddValue("initial_velocity",settings["initial_velocity"][1])
            y_params.AddValue("limit_velocity",settings["limit_velocity"][1])
            y_params.AddValue("velocity_factor",settings["velocity_factor"])
            y_params.AddValue("compression_length",settings["compression_length"])
            y_params.AddValue("young_modulus",settings["young_modulus"])
            y_params.AddValue("stress_increment_tolerance",settings["stress_increment_tolerance"])
            y_params.AddValue("update_stiffness",settings["update_stiffness"])
            y_params.AddValue("start_time",settings["start_time"])
            self.components_process_list.append(DemFem.ControlModuleProcess(model_part, y_params))

        if settings["fixed"][2].GetBool() == True:
            z_params = KratosMultiphysics.Parameters("{}")
            z_params.AddValue("model_part_name",settings["model_part_name"])
            z_params.AddEmptyValue("variable_name").SetString(variable_name+"_Z")
            z_params.AddEmptyValue("reaction_variable_name").SetString(reaction_variable_name+"_Z")
            z_params.AddEmptyValue("target_stress_variable_name").SetString(target_stress_variable_name+"_Z")
            z_params.AddEmptyValue("reaction_stress_variable_name").SetString(reaction_stress_variable_name+"_Z")
            z_params.AddEmptyValue("loading_velocity_variable_name").SetString(loading_velocity_variable_name+"_Z")
            z_params.AddValue("target_stress_table_id",settings["target_stress_table_id"][2])
            z_params.AddValue("initial_velocity",settings["initial_velocity"][2])
            z_params.AddValue("limit_velocity",settings["limit_velocity"][2])
            z_params.AddValue("velocity_factor",settings["velocity_factor"])
            z_params.AddValue("compression_length",settings["compression_length"])
            z_params.AddValue("young_modulus",settings["young_modulus"])
            z_params.AddValue("stress_increment_tolerance",settings["stress_increment_tolerance"])
            z_params.AddValue("update_stiffness",settings["update_stiffness"])
            z_params.AddValue("start_time",settings["start_time"])
            self.components_process_list.append(DemFem.ControlModuleProcess(model_part, z_params))

    def ExecuteInitialize(self):

        for component in self.components_process_list:
            component.ExecuteInitialize()

    def ExecuteInitializeSolutionStep(self):

        for component in self.components_process_list:
            component.ExecuteInitializeSolutionStep()

    def ExecuteFinalizeSolutionStep(self):

        for component in self.components_process_list:
            component.ExecuteFinalizeSolutionStep()
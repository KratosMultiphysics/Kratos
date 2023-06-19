import KratosMultiphysics
import KratosMultiphysics.PoromechanicsApplication as KratosPoro

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return PoromechanicsFaceLoadControlModuleProcess(Model, settings["Parameters"])

## All the processes python should be derived from "Process"

class PoromechanicsFaceLoadControlModuleProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings ):
        KratosMultiphysics.Process.__init__(self)

        # Control module process acting on the imposed direction: 0 (X), 1 (Y), 2 (Z)

        self.model_part = Model[settings["model_part_name"].GetString()]
        self.components_process_list = []

        # Settings string in json format
        default_parameters = KratosMultiphysics.Parameters("{}")
        default_parameters.AddValue("model_part_name",settings["model_part_name"])
        default_parameters.AddValue("initial_velocity",settings["initial_velocity"])
        default_parameters.AddValue("limit_velocity",settings["limit_velocity"])
        default_parameters.AddValue("velocity_factor",settings["velocity_factor"])
        default_parameters.AddValue("initial_stiffness",settings["initial_stiffness"])
        default_parameters.AddValue("force_increment_tolerance",settings["force_increment_tolerance"])
        default_parameters.AddValue("update_stiffness",settings["update_stiffness"])
        default_parameters.AddValue("force_averaging_time",settings["force_averaging_time"])

        if settings["active"][0].GetBool() == True:
            self.x_params = default_parameters
            self.x_params.AddEmptyValue("imposed_direction")
            self.x_params["imposed_direction"].SetInt(0)
            self.x_params.AddValue("face_load_table_id",settings["table"][0])
            self.components_process_list.append(KratosPoro.PoromechanicsFaceLoadControlModuleProcess(self.model_part, self.x_params))
        if settings["active"][1].GetBool() == True:
            self.y_params = default_parameters
            self.y_params.AddEmptyValue("imposed_direction")
            self.y_params["imposed_direction"].SetInt(1)
            self.y_params.AddValue("face_load_table_id",settings["table"][1])
            self.components_process_list.append(KratosPoro.PoromechanicsFaceLoadControlModuleProcess(self.model_part, self.y_params))
        if settings["active"][2].GetBool() == True:
            self.z_params = default_parameters
            self.z_params.AddEmptyValue("imposed_direction")
            self.z_params["imposed_direction"].SetInt(2)
            self.z_params.AddValue("face_load_table_id",settings["table"][2])
            self.components_process_list.append(KratosPoro.PoromechanicsFaceLoadControlModuleProcess(self.model_part, self.z_params))

    def ExecuteInitialize(self):
        for component in self.components_process_list:
            component.ExecuteInitialize()

    def ExecuteInitializeSolutionStep(self):
        for component in self.components_process_list:
            component.ExecuteInitializeSolutionStep()

    def ExecuteFinalizeSolutionStep(self):
        for component in self.components_process_list:
            component.ExecuteFinalizeSolutionStep()
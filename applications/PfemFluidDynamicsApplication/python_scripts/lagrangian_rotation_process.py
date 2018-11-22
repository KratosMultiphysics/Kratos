import KratosMultiphysics
import KratosMultiphysics.PoromechanicsApplication as KratosPoro

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return LagrangianRotationProcess(Model, settings["Parameters"])

## All the processes python should be derived from "Process"

class LagrangianRotationProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings ):
        KratosMultiphysics.Process.__init__(self)

        model_part = Model[settings["model_part_name"].GetString()]
        variable_name = settings["variable_name"].GetString()

        params = KratosMultiphysics.Parameters("{}")
        params.AddValue("model_part_name",settings["model_part_name"])
        params.AddEmptyValue("variable_name_X").SetString(variable_name+"_X")
        params.AddEmptyValue("variable_name_Y").SetString(variable_name+"_Y")
        params.AddEmptyValue("variable_name_Z").SetString(variable_name+"_Z")
        if settings.Has("is_fixed"):
            params.AddValue("is_fixed",settings["is_fixed"][0])

            "model_part_name":"PLEASE_CHOOSE_MODEL_PART_NAME",
            "variable_name": "VELOCITY",
            "is_fixed": false,
            "angular_velocity": 0.0,
            "rotation_axis_initial_point": [0.0,0.0,0.0],
            "rotation_axis_final_point": [0.0,0.0,1.0],
            "initial_time": 0.0

        self.process.append(KratosMultiphysics.ApplyConstantScalarValueProcess(model_part, params))

    def ExecuteInitialize(self):

        self.process.ExecuteInitialize()

    def ExecuteInitializeSolutionStep(self):

        self.process.ExecuteInitializeSolutionStep()

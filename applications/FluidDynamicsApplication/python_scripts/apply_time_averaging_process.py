import KratosMultiphysics
import KratosMultiphysics.FluidDynamicsApplication as KratosFluid

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplyTimeAveragingProcess(Model, settings["Parameters"])

class ApplyTimeAveragingProcess(KratosMultiphysics.Process):

    def __init__(self, Model, settings):

        KratosMultiphysics.Process.__init__(self)

        default_parameters = KratosMultiphysics.Parameters("""
        {
            "mesh_id"                   : 0,
            "model_part_name"           : "",
            "variables"                 : ["VELOCITY", "PRESSURE", "REACTION"]
        }
        """)

        settings.ValidateAndAssignDefaults(default_parameters)

        # Get model part name
        self.fluid_model_part = Model[settings["model_part_name"].GetString()]
        # setting parameters to be averaged
        do_average_velocity = False
        do_average_pressure = False
        do_average_reaction = False

        for i_var in range(0, settings["variables"].size()):
            variable_name = settings["variables"][i_var].GetString()
            print(variable_name)
            if variable_name == "VELOCITY":
                do_average_velocity = True
            elif variable_name == "PRESSURE":
                do_average_pressure = True
            elif variable_name == "REACTION":
                do_average_reaction = True
            else:
                KratosMultiphysics.Logger.PrintWarning("ApplyTimeAveragingProcess", 
                                                        "You are averaging a variable which is not implemented in the averaging process." +
                                                        "Only VELOCITY, PRESSURE and Reaction will be averaged")
        self.TimeAveragingProcess = KratosFluid.TimeAveragingProcess(self.fluid_model_part, 
                                                                     do_average_velocity,
                                                                     do_average_pressure,
                                                                     do_average_reaction)

    def ExecuteInitialize(self):
        self.TimeAveragingProcess.ExecuteInitialize()

    def ExecuteFinalizeSolutionStep(self):
        self.TimeAveragingProcess.ExecuteFinalizeSolutionStep()



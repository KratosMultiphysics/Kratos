import KratosMultiphysics as KratosMultiphysics
import KratosMultiphysics.CableNetApplication as CableNetApplication

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return SlidingEdgeProcess(Model, settings["Parameters"])



class SlidingEdgeProcess(KratosMultiphysics.Process):

    def __init__(self, Model, settings ):
        KratosMultiphysics.Process.__init__(self)
        default_settings = KratosMultiphysics.Parameters("""
        {
            "constraint_name"           : "LinearMasterSlaveConstraint",
            "master_sub_model_part_name"    : "master_connect",
            "slave_sub_model_part_name"     : "slave_connect",
            "model_name"                    : "Structure",
            "variable_names"                : ["DISPLACEMENT_Y","DISPLACEMENT_Z"],
            "reform_every_step"             : true,
            "debug_info"                    : true,
            "angled_initial_line"           : false,
            "follow_line"                   : false
        }
        """)
        default_settings.ValidateAndAssignDefaults(settings)

        # The computing model part
        self.computing_model_part = Model["Structure"].GetSubModelPart("computing_domain")
        self.master_model_part = Model["Structure"].GetSubModelPart(settings["master_sub_model_part_name"].GetString())



        self.sliding_edge_process = CableNetApplication.SlidingEdgeProcess(Model[settings["model_name"].GetString()], settings)



    def ExecuteInitializeSolutionStep(self):
        self.sliding_edge_process.ExecuteInitializeSolutionStep()
        for constraint in self.master_model_part.MasterSlaveConstraints:
            self.computing_model_part.AddMasterSlaveConstraint(constraint)

    def ExecuteFinalizeSolutionStep(self):
        self.sliding_edge_process.ExecuteFinalizeSolutionStep()
import KratosMultiphysics
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD
from KratosMultiphysics.assign_vector_variable_process import AssignVectorVariableProcess

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplyFluidTopologyOptimizationNoSlipProcess(Model, settings["Parameters"])


class ApplyFluidTopologyOptimizationNoSlipProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings):

        KratosMultiphysics.Process.__init__(self)
        
        default_settings = KratosMultiphysics.Parameters("""
            {
                "mesh_id"         : 0,
                "model_part_name" : "please_specify_model_part_name",
                "interval"        : [0.0,1e30]
            }
            """)

        # Check input data
        if (settings["model_part_name"].GetString() == ""):
            raise Exception("Empty inlet model part name string. Set a valid model part name.")
        settings.ValidateAndAssignDefaults(default_settings)

        mesh_id = settings["mesh_id"].GetInt()
        model_part_name = settings["model_part_name"].GetString()
        interval = settings["interval"].PrettyPrintJsonString()

        velocity_settings = KratosMultiphysics.Parameters("""
            {
                "mesh_id"              : """ + str(mesh_id) + """,
                "model_part_name"      : \"""" + model_part_name + """\",
                "variable_name"        : "VELOCITY",
                "interval"             : """ + interval + """,
                "value"                : [0.0, 0.0, 0.0],
                "constrained"          : [true,true,true],
                "local_axes"           : {}
            }
            """
            )

        adjoint_velocity_settings = KratosMultiphysics.Parameters("""
            {
                "mesh_id"              : """ + str(mesh_id) + """,
                "model_part_name"      : \"""" + model_part_name + """\",
                "variable_name"        : "VELOCITY_ADJ",
                "interval"             : """ + interval + """,
                "value"                : [0.0, 0.0, 0.0],
                "constrained"          : [true,true,true],
                "local_axes"           : {}
            }
            """
            )
        
        self.no_slip_model_part = Model[settings["model_part_name"].GetString()]
    
        self.velocity_vector_process = AssignVectorVariableProcess(Model, velocity_settings)
        self.adjoint_velocity_vector_process = AssignVectorVariableProcess(Model, adjoint_velocity_settings)

    def ExecuteInitializeSolutionStep(self):
        # Call the base process ExecuteInitializeSolutionStep()
        if self.IsPhysicsStage():
            self.velocity_vector_process.ExecuteInitializeSolutionStep()
        elif self.IsAdjointStage():
            self.adjoint_velocity_vector_process.ExecuteInitializeSolutionStep()
        else:
            KratosMultiphysics.Logger.PrintError("'ExecuteInitializeSolutionStep' for ApplyFluidTopologyOptimizationNoSlipProcess called during a non valid stage.")

    def ExecuteFinalizeSolutionStep(self):
        # Call the base process ExecuteFinalizeSolutionStep()
        if self.IsPhysicsStage():
            self.velocity_vector_process.ExecuteFinalizeSolutionStep()
        elif self.IsAdjointStage():
            self.adjoint_velocity_vector_process.ExecuteFinalizeSolutionStep()
        else:
            KratosMultiphysics.Logger.PrintError("'ExecuteFinalizeSolutionStep' for ApplyFluidTopologyOptimizationNoSlipProcess called during a non valid stage.")

    def IsPhysicsStage(self):
        top_opt_stage = self.no_slip_model_part.ProcessInfo.GetValue(KratosCFD.FLUID_TOP_OPT_PROBLEM_STAGE)
        if top_opt_stage == 1:
            return True
        else:
            return False
        
    def IsAdjointStage(self):
        top_opt_stage = self.no_slip_model_part.ProcessInfo.GetValue(KratosCFD.FLUID_TOP_OPT_PROBLEM_STAGE)
        if top_opt_stage == 2:
            return True
        else:
            return False

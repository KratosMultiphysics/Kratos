import KratosMultiphysics
import KratosMultiphysics.ConvectionDiffusionApplication as KratosCD
from KratosMultiphysics.assign_scalar_variable_process import AssignScalarVariableProcess

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplyTransportTopologyOptimizationImposeTransportScalarFluxProcess(Model, settings["Parameters"])


class ApplyTransportTopologyOptimizationImposeTransportScalarFluxProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings):

        KratosMultiphysics.Process.__init__(self)
        
        default_settings = KratosMultiphysics.Parameters("""
            {
                "mesh_id"         : 0,
                "model_part_name" : "please_specify_model_part_name",
                "interval"        : [0.0,1e30],
                "transport_scalar_flux_value"         : {},
                "adjoint_transport_scalar_flux_value" : {}
            }
            """)

        # Check input data
        if (settings["model_part_name"].GetString() == ""):
            raise Exception("Empty inlet model part name string. Set a valid model part name.")
        
        # Trick: allows "value" to be a double, a string or a table value (otherwise the ValidateAndAssignDefaults might fail)
        if(settings.Has("transport_scalar_flux_value")):
            if(settings["transport_scalar_flux_value"].IsString()):
                default_settings["transport_scalar_flux_value"].SetString("0.0")
            elif settings["transport_scalar_flux_value"].IsNumber():
                default_settings["transport_scalar_flux_value"].SetDouble(0.0)
        else:
            err_msg = "Provided settings have no 'transport_scalar_flux_value'. This needs to be provided."
            raise Exception(err_msg)
        if(settings.Has("adjoint_transport_scalar_flux_value")):
            if(settings["adjoint_transport_scalar_flux_value"].IsString()):
                default_settings["adjoint_transport_scalar_flux_value"].SetString("0.0")
            elif settings["adjoint_transport_scalar_flux_value"].IsNumber():
                default_settings["adjoint_transport_scalar_flux_value"].SetDouble(0.0)
        else:
            err_msg = "Provided settings have no 'adjoint_transport_scalar_flux_value'. This needs to be provided."
            raise Exception(err_msg)

        settings.ValidateAndAssignDefaults(default_settings)

        mesh_id = settings["mesh_id"].GetInt()
        model_part_name = settings["model_part_name"].GetString()
        interval = settings["interval"].PrettyPrintJsonString()
        if (settings["transport_scalar_flux_value"].IsString()):
            transport_scalar_flux_value = settings["transport_scalar_flux_value"].GetString()
        else:
            transport_scalar_flux_value = str(settings["transport_scalar_flux_value"].GetDouble())
        if (settings["adjoint_transport_scalar_flux_value"].IsString()):
            adjoint_transport_scalar_flux_value = settings["adjoint_transport_scalar_flux_value"].GetString()
        else:
            adjoint_transport_scalar_flux_value = str(settings["adjoint_transport_scalar_flux_value"].GetDouble())

        transport_scalar_flux_settings = KratosMultiphysics.Parameters("""
            {
                "mesh_id"              : """ + str(mesh_id) + """,
                "model_part_name"      : \"""" + model_part_name + """\",
                "variable_name"        : "FACE_HEAT_FLUX",
                "interval"             : """ + interval + """,
                "value"                : """ + transport_scalar_flux_value + """,
                "constrained"          : false,
                "local_axes"           : {}
            }
            """
            )

        adjoint_transport_scalar_flux_settings = KratosMultiphysics.Parameters("""
            {
                "mesh_id"              : """ + str(mesh_id) + """,
                "model_part_name"      : \"""" + model_part_name + """\",
                "variable_name"        : "FACE_HEAT_FLUX_ADJ",
                "interval"             : """ + interval + """,
                "value"                : """ + adjoint_transport_scalar_flux_value + """,
                "constrained"          : false,
                "local_axes"           : {}
            }
            """
            )
        
        self.impose_transport_scalar_flux_model_part = Model[settings["model_part_name"].GetString()]
    
        self.transport_scalar_flux_process = AssignScalarVariableProcess(Model, transport_scalar_flux_settings)
        self.adjoint_transport_scalar_flux_process = AssignScalarVariableProcess(Model, adjoint_transport_scalar_flux_settings)


    def ExecuteInitializeSolutionStep(self):
        # Call the base process ExecuteInitializeSolutionStep()
        if self.IsPhysicsStage():
            self.transport_scalar_flux_process.ExecuteInitializeSolutionStep()
        elif self.IsAdjointStage():
            self.adjoint_transport_scalar_flux_process.ExecuteInitializeSolutionStep()
        else:
            KratosMultiphysics.Logger.PrintError("'ExecuteInitializeSolutionStep' for ApplyTransportTopologyOptimizationImposeTransportScalarProcess called during a non valid stage.")

    def ExecuteFinalizeSolutionStep(self):
        # Call the base process ExecuteFinalizeSolutionStep()
        if self.IsPhysicsStage():
            self.transport_scalar_flux_process.ExecuteFinalizeSolutionStep()
        elif self.IsAdjointStage():
            self.adjoint_transport_scalar_flux_process.ExecuteFinalizeSolutionStep()
        else:
            KratosMultiphysics.Logger.PrintError("'ExecuteFinalizeSolutionStep' for ApplyTransportTopologyOptimizationImposeTransportScalarProcess called during a non valid stage.")

    def IsPhysicsStage(self):
        top_opt_stage = self.impose_transport_scalar_flux_model_part.ProcessInfo.GetValue(KratosCD.TRANSPORT_TOP_OPT_PROBLEM_STAGE)
        if top_opt_stage == 1:
            return True
        else:
            return False
        
    def IsAdjointStage(self):
        top_opt_stage = self.impose_transport_scalar_flux_model_part.ProcessInfo.GetValue(KratosCD.TRANSPORT_TOP_OPT_PROBLEM_STAGE)
        if top_opt_stage == 2:
            return True
        else:
            return False

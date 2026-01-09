import KratosMultiphysics
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD
import KratosMultiphysics.FluidDynamicsBiomedicalApplication as KratosBio

from KratosMultiphysics import assign_vector_by_direction_process
from KratosMultiphysics.FluidDynamicsBiomedicalApplication import apply_parabolic_inlet_process

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplyFluidTopologyOptimization3DParabolicInletProcess(Model, settings["Parameters"])

class ApplyFluidTopologyOptimization3DParabolicInletProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings):
        KratosMultiphysics.Process.__init__(self)

        default_settings = KratosMultiphysics.Parameters("""
            {
                "wall_model_part_name": "",
                "inlet_model_part_name": "",
                "velocity_value"         : {},
                "value_is_average" : false,
                "value_is_flow_rate" : false,
                "adjoint_velocity_value" : "0.0",
                "adjoint_direction": "automatic_inwards_normal",
                "interval"        : [0.0,"End"]
            }
            """)

        # Check user-provided data
        if not settings["inlet_model_part_name"].GetString():
            raise ValueError("'inlet_model_part' not provided.")

        if not settings["wall_model_part_name"].GetString():
            raise ValueError("'wall_model_part' needs to be provided for 'parabolic' inlet distribution.")

        if settings.Has("velocity_value"):
            if settings["velocity_value"].IsString():
                default_settings["velocity_value"].SetString("0.0")
            elif settings["velocity_value"].IsDouble():
                default_settings["velocity_value"].SetDouble(0.0)
        else:
            raise Exception("'velocity_value' not found. It needs to be user-provided.")
        
        if not settings.Has("adjoint_velocity_value"):
            raise Exception("'adjoint_velocity_value' not found. It needs to be user-provided.")
        
        settings.ValidateAndAssignDefaults(default_settings)

        wall_model_part_name = settings["wall_model_part_name"].GetString()
        inlet_model_part_name = settings["inlet_model_part_name"].GetString()
        # Set the maximum velocity value string
        if settings["velocity_value"].IsString():
            velocity_value = "\"" + settings["velocity_value"].GetString() + "\""
        else:
            velocity_value = settings["velocity_value"].PrettyPrintJsonString()
        value_is_average = settings["value_is_average"].GetBool()
        if value_is_average:
            value_is_average_str = "true"
        else:
            value_is_average_str = "false"
        value_is_flow_rate = settings["value_is_flow_rate"].GetBool()
        if value_is_flow_rate:
            value_is_flow_rate_str = "true"
        else:
            value_is_flow_rate_str = "false"
        adjoint_velocity_value = settings["adjoint_velocity_value"].GetString()
        adjoint_direction = settings["adjoint_direction"].GetString()
        interval = settings["interval"].PrettyPrintJsonString()

        velocity_settings = KratosMultiphysics.Parameters("""
            {
                "wall_model_part_name" : \"""" + wall_model_part_name + """\",
                "inlet_model_part_name" : \"""" + inlet_model_part_name + """\",
                "value"         : """ + velocity_value + """,
                "value_is_average"  : """ + value_is_average_str + """,
                "value_is_flow_rate" : """ + value_is_flow_rate_str + """,
                "interval"        : """ + interval + """
            }
            """)

        adjoint_velocity_settings = KratosMultiphysics.Parameters("""
            {
                "model_part_name" : \"""" + inlet_model_part_name + """\",
                "variable_name"   : "VELOCITY_ADJ",
                "modulus"         : \"""" + adjoint_velocity_value + """\",
                "constrained"     : true,
                "direction"       : \"""" + adjoint_direction + """\",
                "interval"        : """ + interval + """
            }
            """)
            
        # Set the INLET flag in the inlet model part nodes and conditions
        self.inlet_model_part = Model[settings["inlet_model_part_name"].GetString()]
        for node in self.inlet_model_part.Nodes:
            node.Set(KratosMultiphysics.INLET, True)
        for condition in self.inlet_model_part.Conditions:
            condition.Set(KratosMultiphysics.INLET, True)

        # Construct the base process AssignVectorByDirectionProcess
        self.aux_process         = apply_parabolic_inlet_process.ApplyParabolicInletProcess(Model, velocity_settings)
        self.aux_process_adjoint = assign_vector_by_direction_process.AssignVectorByDirectionProcess(Model, adjoint_velocity_settings)

    def ExecuteBeforeSolutionLoop(self):
        self.aux_process.ExecuteBeforeSolutionLoop()

    def ExecuteInitializeSolutionStep(self):
        # Call the base process ExecuteInitializeSolutionStep()
        if self.IsPhysicsStage():
            self.aux_process.ExecuteInitializeSolutionStep()
        elif self.IsAdjointStage():
            self.aux_process_adjoint.ExecuteInitializeSolutionStep()
        else:
            KratosMultiphysics.Logger.PrintError("'ExecuteInitializeSolutionStep' for ApplyFluidTopologyOptimization3DParabolicInletProcess called during a non valid stage.")

    def ExecuteFinalizeSolutionStep(self):
        # Call the base process ExecuteFinalizeSolutionStep()
        if self.IsPhysicsStage():
            self.aux_process.ExecuteFinalizeSolutionStep()
        elif self.IsAdjointStage():
            self.aux_process_adjoint.ExecuteFinalizeSolutionStep()
        else:
            KratosMultiphysics.Logger.PrintError("'ExecuteFinalizeSolutionStep' for ApplyFluidTopologyOptimization3DParabolicInletProcess called during a non valid stage.")

    def IsPhysicsStage(self):
        top_opt_stage = self.inlet_model_part.ProcessInfo.GetValue(KratosCFD.FLUID_TOP_OPT_PROBLEM_STAGE)
        if top_opt_stage == 1:
            return True
        else:
            return False
        
    def IsAdjointStage(self):
        top_opt_stage = self.inlet_model_part.ProcessInfo.GetValue(KratosCFD.FLUID_TOP_OPT_PROBLEM_STAGE)
        if top_opt_stage == 2:
            return True
        else:
            return False

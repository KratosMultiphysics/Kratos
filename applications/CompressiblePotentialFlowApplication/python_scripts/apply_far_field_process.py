import KratosMultiphysics
import KratosMultiphysics.CompressiblePotentialFlowApplication as CPFApp
import math

def Factory(settings, Model):
    if( not isinstance(settings,KratosMultiphysics.Parameters) ):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplyFarFieldProcess(Model, settings["Parameters"])

## All the processes python should be derived from "Process"
class ApplyFarFieldProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings ):
        KratosMultiphysics.Process.__init__(self)

        default_parameters = KratosMultiphysics.Parameters( """
            {
                "model_part_name":"",
                "angle_of_attack": 0.0,
                "mach_infinity": 0.02941176471,
                "free_stream_density"  : 1.0,
                "speed_of_sound": 340,
                "heat_capacity_ratio": 1.4,
                "inlet_potential": 1.0,
                "free_stream_velocity": [1.0,0.0,0],
                "initialize_flow_field": true
            }  """ );


        settings.ValidateAndAssignDefaults(default_parameters);

        self.model_part = Model[settings["model_part_name"].GetString()]
        self.fluid_model_part = self.model_part.GetRootModelPart()

        self.angle_of_attack = settings["angle_of_attack"].GetDouble()
        self.free_stream_mach = settings["mach_infinity"].GetDouble()
        self.density_inf = settings["free_stream_density"].GetDouble()
        self.free_stream_speed_of_sound = settings["speed_of_sound"].GetDouble()
        self.heat_capacity_ratio = settings["heat_capacity_ratio"].GetDouble()
        self.inlet_potential = settings["inlet_potential"].GetDouble()
        self.initialize_flow_field = settings["initialize_flow_field"].GetBool()

        # Computing free stream velocity
        self.u_inf = self.free_stream_mach * self.free_stream_speed_of_sound
        self.free_stream_velocity = KratosMultiphysics.Vector(3)
        self.free_stream_velocity[0] = round(self.u_inf*math.cos(self.angle_of_attack),8)
        self.free_stream_velocity[1] = round(self.u_inf*math.sin(self.angle_of_attack),8)
        self.free_stream_velocity[2] = 0.0

        self.fluid_model_part.ProcessInfo.SetValue(CPFApp.FREE_STREAM_MACH,self.free_stream_mach)
        self.fluid_model_part.ProcessInfo.SetValue(CPFApp.FREE_STREAM_VELOCITY,self.free_stream_velocity)
        self.fluid_model_part.ProcessInfo.SetValue(CPFApp.FREE_STREAM_DENSITY,self.density_inf)
        self.fluid_model_part.ProcessInfo.SetValue(KratosMultiphysics.SOUND_VELOCITY,self.free_stream_speed_of_sound)
        self.fluid_model_part.ProcessInfo.SetValue(CPFApp.HEAT_CAPACITY_RATIO,self.heat_capacity_ratio)

    def Execute(self):
        #KratosMultiphysics.VariableUtils().SetVectorVar(CPFApp.FREE_STREAM_VELOCITY, self.velocity_infinity, self.model_part.Conditions)
        for cond in self.model_part.Conditions:
            cond.SetValue(CPFApp.FREE_STREAM_VELOCITY, self.free_stream_velocity)

        #select the first node
        for node in self.model_part.Nodes:
            node1 = node
            break

        #find the node with the minimal x
        x0 = node1.X
        y0 = node1.X
        z0 = node1.X

        pos = 1e30
        for node in self.model_part.Nodes:
            dx = node.X - x0
            dy = node.Y - y0
            dz = node.Z - z0

            tmp = dx*self.free_stream_velocity[0] + dy*self.free_stream_velocity[1] + dz*self.free_stream_velocity[2]

            if(tmp < pos):
                pos = tmp
                self.reference_inlet_node = node

        for node in self.model_part.Nodes:
            dx = node.X - x0
            dy = node.Y - y0
            dz = node.Z - z0

            tmp = dx*self.free_stream_velocity[0] + dy*self.free_stream_velocity[1] + dz*self.free_stream_velocity[2]

            if(tmp < pos+1e-9):
                node.Fix(CPFApp.VELOCITY_POTENTIAL)

                node.SetSolutionStepValue(CPFApp.VELOCITY_POTENTIAL,0,self.inlet_potential)
                if self.model_part.HasNodalSolutionStepVariable(CPFApp.ADJOINT_VELOCITY_POTENTIAL):
                    node.Fix(CPFApp.ADJOINT_VELOCITY_POTENTIAL)
                    node.SetSolutionStepValue(CPFApp.ADJOINT_VELOCITY_POTENTIAL,0,0.0)

        if(self.initialize_flow_field):
            for node in self.fluid_model_part.Nodes:
                # Computing distance to reference
                dx = node.X - self.reference_inlet_node.X
                dy = node.Y - self.reference_inlet_node.Y
                dz = node.Z - self.reference_inlet_node.Z

                initial_potential = dx*self.free_stream_velocity[0] + dy*self.free_stream_velocity[1] + dz*self.free_stream_velocity[2]
                node.SetSolutionStepValue(CPFApp.VELOCITY_POTENTIAL,0,initial_potential + self.inlet_potential)
                node.SetSolutionStepValue(CPFApp.AUXILIARY_VELOCITY_POTENTIAL,0,initial_potential + self.inlet_potential)

    def ExecuteInitializeSolutionStep(self):
        self.Execute()


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
                "model_part_name":"PLEASE_CHOOSE_MODEL_PART_NAME",
                "angle_of_attack": 0.0,
                "mach_infinity": 0.02941176471,
                "density_infinity"  : 1.0,
                "speed_of_sound": 340,
                "heat_capacity_ratio": 1.4,
                "mesh_id": 0,
                "inlet_phi": 1.0,
                "velocity_infinity": [1.0,0.0,0]
            }  """ );


        settings.ValidateAndAssignDefaults(default_parameters);

        self.model_part = Model[settings["model_part_name"].GetString()]
        self.fluid_model_part = self.model_part.GetRootModelPart()

        self.angle_of_attack = settings["angle_of_attack"].GetDouble()
        self.mach_inf = settings["mach_infinity"].GetDouble()
        self.density_inf = settings["density_infinity"].GetDouble()
        self.seed_of_sound_inf = settings["speed_of_sound"].GetDouble()
        self.heat_capacity_ratio = settings["heat_capacity_ratio"].GetDouble()
        self.inlet_phi = settings["inlet_phi"].GetDouble()

        # Computing free stream velocity
        self.u_inf = self.mach_inf * self.seed_of_sound_inf
        self.velocity_infinity = KratosMultiphysics.Vector(3)
        self.velocity_infinity[0] = round(self.u_inf*math.cos(self.angle_of_attack),8)
        self.velocity_infinity[1] = round(self.u_inf*math.sin(self.angle_of_attack),8)
        self.velocity_infinity[2] = 0.0

        self.fluid_model_part.ProcessInfo.SetValue(CPFApp.MACH_INFINITY,self.mach_inf)
        self.fluid_model_part.ProcessInfo.SetValue(CPFApp.VELOCITY_INFINITY,self.velocity_infinity)
        self.fluid_model_part.ProcessInfo.SetValue(CPFApp.DENSITY_INFINITY,self.density_inf)
        self.fluid_model_part.ProcessInfo.SetValue(KratosMultiphysics.SOUND_VELOCITY,self.seed_of_sound_inf)
        self.fluid_model_part.ProcessInfo.SetValue(CPFApp.HEAT_CAPACITY_RATIO,self.heat_capacity_ratio)

    def Execute(self):
        #KratosMultiphysics.VariableUtils().SetVectorVar(CPFApp.VELOCITY_INFINITY, self.velocity_infinity, self.model_part.Conditions)
        for cond in self.model_part.Conditions:
            cond.SetValue(CPFApp.VELOCITY_INFINITY, self.velocity_infinity)

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

            tmp = dx*self.velocity_infinity[0] + dy*self.velocity_infinity[1] + dz*self.velocity_infinity[2]

            if(tmp < pos):
                pos = tmp

        for node in self.model_part.Nodes:
            dx = node.X - x0
            dy = node.Y - y0
            dz = node.Z - z0

            tmp = dx*self.velocity_infinity[0] + dy*self.velocity_infinity[1] + dz*self.velocity_infinity[2]

            if(tmp < pos+1e-9):
                node.Fix(CPFApp.VELOCITY_POTENTIAL)
                node.SetSolutionStepValue(CPFApp.VELOCITY_POTENTIAL,0,self.inlet_phi)
                if self.model_part.HasNodalSolutionStepVariable(CPFApp.ADJOINT_VELOCITY_POTENTIAL):
                    node.Fix(CPFApp.ADJOINT_VELOCITY_POTENTIAL)
                    node.SetSolutionStepValue(CPFApp.ADJOINT_VELOCITY_POTENTIAL,0,0.0)

    def ExecuteInitializeSolutionStep(self):
        self.Execute()


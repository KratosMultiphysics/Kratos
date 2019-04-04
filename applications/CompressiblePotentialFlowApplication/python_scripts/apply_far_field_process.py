import KratosMultiphysics
import KratosMultiphysics.CompressiblePotentialFlowApplication as CompressiblePotentialFlowApplication
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
                "inlet_phi": 1.0,
                "velocity_infinity": [3.4,0.0,0],
                "density_infinity"  : 1.0,
                "mach_infinity": 0.01,
                "gamma": 1.4,
                "pressure_infinity": 101325
            }  """)

        settings.ValidateAndAssignDefaults(default_parameters)

        self.model_part = Model[settings["model_part_name"].GetString()]
        self.fluid_model_part = self.model_part.GetRootModelPart()

        self.inlet_phi = settings["inlet_phi"].GetDouble()
        self.velocity_infinity = KratosMultiphysics.Vector(3)
        self.velocity_infinity[0] = settings["velocity_infinity"][0].GetDouble(
        )
        self.velocity_infinity[1] = settings["velocity_infinity"][1].GetDouble(
        )
        self.velocity_infinity[2] = settings["velocity_infinity"][2].GetDouble(
        )
        self.density_infinity = settings["density_infinity"].GetDouble()
        self.mach_infinity = settings["mach_infinity"].GetDouble()
        self.gamma = settings["gamma"].GetDouble()
        self.pressure_infinity = settings["pressure_infinity"].GetDouble()

        self.u_infinity = math.sqrt(
            self.velocity_infinity[0]**2 + self.velocity_infinity[1]**2 + self.velocity_infinity[2]**2)
        self.a_infinity = self.u_infinity / self.mach_infinity

        # For the model part
        self.model_part.ProcessInfo.SetValue(
            CompressiblePotentialFlowApplication.VELOCITY_INFINITY, self.velocity_infinity)

        # For the conditions
        self.fluid_model_part.GetProperties()[0].SetValue(
            CompressiblePotentialFlowApplication.DENSITY_INFINITY, self.density_infinity)

        # For the elements
        self.fluid_model_part.GetProperties()[1].SetValue(
            CompressiblePotentialFlowApplication.VELOCITY_INFINITY, self.velocity_infinity)
        self.fluid_model_part.GetProperties()[1].SetValue(
            CompressiblePotentialFlowApplication.DENSITY_INFINITY, self.density_infinity)
        self.fluid_model_part.GetProperties()[1].SetValue(
            CompressiblePotentialFlowApplication.MACH_INFINITY, self.mach_infinity)
        self.fluid_model_part.GetProperties()[1].SetValue(
            CompressiblePotentialFlowApplication.GAMMA, self.gamma)
        self.fluid_model_part.GetProperties()[1].SetValue(
            KratosMultiphysics.SOUND_VELOCITY, self.a_infinity)
        self.fluid_model_part.GetProperties()[1].SetValue(
            CompressiblePotentialFlowApplication.PRESSURE_INFINITY, self.pressure_infinity)

    def Execute(self):
        #KratosMultiphysics.VariableUtils().SetVectorVar(CompressiblePotentialFlowApplication.VELOCITY_INFINITY, self.velocity_infinity, self.model_part.Conditions)
        for cond in self.model_part.Conditions:
            cond.SetValue(CompressiblePotentialFlowApplication.VELOCITY_INFINITY, self.velocity_infinity)

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
                node.Fix(CompressiblePotentialFlowApplication.VELOCITY_POTENTIAL)
                node.SetSolutionStepValue(CompressiblePotentialFlowApplication.VELOCITY_POTENTIAL,0,self.inlet_phi)
                if self.model_part.HasNodalSolutionStepVariable(CompressiblePotentialFlowApplication.ADJOINT_VELOCITY_POTENTIAL):
                    node.Fix(CompressiblePotentialFlowApplication.ADJOINT_VELOCITY_POTENTIAL)
                    node.SetSolutionStepValue(CompressiblePotentialFlowApplication.ADJOINT_VELOCITY_POTENTIAL,0,0.0)
        
    def ExecuteInitializeSolutionStep(self):
        self.Execute()
        
        

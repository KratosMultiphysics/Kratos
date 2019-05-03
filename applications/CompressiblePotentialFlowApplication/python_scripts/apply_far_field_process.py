import KratosMultiphysics
import KratosMultiphysics.CompressiblePotentialFlowApplication as CompressiblePotentialFlowApplication
#from CompressiblePotentialFlowApplication import*

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
                "mesh_id": 0,
                "inlet_phi": 1.0,
                "free_stream_velocity": [1.0,0.0,0]
            }  """ );


        settings.ValidateAndAssignDefaults(default_parameters);

        self.model_part = Model[settings["model_part_name"].GetString()]
        self.free_stream_velocity = settings["free_stream_velocity"].GetVector()
        self.inlet_phi = settings["inlet_phi"].GetDouble()
        self.model_part.ProcessInfo.SetValue(CompressiblePotentialFlowApplication.FREE_STREAM_VELOCITY,self.free_stream_velocity)

    def Execute(self):
        #KratosMultiphysics.VariableUtils().SetVectorVar(CompressiblePotentialFlowApplication.FREE_STREAM_VELOCITY, self.free_stream_velocity, self.model_part.Conditions)
        for cond in self.model_part.Conditions:
            cond.SetValue(CompressiblePotentialFlowApplication.FREE_STREAM_VELOCITY, self.free_stream_velocity)

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

        for node in self.model_part.Nodes:
            dx = node.X - x0
            dy = node.Y - y0
            dz = node.Z - z0

            tmp = dx*self.free_stream_velocity[0] + dy*self.free_stream_velocity[1] + dz*self.free_stream_velocity[2]

            if(tmp < pos+1e-9):
                node.Fix(CompressiblePotentialFlowApplication.VELOCITY_POTENTIAL)
                node.SetSolutionStepValue(CompressiblePotentialFlowApplication.VELOCITY_POTENTIAL,0,self.inlet_phi)
                if self.model_part.HasNodalSolutionStepVariable(CompressiblePotentialFlowApplication.ADJOINT_VELOCITY_POTENTIAL):
                    node.Fix(CompressiblePotentialFlowApplication.ADJOINT_VELOCITY_POTENTIAL)
                    node.SetSolutionStepValue(CompressiblePotentialFlowApplication.ADJOINT_VELOCITY_POTENTIAL,0,0.0)

    def ExecuteInitializeSolutionStep(self):
        self.Execute()


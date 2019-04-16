import KratosMultiphysics
import KratosMultiphysics.CompressiblePotentialFlowApplication as CPFApp

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
                "velocity_infinity": [1.0,0.0,0],
                "initialize_flow_field": true
            }  """ );


        settings.ValidateAndAssignDefaults(default_parameters);

        self.model_part = Model[settings["model_part_name"].GetString()]
        self.fluid_model_part = self.model_part.GetRootModelPart()
        self.velocity_infinity = settings["velocity_infinity"].GetVector()
        self.velocity_infinity[0] = settings["velocity_infinity"][0].GetDouble()
        self.velocity_infinity[1] = settings["velocity_infinity"][1].GetDouble()
        self.velocity_infinity[2] = settings["velocity_infinity"][2].GetDouble()
        self.inlet_potential = settings["inlet_phi"].GetDouble()
        self.model_part.ProcessInfo.SetValue(CPFApp.VELOCITY_INFINITY,self.velocity_infinity)
        self.initialize_flow_field = settings["initialize_flow_field"].GetBool()

    def Execute(self):
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
                self.reference_inlet_node = node

        for node in self.model_part.Nodes:
            dx = node.X - x0
            dy = node.Y - y0
            dz = node.Z - z0

            tmp = dx*self.velocity_infinity[0] + dy*self.velocity_infinity[1] + dz*self.velocity_infinity[2]

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

                initial_potential = dx*self.velocity_infinity[0] + dy*self.velocity_infinity[1] + dz*self.velocity_infinity[2]
                node.SetSolutionStepValue(CPFApp.VELOCITY_POTENTIAL,0,initial_potential + self.inlet_potential)
                node.SetSolutionStepValue(CPFApp.AUXILIARY_VELOCITY_POTENTIAL,0,initial_potential + self.inlet_potential)

    def ExecuteInitializeSolutionStep(self):
        self.Execute()


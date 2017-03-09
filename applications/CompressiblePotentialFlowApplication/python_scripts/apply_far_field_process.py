import KratosMultiphysics

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplyFarFieldProcess(Model, settings["Parameters"])

##all the processes python processes should be derived from "python_process"
class ApplyFarFieldProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings ):
        KratosMultiphysics.Process.__init__(self)
        
        default_parameters = KratosMultiphysics.Parameters( """
            {
                "model_part_name":"PLEASE_CHOOSE_MODEL_PART_NAME",
                "mesh_id": 0,
                "inlet_phi": 1.0,
                "velocity_infinity": [1.0,0.0,0],
                "density_infinity": 1.2,
                "far_field_model_part" : "MainModelPart"
            }  """ );
        
            
        settings.ValidateAndAssignDefaults(default_parameters);
        
        self.model_part = Model[settings["model_part_name"].GetString()]
        self.velocity_infinity = [0,0,0]#array('d', [1.0, 2.0, 3.14])#np.array([0,0,0])#np.zeros(3)#vector(3)
        self.velocity_infinity[0] = settings["velocity_infinity"][0].GetDouble()
        self.velocity_infinity[1] = settings["velocity_infinity"][1].GetDouble()
        self.velocity_infinity[2] = settings["velocity_infinity"][2].GetDouble()
        self.density_infinity = settings["density_infinity"].GetDouble()
        self.inlet_phi = settings["inlet_phi"].GetDouble()
        
        self.model_part = Model[settings["far_field_model_part"].GetString()]
        
    def Execute(self):
        for cond in self.model_part.Conditions:
            cond.SetValue(KratosMultiphysics.VELOCITY, self.velocity_infinity)

        #find the node with the minimal x
        node1 = self.model_part.Nodes[1]
        x0 = node1.X
        y0 = node1.X
        z0 = node1.X
        
        ref_id = 1
        pos = 1e30
        for node in self.model_part.Nodes:
            dx = node.X - x0
            dy = node.Y - y0
            dz = node.Z - z0
            
            tmp = dx*self.velocity_infinity[0] + dy*self.velocity_infinity[1] + dz*self.velocity_infinity[2]
            if(tmp < pos):
              ref_id = node.Id
              
        self.model_part.Nodes[ref_id].Fix(KratosMultiphysics.POSITIVE_FACE_PRESSURE)
        self.model_part.Nodes[ref_id].SetSolutionStepValue(KratosMultiphysics.POSITIVE_FACE_PRESSURE,0,self.inlet_phi)
        
    def ExecuteInitializeSolutionStep(self):
        self.Execute()
        
        
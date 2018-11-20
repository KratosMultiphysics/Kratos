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
                "velocity_infinity": [1.0,0.0,0]
            }  """ );
        
            
        settings.ValidateAndAssignDefaults(default_parameters);
        self.model_part = Model.GetModelPart(settings["model_part_name"].GetString())
        self.velocity_infinity = KratosMultiphysics.Vector(3)#array('d', [1.0, 2.0, 3.14])#np.array([0,0,0])#np.zeros(3)#vector(3)
        self.velocity_infinity[0] = settings["velocity_infinity"][0].GetDouble()
        self.velocity_infinity[1] = settings["velocity_infinity"][1].GetDouble()
        self.velocity_infinity[2] = settings["velocity_infinity"][2].GetDouble()
        #self.density_infinity = settings["density_infinity"].GetDouble() #TODO: must read this from the properties
        self.inlet_phi = settings["inlet_phi"].GetDouble()
        
        self.model_part.ProcessInfo.SetValue(CompressiblePotentialFlowApplication.VELOCITY_INFINITY,self.velocity_infinity)
        
        
        
    def Execute(self):
        #KratosMultiphysics.VariableUtils().SetVectorVar(CompressiblePotentialFlowApplication.VELOCITY_INFINITY, self.velocity_infinity, self.model_part.Conditions)
        for cond in self.model_part.Conditions:
            cond.SetValue(CompressiblePotentialFlowApplication.VELOCITY_INFINITY, self.velocity_infinity)
            npos=0
            nneg=0
            # for node in cond.GetNodes():
            #     distance=node.GetSolutionStepValue(CompressiblePotentialFlowApplication.WAKE_DISTANCE)
            #     if distance<0:
            #         npos += 1
            #     else:
            #         nneg += 1
            # if not (npos*nneg==0):
            #     cond.Set(KratosMultiphysics.STRUCTURE,True)
            # else:
            #     cond.Set(KratosMultiphysics.STRUCTURE,False)
                

        #select the first node
        for node in self.model_part.Nodes:
            node1 = node
            break
        
        #find the node with the minimal x
        x0 = node1.X
        y0 = node1.Y
        z0 = node1.Z

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
                    node.Set(KratosMultiphysics.INLET)
                    node.Fix(CompressiblePotentialFlowApplication.POSITIVE_POTENTIAL)
                    node.SetSolutionStepValue(CompressiblePotentialFlowApplication.POSITIVE_POTENTIAL,0,self.inlet_phi)
        
    def ExecuteInitializeSolutionStep(self):
        self.Execute()
        
        
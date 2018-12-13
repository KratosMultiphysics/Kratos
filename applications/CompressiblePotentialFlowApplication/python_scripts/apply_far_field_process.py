import KratosMultiphysics
import KratosMultiphysics.CompressiblePotentialFlowApplication as CompressiblePotentialFlowApplication
import numpy as np
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
        self.main_model_part = self.model_part.GetRootModelPart()
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
            for node in cond.GetNodes():
                distance=node.GetSolutionStepValue(CompressiblePotentialFlowApplication.WAKE_DISTANCE)
                if distance>0:
                    npos += 1
                elif distance<0:
                    nneg += 1
            if (npos>0 and nneg>0):
                cond.Set(KratosMultiphysics.STRUCTURE,True)
            else:
                cond.Set(KratosMultiphysics.STRUCTURE,False)
                
        
        # #select the first node
        # for node in self.model_part.Nodes:
        #     node1 = node
        #     break
        
        # #find the node with the minimal x
        # x0 = node1.X
        # y0 = node1.Y
        # z0 = node1.Z

        # pos = 1e30
        # for node in self.model_part.Nodes:
        #         dx = node.X - x0
        #         dy = node.Y - y0
        #         dz = node.Z - z0

        #         tmp = dx*self.velocity_infinity[0] + dy*self.velocity_infinity[1] + dz*self.velocity_infinity[2]

        #         if(tmp < pos):
        #             pos = tmp
      
        # for node in self.model_part.Nodes:
        #         dx = node.X - x0
        #         dy = node.Y - y0
        #         dz = node.Z - z0
            
        #         tmp = dx*self.velocity_infinity[0] + dy*self.velocity_infinity[1] + dz*self.velocity_infinity[2]
                
        #         if(tmp < pos+1e-9):
        #             node.Set(KratosMultiphysics.INLET)
        #             node.Fix(CompressiblePotentialFlowApplication.POSITIVE_POTENTIAL)
        #             node.SetSolutionStepValue(CompressiblePotentialFlowApplication.POSITIVE_POTENTIAL,0,self.inlet_phi)
        for cond in self.model_part.Conditions:
            normal=cond.GetValue(KratosMultiphysics.NORMAL)
            v_inf=cond.GetValue(CompressiblePotentialFlowApplication.VELOCITY_INFINITY)
           
            value = np.dot(normal,v_inf)

            if value<0:
                for node in cond.GetNodes():
                    inlet_phi=node.X*self.velocity_infinity[0] + node.Y*self.velocity_infinity[1] + node.Z*self.velocity_infinity[2]
                    node.Fix(CompressiblePotentialFlowApplication.POSITIVE_POTENTIAL)
                    node.SetSolutionStepValue(CompressiblePotentialFlowApplication.POSITIVE_POTENTIAL,0,inlet_phi)



        
        #             inlet_phi_min=node.X*self.velocity_infinity[0] + node.Y*self.velocity_infinity[1] + node.Z*self.velocity_infinity[2]
        for node in self.main_model_part.Nodes:
            initial_phi=node.X*self.velocity_infinity[0] + node.Y*self.velocity_infinity[1] + node.Z*self.velocity_infinity[2]
            node.SetSolutionStepValue(CompressiblePotentialFlowApplication.POSITIVE_POTENTIAL,0,initial_phi)
            node.SetSolutionStepValue(CompressiblePotentialFlowApplication.NEGATIVE_POTENTIAL,0,initial_phi)
        
        # for cond in self.model_part.Conditions:
        #     n_inlet=0
        #     n_outlet=0
        #     n_wall=0
        #     for node in cond.GetNodes():
        #         if node.X == -50.0:
        #             n_inlet += 1
        #         elif node.X == 50:
        #             n_outlet += 1
        #         else:
        #             n_wall +=1
        #     if n_outlet==0 and n_wall==0:
               
        #         cond.SetValue(CompressiblePotentialFlowApplication.VELOCITY_INFINITY, self.velocity_infinity)
        #         for node in cond.GetNodes():
        #             potential=node.X*self.velocity_infinity[0] + node.Y*self.velocity_infinity[1] + node.Z*self.velocity_infinity[2]
        #             node.Fix(CompressiblePotentialFlowApplication.POSITIVE_POTENTIAL)
        #             node.SetSolutionStepValue(CompressiblePotentialFlowApplication.POSITIVE_POTENTIAL,0,potential)
        #             # node.SetSolutionStepValue(CompressiblePotentialFlowApplication.POSITIVE_POTENTIAL,0,1.0)

        #     elif n_inlet==0 and n_wall==0:
        #         cond.SetValue(CompressiblePotentialFlowApplication.VELOCITY_INFINITY, self.velocity_infinity)
        #         npos=0
        #         nneg=0
        #         for node in cond.GetNodes():
        #             distance=node.GetSolutionStepValue(CompressiblePotentialFlowApplication.WAKE_DISTANCE)
        #             if distance>0:
        #                 npos += 1
        #             elif distance<0:
        #                 nneg += 1
        #         if (npos>0 and nneg>0):
        #             cond.Set(KratosMultiphysics.STRUCTURE,True)
        #         else:
        #             cond.Set(KratosMultiphysics.STRUCTURE,False)
        #     else:
        #         zero_velocity=[0,0,0]
        #         cond.SetValue(CompressiblePotentialFlowApplication.VELOCITY_INFINITY, self.velocity_infinity)

        
    def ExecuteInitializeSolutionStep(self):
        self.Execute()
        
        
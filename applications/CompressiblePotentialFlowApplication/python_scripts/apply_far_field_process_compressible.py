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
            }  """ )


        settings.ValidateAndAssignDefaults(default_parameters);
        self.domain_model_part = Model.GetModelPart(settings["model_part_name"].GetString())
        self.velocity_infinity = KratosMultiphysics.Vector(3)#array('d', [1.0, 2.0, 3.14])#np.array([0,0,0])#np.zeros(3)#vector(3)
        self.main_model_part = self.domain_model_part.GetRootModelPart()
        self.velocity_infinity[0] = settings["velocity_infinity"][0].GetDouble()
        self.velocity_infinity[1] = settings["velocity_infinity"][1].GetDouble()
        self.velocity_infinity[2] = settings["velocity_infinity"][2].GetDouble()
                # For the elements
        print(self.velocity_infinity)
        self.main_model_part.GetProperties()[1].SetValue(
            CompressiblePotentialFlowApplication.VELOCITY_INFINITY, self.velocity_infinity)
        #self.density_infinity = settings["density_infinity"].GetDouble() #TODO: must read this from the properties
        self.inlet_phi = settings["inlet_phi"].GetDouble()
        self.domain_model_part.ProcessInfo.SetValue(CompressiblePotentialFlowApplication.VELOCITY_INFINITY,self.velocity_infinity)



    def Execute(self):
        for cond in self.domain_model_part.Conditions:
            cond.SetValue(CompressiblePotentialFlowApplication.VELOCITY_INFINITY, self.velocity_infinity)
            npos=0
            nneg=0

        for cond in self.domain_model_part.Conditions:
            normal=cond.GetValue(KratosMultiphysics.NORMAL)
            v_inf=cond.GetValue(CompressiblePotentialFlowApplication.VELOCITY_INFINITY)

            value = np.dot(normal,v_inf)

            if value<0:
                for node in cond.GetNodes():
                    inlet_phi=node.X*self.velocity_infinity[0] + node.Y*self.velocity_infinity[1] + node.Z*self.velocity_infinity[2]
                    node.Fix(CompressiblePotentialFlowApplication.VELOCITY_POTENTIAL)
                    node.Set(KratosMultiphysics.INLET)
                    node.SetSolutionStepValue(CompressiblePotentialFlowApplication.VELOCITY_POTENTIAL,0,inlet_phi)

        for node in self.main_model_part.Nodes:
            initial_phi=node.X*self.velocity_infinity[0] + node.Y*self.velocity_infinity[1] + node.Z*self.velocity_infinity[2]
            node.SetSolutionStepValue(CompressiblePotentialFlowApplication.VELOCITY_POTENTIAL,0,initial_phi)
            node.SetSolutionStepValue(CompressiblePotentialFlowApplication.AUXILIARY_VELOCITY_POTENTIAL,0,initial_phi)

    def ExecuteInitializeSolutionStep(self):
        self.Execute()


import KratosMultiphysics
import KratosMultiphysics.CompressiblePotentialFlowApplication as CompressiblePotentialFlowApplication
import numpy as np
import math
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
                "velocity_infinity": [3.4,0.0,0],
                "density_infinity"  : 1.0,
                "mach_infinity": 0.01,
                "gamma": 1.4,
                "pressure_infinity": 101325
            }  """ )


        settings.ValidateAndAssignDefaults(default_parameters);
        self.domain_model_part = Model.GetModelPart(settings["model_part_name"].GetString())
        self.fluid_model_part = self.domain_model_part.GetRootModelPart()
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
        self.domain_model_part.ProcessInfo.SetValue(
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

        for node in self.fluid_model_part.Nodes:
            initial_phi=node.X*self.velocity_infinity[0] + node.Y*self.velocity_infinity[1] + node.Z*self.velocity_infinity[2]
            node.SetSolutionStepValue(CompressiblePotentialFlowApplication.VELOCITY_POTENTIAL,0,initial_phi)
            node.SetSolutionStepValue(CompressiblePotentialFlowApplication.AUXILIARY_VELOCITY_POTENTIAL,0,initial_phi)

    def ExecuteInitializeSolutionStep(self):
        self.Execute()


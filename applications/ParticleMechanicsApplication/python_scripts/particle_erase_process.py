import KratosMultiphysics
import KratosMultiphysics.ParticleMechanicsApplication as KratosParticle

from math import sqrt

def Factory(settings, Model):
    if(not isinstance(settings, KratosMultiphysics.Parameters)):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ParticleEraseProcess(Model, settings["Parameters"])

## All the processes python should be derived from "Process"
class ParticleEraseProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings ):
        KratosMultiphysics.Process.__init__(self)

        default_parameters = KratosMultiphysics.Parameters( """
            {
                "model_part_name" : "PLEASE_CHOOSE_MODEL_PART_NAME",
                "eraser_type": "auto",
                "eraser_settings": {
                    "volume_center"     : [],
                    "sphere_radius"     : 0.0,
                    "bounding_box_size" : []
                }
            }  """ )

        settings.ValidateAndAssignDefaults(default_parameters)

        self.model_part = Model[settings["model_part_name"].GetString()]

        # Check type
        self.eraser_type = settings["eraser_type"].GetString()
        if (self.eraser_type in ["auto","sphere","bounding_box"]):
            if (self.eraser_type == "sphere"):
                self.volume_center = settings["eraser_settings"]["volume_center"].GetVector()
                self.sphere_radius = settings["eraser_settings"]["sphere_radius"].GetDouble()
            elif (self.eraser_type == "bounding_box"):
                self.volume_center = settings["eraser_settings"]["volume_center"].GetVector()
                self.bounding_box_size = settings["eraser_settings"]["bounding_box_size"].GetVector()
        else:
            err_msg =  "The requested type of particle erase feature: \"" + self.eraser_type + "\" is not available!\n"
            err_msg += "Available option is: \"auto\", \"sphere\", and \"bounding_box\"."
            raise Exception(err_msg)

        # Create process
        self.process = KratosParticle.ParticleEraseProcess(self.model_part)

    def ExecuteInitializeSolutionStep(self):
        ## Check
        if (self.eraser_type in ["auto","sphere","bounding_box"]):
            if (self.eraser_type == "sphere"):
                for mpm in self.model_part.Elements:
                    xg = mpm.GetValue(KratosParticle.MP_COORD)
                    dist = xg - self.volume_center
                    distance = sqrt(dist[0]*dist[0]+dist[1]*dist[1]+dist[2]*dist[2])
                    if (distance > self.sphere_radius):
                        mpm.Set(KratosMultiphysics.TO_ERASE, True)

                for mpc in self.model_part.Conditions:
                    xg = mpc.GetValue(KratosParticle.MPC_COORD)
                    dist = xg - self.volume_center
                    distance = sqrt(dist[0]*dist[0]+dist[1]*dist[1]+dist[2]*dist[2])
                    if (distance > self.sphere_radius):
                        mpc.Set(KratosMultiphysics.TO_ERASE, True)

            elif (self.eraser_type == "bounding_box"):
                bb_min = KratosMultiphysics.Vector(3)
                bb_max = KratosMultiphysics.Vector(3)
                for i in range(len(self.volume_center)):
                    bb_min[i] = self.volume_center[i] - 0.5 * self.bounding_box_size[i]
                    bb_max[i] = self.volume_center[i] + 0.5 * self.bounding_box_size[i]

                for mpm in self.model_part.Elements:
                    xg = mpm.GetValue(KratosParticle.MP_COORD)
                    for i in range(len(self.volume_center)):
                        if (xg[i] < bb_min[i] or xg[i] > bb_max[i]):
                            mpm.Set(KratosMultiphysics.TO_ERASE, True)
                            break

                for mpc in self.model_part.Conditions:
                    xg = mpc.GetValue(KratosParticle.MPC_COORD)
                    for i in range(len(self.volume_center)):
                        if (xg[i] < bb_min[i] or xg[i] > bb_max[i]):
                            mpc.Set(KratosMultiphysics.TO_ERASE, True)
                            break

        self.process.Execute()
from KratosMultiphysics import *
import KratosMultiphysics
from numpy import *
import itertools

def Factory(settings, Model):
    if( not isinstance(settings,KratosMultiphysics.Parameters) ):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ComputeLiftProcess(Model, settings["Parameters"])

##all the processes python processes should be derived from "python_process"
class ComputeLiftProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings ):
        KratosMultiphysics.Process.__init__(self)

        default_parameters = KratosMultiphysics.Parameters("""
            {
                "model_part_name":"PLEASE_CHOOSE_MODEL_PART_NAME",
                "upper_surface_model_part_name" : "please specify the model part that contains the upper surface nodes",
                "lower_surface_model_part_name" : "please specify the model part that contains the lower surface nodes",
                "mesh_id": 0,
                "velocity_infinity": [1.0,0.0,0],
                "reference_area": 1
            }  """)

        settings.ValidateAndAssignDefaults(default_parameters)

        self.upper_surface_model_part = Model[settings["upper_surface_model_part_name"].GetString()]
        self.lower_surface_model_part = Model[settings["lower_surface_model_part_name"].GetString()]
        self.velocity_infinity = [0,0,0]
        self.velocity_infinity[0] = settings["velocity_infinity"][0].GetDouble()
        self.velocity_infinity[1] = settings["velocity_infinity"][1].GetDouble()
        self.velocity_infinity[2] = settings["velocity_infinity"][2].GetDouble()
        self.reference_area =  settings["reference_area"].GetDouble()

    def ExecuteFinalize(self):
         print('COMPUTE LIFT')

         rx = 0.0
         ry = 0.0
         rz = 0.0

         for cond in itertools.chain(self.upper_surface_model_part.Conditions, self.lower_surface_model_part.Conditions):
           n = cond.GetValue(NORMAL)
           cp = cond.GetValue(PRESSURE)

           rx += n[0]*cp
           ry += n[1]*cp
           rz += n[2]*cp

         RZ = rz/self.reference_area
         RX = rx/self.reference_area
         RY = ry/self.reference_area

         Cl = RY
         Cd = RX

         print('RZ = ', RZ)

         print('Cl = ', Cl)
         print('Cd = ', Cd)
         print('Mach = ', self.velocity_infinity[0]/340)

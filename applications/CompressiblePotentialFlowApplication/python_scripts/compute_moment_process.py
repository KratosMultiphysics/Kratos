import KratosMultiphysics
import itertools
import numpy as np

def Factory(settings, Model):
    if( not isinstance(settings,KratosMultiphysics.Parameters) ):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ComputeMomentProcess(Model, settings["Parameters"])

##all the processes python processes should be derived from "python_process"
class ComputeMomentProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings ):
        KratosMultiphysics.Process.__init__(self)

        default_parameters = KratosMultiphysics.Parameters("""
            {
                "model_part_name":"PLEASE_CHOOSE_MODEL_PART_NAME",
                "upper_surface_model_part_name" : "please specify the model part that contains the upper surface nodes",
                "lower_surface_model_part_name" : "please specify the model part that contains the lower surface nodes",
                "reference_point" : [0.0,0.0,0.0],
                "mesh_id": 0,
                "velocity_infinity": [1.0,0.0,0.0],
                "reference_area": 1,
                "create_output_file": false
            }  """)

        settings.ValidateAndAssignDefaults(default_parameters)

        self.upper_surface_model_part = Model[settings["upper_surface_model_part_name"].GetString()]
        self.lower_surface_model_part = Model[settings["lower_surface_model_part_name"].GetString()]
        self.velocity_infinity = [0,0,0]
        self.velocity_infinity[0] = settings["velocity_infinity"][0].GetDouble()
        self.velocity_infinity[1] = settings["velocity_infinity"][1].GetDouble()
        self.velocity_infinity[2] = settings["velocity_infinity"][2].GetDouble()
        self.reference_area = settings["reference_area"].GetDouble()
        self.reference_point = np.array([0,0,0])
        self.create_output_file = settings["create_output_file"].GetBool()
        self.reference_point[0] = settings["reference_point"][0].GetDouble()
        self.reference_point[1] = settings["reference_point"][1].GetDouble()
        self.reference_point[2] = settings["reference_point"][2].GetDouble()

    def ExecuteFinalizeSolutionStep(self, refPoint = [0,0,0]):
         print('COMPUTE MOMENT')

         m = [0.0,0.0,0.0]

         for cond in itertools.chain(self.upper_surface_model_part.Conditions, self.lower_surface_model_part.Conditions):
           n =  np.array(cond.GetValue(KratosMultiphysics.NORMAL)) #Inigo: normal direction? (ass. inward of domain)
           cp = cond.GetValue(KratosMultiphysics.PRESSURE)
           mid_point = cond.GetGeometry().Center()
           mid_point =  np.array([mid_point.X, mid_point.Y, mid_point.Z])
           lever = mid_point-self.reference_point
           m += np.cross(lever, (-n*cp))

         Cm = m[2]/self.reference_area

         print('moment = ', m)
         print('Cm = ', Cm)
         print('Mach = ', self.velocity_infinity[0]/340)

         if self.create_output_file and False: #Todo
             with open("cl.dat", 'w') as cl_file:
                 cl_file.write('{0:15.12f}'.format(Cl))

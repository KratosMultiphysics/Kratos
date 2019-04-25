import KratosMultiphysics
from KratosMultiphysics import Logger

def crossProduct(A, B):
    C = KratosMultiphysics.Vector(3)
    C[0] = A[1]*B[2]-A[2]*B[1]
    C[1] = A[2]*B[0]-A[0]*B[2]
    C[2] = A[0]*B[1]-A[1]*B[0]
    return C

def Factory(settings, Model):
    if not isinstance(settings,KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ComputeMomentProcess(Model, settings["Parameters"])

##all the processes python processes should be derived from "python_process"
class ComputeMomentProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings ):
        KratosMultiphysics.Process.__init__(self)

        default_parameters = KratosMultiphysics.Parameters("""{
                "model_part_name": "please specify the model part that contains the surface nodes",
                "reference_point" : [0.0,0.0,0.0],
                "velocity_infinity": [1.0,0.0,0.0],
                "reference_area": 1,
                "create_output_file": false
            }  """)

        settings.ValidateAndAssignDefaults(default_parameters)

        self.body_model_part = Model[settings["model_part_name"].GetString()]
        self.velocity_infinity = KratosMultiphysics.Vector(3)
        self.velocity_infinity = settings["velocity_infinity"].GetVector()
        self.reference_area = settings["reference_area"].GetDouble()
        self.create_output_file = settings["create_output_file"].GetBool()
        self.reference_point = KratosMultiphysics.Vector(3)
        self.reference_point = settings["reference_point"].GetVector()

    def ExecuteFinalizeSolutionStep(self):
        Logger.PrintInfo('COMPUTE MOMENT')

        m = KratosMultiphysics.Vector(3)

        for cond in self.body_model_part.Conditions:
            n =  KratosMultiphysics.Vector(cond.GetValue(KratosMultiphysics.NORMAL)) #normal direction assumed outward of domain
            cp = cond.GetValue(KratosMultiphysics.PRESSURE)
            mid_point = cond.GetGeometry().Center()
            lever = mid_point-self.reference_point
            m += crossProduct(lever, n*(-cp))

        Cm = m[2]/self.reference_area

        Logger.PrintInfo('moment', m[2])
        Logger.PrintInfo('Cm', Cm)
        Logger.PrintInfo('Mach', self.velocity_infinity[0]/340)

        if self.create_output_file:
            with open("moment.dat", 'w') as mom_file:
                mom_file.write('{0:15.12f}'.format(Cm))

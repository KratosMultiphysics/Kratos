import KratosMultiphysics

def Factory(settings, Model):
    if( not isinstance(settings,KratosMultiphysics.Parameters) ):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ComputeLiftProcess(Model, settings["Parameters"])

##all the processes python processes should be derived from "python_process"
class ComputeLiftProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings ):
        KratosMultiphysics.Process.__init__(self)

        default_parameters = KratosMultiphysics.Parameters(r'''{
            "model_part_name": "please specify the model part that contains the surface nodes",
            "velocity_infinity": [1.0,0.0,0.0],
            "reference_area": 1.0,
            "create_output_file": false
        }''')

        settings.ValidateAndAssignDefaults(default_parameters)

        self.body_model_part = Model[settings["model_part_name"].GetString()]
        self.velocity_infinity = [0,0,0]
        self.velocity_infinity[0] = settings["velocity_infinity"][0].GetDouble()
        self.velocity_infinity[1] = settings["velocity_infinity"][1].GetDouble()
        self.velocity_infinity[2] = settings["velocity_infinity"][2].GetDouble()
        self.reference_area =  settings["reference_area"].GetDouble()
        self.create_output_file = settings["create_output_file"].GetBool()

    def ExecuteFinalizeSolutionStep(self):
        KratosMultiphysics.Logger.PrintInfo('ComputeLiftProcess','COMPUTE LIFT')

        rx = 0.0
        ry = 0.0
        rz = 0.0

        for cond in self.body_model_part.Conditions:
            n = cond.GetValue(KratosMultiphysics.NORMAL)
            cp = cond.GetValue(KratosMultiphysics.PRESSURE)

            rx += n[0]*cp
            ry += n[1]*cp
            rz += n[2]*cp

        RZ = rz/self.reference_area
        RX = rx/self.reference_area
        RY = ry/self.reference_area

        Cl = RY
        Cd = RX

        KratosMultiphysics.Logger.PrintInfo('ComputeLiftProcess','Cl = ', Cl)
        KratosMultiphysics.Logger.PrintInfo('ComputeLiftProcess','Cd = ', Cd)
        KratosMultiphysics.Logger.PrintInfo('ComputeLiftProcess','RZ = ', RZ)
        KratosMultiphysics.Logger.PrintInfo('ComputeLiftProcess','Mach = ', self.velocity_infinity[0]/340)

        if self.create_output_file:
            with open("cl.dat", 'w') as cl_file:
                cl_file.write('{0:15.12f}'.format(Cl))

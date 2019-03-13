from KratosMultiphysics import *
import KratosMultiphysics
from numpy import *
import itertools
import matplotlib.pyplot as plt
from exaqute.ExaquteTaskPyCOMPSs import *   # to exequte with pycompss
# from exaqute.ExaquteTaskHyperLoom import *  # to exequte with the IT4 scheduler
# from exaqute.ExaquteTaskLocal import *      # to execute with python3
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
                "reference_area": 1,
                "create_output_file": true,
                "path": "cp_distribution.dat"
            }  """)

        settings.ValidateAndAssignDefaults(default_parameters)
        self.upper_surface_model_part =Model.GetModelPart(settings["upper_surface_model_part_name"].GetString())
        self.lower_surface_model_part =Model.GetModelPart(settings["lower_surface_model_part_name"].GetString())
        for cond in itertools.chain(self.upper_surface_model_part.Conditions, self.lower_surface_model_part.Conditions):
            cond.Set(KratosMultiphysics.SOLID)
        self.fluid_model_part=self.upper_surface_model_part.GetRootModelPart()
        self.velocity_infinity = [0,0,0]
        self.velocity_infinity[0] = settings["velocity_infinity"][0].GetDouble()
        self.velocity_infinity[1] = settings["velocity_infinity"][1].GetDouble()
        self.velocity_infinity[2] = settings["velocity_infinity"][2].GetDouble()
        self.reference_area =  settings["reference_area"].GetDouble()
        self.create_output_file = settings["create_output_file"].GetBool()
        self.path  = settings["path"].GetString()

    def ExecuteFinalizeSolutionStep(self):
        print('COMPUTE LIFT')
        # self.process.Execute()
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
        self.fluid_model_part.SetValue(KratosMultiphysics.FRICTION_COEFFICIENT,Cl)

        print('Cl = ', Cl)
        print('Cd = ', Cd)
        print('RZ = ', RZ)
        print('Mach = ', self.velocity_infinity[0]/340)
        self.Output(self.path)


    @ExaquteTask(filepath=FILE_OUT)
    def Output(self,filepath):
            x_list=[]
            cp_list=[]

            for cond in self.upper_surface_model_part.Conditions:
                cp = cond.GetValue(PRESSURE)
                x = cond.GetGeometry().Center().X
                x_list.append(x)
                cp_list.append(cp)
            for cond in self.lower_surface_model_part.Conditions:
                cp = cond.GetValue(PRESSURE)
                x = cond.GetGeometry().Center().X
                x_list.append(x)
                cp_list.append(cp)  
            max_x=max(x_list)
            min_x=min(x_list)
            for i in range(0,len(x_list)):
                x_list[i]=(x_list[i]-min_x)/abs(max_x-min_x)

            with open(filepath, 'w') as cp_file:  
                for i in range(len(x_list)):
                    cp_file.write('%f %f\n' % (x_list[i], cp_list[i]))
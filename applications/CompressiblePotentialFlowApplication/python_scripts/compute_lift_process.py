from KratosMultiphysics import *
import KratosMultiphysics
from numpy import *
import itertools
import matplotlib.pyplot as plt
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
                "problem_name": "please specify problem name"
            }  """)

        settings.ValidateAndAssignDefaults(default_parameters)
        self.problem_name=settings["problem_name"].GetString()
        self.upper_surface_model_part =Model.GetModelPart(settings["upper_surface_model_part_name"].GetString()) 
        self.lower_surface_model_part =Model.GetModelPart(settings["lower_surface_model_part_name"].GetString())
        self.velocity_infinity = [0,0,0]
        self.velocity_infinity[0] = settings["velocity_infinity"][0].GetDouble()
        self.velocity_infinity[1] = settings["velocity_infinity"][1].GetDouble()
        self.velocity_infinity[2] = settings["velocity_infinity"][2].GetDouble()
        self.reference_area =  settings["reference_area"].GetDouble()
        
        NormalUpper=KratosMultiphysics.NormalCalculationUtils()
        NormalUpper.CalculateOnSimplex(self.upper_surface_model_part,2)
        NormalLower=KratosMultiphysics.NormalCalculationUtils()
        NormalLower.CalculateOnSimplex(self.lower_surface_model_part,2)


    def ExecuteFinalizeSolutionStep(self):
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
        x_upper=[]
        cp_upper=[]
        x_lower=[]
        cp_lower=[]
        for cond in self.upper_surface_model_part.Conditions:
            x_average=0
            cp = cond.GetValue(PRESSURE)
            for node in cond.GetNodes():
                x_average += 0.5*node.X
            x_upper.append(x_average)
            cp_upper.append(cp)
        for cond in self.lower_surface_model_part.Conditions:
            x_average=0
            cp = cond.GetValue(PRESSURE)
            for node in cond.GetNodes():
                x_average += 0.5*node.X
            x_lower.append(x_average)
            cp_lower.append(cp)  
        max_x=max(max(x_upper),max(x_lower))
        min_x=min(min(x_upper),min(x_lower))
        for i in range(0,len(x_upper)):
            x_upper[i]=(x_upper[i]-min_x)/abs(max_x-min_x)
        for i in range(0,len(x_lower)):
            x_lower[i]=(x_lower[i]-min_x)/abs(max_x-min_x) 

        plt.plot(x_upper,cp_upper,'o',x_lower,cp_lower,'ro')
        title="Cl: %.5f, Cd: %.5f" % (Cl,Cd)
        plt.title(title)
        plt.gca().invert_yaxis()
        plt.savefig('./Figures/'+self.problem_name+'.png', bbox_inches='tight')
        plt.close('all')
        print('RZ = ', RZ)

        print('Cl = ', Cl)
        print('Cd = ', Cd)
        print('Mach = ', self.velocity_infinity[0]/340)

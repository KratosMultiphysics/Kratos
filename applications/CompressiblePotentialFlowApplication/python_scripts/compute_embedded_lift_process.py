import KratosMultiphysics
import KratosMultiphysics.CompressiblePotentialFlowApplication
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
                
                "model_part_name" : "please specify the main model part",
                "reference_area": 1
            }  """)

        settings.ValidateAndAssignDefaults(default_parameters)
        self.main_model_part=Model.GetModelPart(settings["model_part_name"].GetString()).GetRootModelPart() 
        self.result_force=KratosMultiphysics.Vector(3)  
        self.process=KratosMultiphysics.CompressiblePotentialFlowApplication.ComputeEmbeddedLiftProcess(self.main_model_part,self.result_force)


    def ExecuteFinalizeSolutionStep(self):
        print("wip_compute_lift_level_set_process")
        self.process.Execute()
        x_upper=[]
        cp_upper=[]
        x_lower=[]
        cp_lower=[]
        for element in self.main_model_part.Elements:
            if element.Is(KratosMultiphysics.TO_SPLIT):
                gp=element.GetValue(KratosMultiphysics.BODY_FORCE)#provisional name
                pressure=element.GetValue(KratosMultiphysics.PRESSURE_COEFFICIENT)
                normal=element.GetValue(KratosMultiphysics.NORMAL)
                if normal[1]<=0:
                    x_upper.append(gp[0])
                    cp_upper.append(pressure)
                else:
                    x_lower.append(gp[0])
                    cp_lower.append(pressure)
        max_x=max(max(x_upper),max(x_lower))
        min_x=min(min(x_upper),min(x_lower))
        for i in range(0,len(x_upper)):
            x_upper[i]=(x_upper[i]-min_x)/abs(max_x-min_x)
        for i in range(0,len(x_lower)):
            x_lower[i]=(x_lower[i]-min_x)/abs(max_x-min_x)
        print("Cl:",self.result_force[1]) 
        print("Cd:",self.result_force[0])   
        plt.plot(x_upper,cp_upper,'o',label='Upper surface')
        plt.plot(x_lower,cp_lower,'ro',label='Lower surface')
        title="Cl: %.5f, Cd: %.5f" % (self.result_force[1],self.result_force[0])
        plt.title(title)
        plt.legend()
        plt.gca().invert_yaxis()
        plt.savefig('embedded_lift.png', bbox_inches='tight')
        plt.close('all')
    

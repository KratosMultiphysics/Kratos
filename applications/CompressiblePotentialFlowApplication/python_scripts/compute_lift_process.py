from KratosMultiphysics import *

import KratosMultiphysics

from array import array
from numpy import *
from math import *
import numpy as np
import itertools
#import matplotlib.pyplot as plt



def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ComputeLiftProcess(Model, settings["Parameters"])

##all the processes python processes should be derived from "python_process"
class ComputeLiftProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings ):
        KratosMultiphysics.Process.__init__(self) 
        
        default_parameters = KratosMultiphysics.Parameters( """
            {
                "model_part_name":"PLEASE_CHOOSE_MODEL_PART_NAME",
                "upper_surface_model_part_name" : "please specify the model part that contains the upper surface nodes",
                "lower_surface_model_part_name" : "please specify the model part that contains the lower surface nodes",
                "mesh_id": 0,
                "velocity_infinity": [1.0,0.0,0]
            }  """ );
        
        settings.ValidateAndAssignDefaults(default_parameters);
        
        
        #self.model_part = Model[settings["model_part_name"].GetString()]
        self.upper_surface_model_part = Model[settings["upper_surface_model_part_name"].GetString()]
        self.lower_surface_model_part = Model[settings["lower_surface_model_part_name"].GetString()]
        #self.model_part = Model.get('model_part_name',None)
        self.velocity_infinity = [0,0,0]#array('d', [1.0, 2.0, 3.14])#np.array([0,0,0])#np.zeros(3)#vector(3)
        self.velocity_infinity[0] = settings["velocity_infinity"][0].GetDouble()
        self.velocity_infinity[1] = settings["velocity_infinity"][1].GetDouble()
        self.velocity_infinity[2] = settings["velocity_infinity"][2].GetDouble()
        
        self.reference_area =  1#8#383 #m² # WRONG 489.89 m²   
        
        self.AOAdeg                 = 5#°
        
        #convert angle from degrees to radians
        self.AOArad = self.AOAdeg*pi/180  

  
    def ExecuteFinalize(self):
         print('COMPUTE LIFT')

         rx = 0.0
         ry = 0.0
         rz = 0.0

         for cond in itertools.chain(self.upper_surface_model_part.Conditions,self.lower_surface_model_part.Conditions):
           n = cond.GetValue(NORMAL)
           cp = cond.GetValue(PRESSURE)
           #print(cp)

           rx += n[0]*cp
           ry += n[1]*cp
           rz += n[2]*cp

         RZ = rz/self.reference_area
         RX = rx/self.reference_area
         RY = ry/self.reference_area
         
         Cl = RY*cos(self.AOArad) - RX*sin(self.AOArad)
         Cd = RY*sin(self.AOArad) + RX*cos(self.AOArad)
         
         print('RZ = ', RZ)
         print('RX = ', RX)
         print('RY = ', RY)
         
         print('Cl = ', Cl) 
         print('Cd = ', Cd)
         print('Mach = ', self.velocity_infinity[0]/340) 
         
         '''
         loads_file = open("loads.dat",'w') 
         loads_file.write("FULL POTENTIAL APPLICATION LOADS FILE\n\n")
         loads_file.write("UInf {0:13f} \n".format(self.velocity_infinity[0]))
         loads_file.write("Mach {0:13f} \n\n".format(self.velocity_infinity[0]/340))
         
         loads_file.write("Cl {0:15f} \n".format(Cl))
         loads_file.write("Cd {0:15f} \n\n".format(Cd))
         
         loads_file.write("Lift {0:15f} \n".format(Lift))
         loads_file.write("Drag {0:15f} \n".format(Drag))  
         
         loads_file.flush()
         '''
from KratosMultiphysics import *

import KratosMultiphysics

from array import array
from numpy import *
from math import *
import numpy as np
#import matplotlib.pyplot as plt



def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ComputeLiftProcess3D(Model, settings["Parameters"])

##all the processes python processes should be derived from "python_process"
class ComputeLiftProcess3D(KratosMultiphysics.Process):
    def __init__(self, Model, settings ):
        KratosMultiphysics.Process.__init__(self) 
        
        default_parameters = KratosMultiphysics.Parameters( """
            {
                "model_part_name"   : "PLEASE_CHOOSE_MODEL_PART_NAME",
                "skin_parts"        :["UpperSurface","LowerSurface","Fuselage"],
                "mesh_id"           : 0,
                "velocity_infinity" : [1.0,0.0,0]
            }  """ );
        
        settings.ValidateAndAssignDefaults(default_parameters);
        
        
        self.model_part = Model[settings["model_part_name"].GetString()]
        print('Creating aircraft modelpart...')
        self.aircraft_modelpart = self.model_part.CreateSubModelPart("aircraft_modelpart")
        for i in range(settings["skin_parts"].size()):
            mp = Model[settings["skin_parts"][i].GetString()]
            for node in mp.Nodes:
                self.aircraft_modelpart.Nodes.append(node)
            for cond in mp.Conditions:
                self.aircraft_modelpart.Conditions.append(cond)
        print('Finished creating aircraft modelpart')
        
        #KratosMultiphysics.NormalCalculationUtils().CalculateOnSimplex(self.aircraft_modelpart,self.model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE])  
        
        #self.model_part = Model.get('model_part_name',None)
        self.velocity_infinity = [0,0,0]#array('d', [1.0, 2.0, 3.14])#np.array([0,0,0])#np.zeros(3)#vector(3)
        self.velocity_infinity[0] = settings["velocity_infinity"][0].GetDouble()
        self.velocity_infinity[1] = settings["velocity_infinity"][1].GetDouble()
        self.velocity_infinity[2] = settings["velocity_infinity"][2].GetDouble()
        
        self.reference_area =  8#383 #m² # WRONG 489.89 m²   
        
        self.AOAdeg                 = 5#°
        
        #convert angle from degrees to radians
        self.AOArad = self.AOAdeg*pi/180  

  
    def ExecuteFinalize(self):
         print('COMPUTE LIFT')

         rx = 0.0
         ry = 0.0
         rz = 0.0
         
         counter = 1        
         for cond in self.aircraft_modelpart.Conditions:
           n = cond.GetValue(NORMAL)
           cp = cond.GetValue(PRESSURE)
           #print(cp)

           rx += n[0]*cp
           ry += n[1]*cp
           rz += n[2]*cp
           counter +=1

         
         print('Looped over ', counter, ' lifting conditions.')
         Lift = rz/self.reference_area
         Drag = rx/self.reference_area
         Side = ry/self.reference_area
         
         Cl = Lift*cos(self.AOArad) - Drag*sin(self.AOArad)
         Cd = Lift*sin(self.AOArad) + Drag*cos(self.AOArad)
         
         print('Lift = ', Lift)
         print('Drag = ', Drag)
         print('Side = ', Side)
         
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
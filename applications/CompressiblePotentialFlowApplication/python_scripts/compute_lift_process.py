from KratosMultiphysics import *

import KratosMultiphysics


from array import array
from numpy import *
from math import *
import numpy as np
import matplotlib.pyplot as plt



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
                "mesh_id": 0,
                "velocity_infinity": [1.0,0.0,0]
            }  """ );
        
        settings.ValidateAndAssignDefaults(default_parameters);
        
        self.model_part = Model[settings["model_part_name"].GetString()]
        #self.model_part = Model.get('model_part_name',None)
        self.velocity_infinity = [0,0,0]#array('d', [1.0, 2.0, 3.14])#np.array([0,0,0])#np.zeros(3)#vector(3)
        self.velocity_infinity[0] = settings["velocity_infinity"][0].GetDouble()
        self.velocity_infinity[1] = settings["velocity_infinity"][1].GetDouble()
        self.velocity_infinity[2] = settings["velocity_infinity"][2].GetDouble()
        
        #My stuff           
        #self.airfoillift_modelpart = self.model_part.CreateSubModelPart("airfoillift_modelpart")
        #mp = Model[settings["skin_parts"][3].GetString()]
        #for node in mp.Nodes:
        #   self.airfoillift_modelpart.Nodes.append(node)
        #for cond in mp.Conditions:
        #   self.airfoillift_modelpart.Conditions.append(cond)
        #   
        #self.upperwakelift_modelpart = self.model_part.CreateSubModelPart("upperwakelift_modelpart")
        #mp = Model[settings["skin_parts"][4].GetString()]
        #for node in mp.Nodes:
        #   self.upperwakelift_modelpart.Nodes.append(node)
        #for cond in mp.Conditions:
        #   self.upperwakelift_modelpart.Conditions.append(cond)
        

  
    def ExecuteFinalize(self):
        #compute the normal on the nodes of interest
         #KratosMultiphysics.NormalCalculationUtils().CalculateOnSimplex(self.airfoillift_modelpart, self.model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE])
         
         print('COMPUTE LIFT')
         
         #angle of attack
         alphadeg = 10.0
         
         #convert angle to radians
         alpharad = alphadeg*pi/180

         #center of rotation
         X0 = 1.0
         Y0 = 0.0
         
         rx = 0.0
         ry = 0.0
         rz = 0.0
         
         CP = []
         X = []
         
         cp_file = open("cp.dat",'w') 
         cp_file.write("FULL POTENTIAL APPLICATION CP FILE\n\n")
         cp_file.write("UInf {0:13f} \n".format(self.velocity_infinity[0]))
         #cp_file.write("PInf {0:13f} \n".format(self.pressure_infinity))
         #cp_file.write("DInf {0:13f} \n".format(self.density_infinity))
         cp_file.write("Mach {0:13f} \n\n".format(self.velocity_infinity[0]/340))
         
         cp_file.write("  pressure              cp               x               y               nx              ny              nz           length  \n")
         
         counter = 0
         for cond in self.model_part.Conditions:
           n = cond.GetValue(NORMAL)           
           norm = linalg.norm(n)
           n = n/norm
           edge = [ cond.GetNodes()[1].X0-cond.GetNodes()[0].X0 , cond.GetNodes()[1].Y0-cond.GetNodes()[0].Y0 , cond.GetNodes()[1].Z0-cond.GetNodes()[0].Z0]
           length = linalg.norm(edge)
           #print(cond)
           
           # Previous way of computing pressure
           #calc_pressure = CalculatePressure(cond)
           #pressure = calc_pressure.Execute()
           #print(pressure, cond.GetValue(PRESSURE))
           
           
           
           #pressure = 1
           counter +=1
           pressure = cond.GetValue(PRESSURE)
           
           
           cp = pressure# 2*(pressure - self.pressure_infinity)/(self.density_infinity*self.velocity_infinity[0]*self.velocity_infinity[0])
           CP.append(cp)
           x = 0.5*(cond.GetNodes()[1].X0+cond.GetNodes()[0].X0)
           y = 0.5*(cond.GetNodes()[1].Y0+cond.GetNodes()[0].Y0)
           
           #rotation of the points
           #translation
           x = x - X0
           y = y - Y0
           
           #rotation
           xnew =  x*cos(alpharad) - y*sin(alpharad)
           ynew =  x*sin(alpharad) + y*cos(alpharad)  
            
           #tranlation
           x = xnew + X0
           y = ynew + Y0           
           
           X.append(x)
           
           #if cond.GetNodes()[0].X0 > 0.6 and cond.GetNodes()[0].X0 < 0.8 and cond.GetNodes()[0].Y0 > 0:
           #   print(' ')
           #   print('X0 = ', cond.GetNodes()[0].X0, 'Y0 = ', cond.GetNodes()[0].Y0)
           #   print('X1 = ', cond.GetNodes()[1].X0, 'Y1 = ', cond.GetNodes()[1].Y0)
           #   print('n =', n)
           #   print('length =', length)
           #   print('pressure =', pressure)
           #   print('rx = ', n[0]*pressure*length, 'ry = ', n[1]*pressure*length, 'rz = ', n[2]*pressure*length)
           #   if cond.GetNodes()[0].Y0 < 0.0001 or cond.GetNodes()[1].Y0 < 0.0001:
           #     print('LOWERCAMBER')
           #   else:
           #     print('UPPERCAMBER')
           #   print('X0 = ', cond.GetNodes()[0].X0, 'Y0 = ', cond.GetNodes()[0].Y0)
           #   print('X1 = ', cond.GetNodes()[1].X0, 'Y1 = ', cond.GetNodes()[1].Y0)
           #   print('Pressure = ', pressure)
           #   print('Cp = ', cp)
           #   print(' ')
           
           
           rx += n[0]*pressure*length
           ry += n[1]*pressure*length
           rz += n[2]*pressure*length
           #print(ry)
           #print("  counter =", counter,"   X =", cond.GetNodes()[0].X0 ,"  An =", n*length ,"  pressure =", pressure ,"    lift =", ry)
           
           cp_file.write('{0:13f} {1:15f} {2:15f} {3:15f} {4:15f} {5:15f} {6:15f} {7:15f}\n'.format(pressure, cp, x, y, n[0], n[1], n[2], length))
           
           #print(cond.GetNodes()[0])
           #print(n)
           #print(ry)
         
         cp_file.flush()
         
         Lift = ry
         Drag = rx
         
         c = 1 #chord
         Cl = Lift #2*Lift/(self.density_infinity*self.velocity_infinity[0]*self.velocity_infinity[0]*c)
         Cd = Drag #2*Drag/(self.density_infinity*self.velocity_infinity[0]*self.velocity_infinity[0]*c)        
         
         
         loads_file = open("loads.dat",'w') 
         loads_file.write("FULL POTENTIAL APPLICATION LOADS FILE\n\n")
         loads_file.write("UInf {0:13f} \n".format(self.velocity_infinity[0]))
         #loads_file.write("PInf {0:13f} \n".format(self.pressure_infinity))
         #loads_file.write("DInf {0:13f} \n".format(self.density_infinity))
         loads_file.write("Mach {0:13f} \n\n".format(self.velocity_infinity[0]/340))
         
         loads_file.write("Cl {0:15f} \n".format(Cl))
         loads_file.write("Cd {0:15f} \n\n".format(Cd))
         
         loads_file.write("Lift {0:15f} \n".format(Lift))
         loads_file.write("Drag {0:15f} \n".format(Drag))        
        
         
         loads_file.flush()
         
         print('Lift = ', Lift)
         print('Drag = ', Drag)
         
         print('Cl = ', Cl) 
         print('Cd = ', Cd)
         print('Mach = ', self.velocity_infinity[0]/340)        
        
         
         plt.plot(X,CP,'ro')
         plt.title('alpha = 15 deg', fontsize=25)
         plt.ylabel(r'$c_p$', fontsize=25)
         plt.xlabel(r'$x/c$', fontsize=25)
         #plt.axis([0,1,-0.6,1.2])
         plt.tick_params(axis='x', labelsize=20)
         plt.tick_params(axis='y', labelsize=20)
         #plt.xticks(np.arange(min(X), max(X)+0.05, 0.05))
         plt.gca().invert_yaxis()
         plt.grid(True)         
         plt.show()
         
         
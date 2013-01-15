from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *
#from KratosMultiphysics.MetisApplication import *
#from KratosMultiphysics.mpi import *

from numpy import *

from DEM_explicit_solver_var import *


#PRESSURE CALCULATION

def ApplyPressure(Pressure,model_part,solver,SKIN,BOT,TOP,LAT,XLAT,XBOT,XBOTCORNER,XTOP,XTOPCORNER): 

  print("")
  print("Applying Pressure")
  print("")
  
  
  skin_list = list()
  top_nodes_list = list()
  bot_nodes_list = list()
  total_cross_section = 0.0
  
  #Cylinder dimensions
  
  h   = 0.3
  d   = 0.15
  eps = 2
  
  surface = 2*(3.141592*d*d*0.25)+(3.141592*d*h)
  
  top_pressure = 0.0
  bot_pressure = 0.0
  
 
  for element in model_part.Elements:
 
    node = element.GetNode(0)
    r = node.GetSolutionStepValue(RADIUS,0)
    x = node.X
    y = node.Y
    z = node.Z
    
    values = Array3()
    values[0] = 0.0
    values[1] = 0.0
    values[2] = 0.0
      
    cross_section = 3.141592*r*r
    
    vect = zeros(3, double) 

    if ( (x*x+z*z)>=((d/2-eps*r)*(d/2-eps*r)) ): 
 
      skin_list.append(element)
      total_cross_section = total_cross_section + cross_section 
      
      #vector normal al centre:
      vect_moduli = sqrt(x*x+z*z)
      #print(vect_moduli)
      if(vect_moduli>0.0):
        vect[0]=-x/vect_moduli
        vect[1]=0
        vect[2]=-z/vect_moduli
      
      values[0]=cross_section*Pressure*vect[0]
      values[1]= 0.0
      values[2]=cross_section*Pressure*vect[2]
      
      
    if ( (y<=eps*r ) or (y>=(h-eps*r)) ): 
   
        values[0]=0.0
        values[2]=0.0
        if ( y>h/2 ):
            values[1]=-cross_section*Pressure
        else:
            values[1]= cross_section*Pressure

        if ( (x*x+z*z) >= ((d/2-eps*r)*(d/2-eps*r) ) ) :
            #vector normal al centre:
            vect_moduli = sqrt(x*x+z*z)
            
            if ( vect_moduli>0.0 ) :
                vect[0]=-x/vect_moduli
                vect[1]=0.0
                vect[2]=-z/vect_moduli

            values[0]= cross_section*Pressure*vect[0]*0.70710678
            values[2]= cross_section*Pressure*vect[2]*0.70710678
            
            if ( y>h/2 ):
                values[1]=-cross_section*Pressure*0.70710678
                top_pressure = top_pressure + values[1]
                top_nodes_list.append(node)
                
            else:
                values[1]=cross_section*Pressure*0.70710678 
                bot_pressure = bot_pressure + values[1]
                bot_nodes_list.append(node)
        else:
            skin_list.append(element)
            total_cross_section = total_cross_section + cross_section 

    node.SetSolutionStepValue(EXTERNAL_APPLIED_FORCE,values)

  correction_factor = surface/total_cross_section
  
  print("Correction Factor: ")
  print( correction_factor )
  
  top_bot_mean = 0.5*(top_pressure + bot_pressure)
  
  top_corrected = top_pressure + top_bot_mean
  bot_corrected = bot_pressure + top_bot_mean
  
  residual = top_corrected+bot_corrected
  
  print(top_pressure)
  print(bot_pressure)
  print(residual)
  
  #rebalancing the total pressure before applying:
  
  for element in skin_list:
    
    node = element.GetNode(0)

    values[0] = correction_factor*node.GetSolutionStepValue(EXTERNAL_APPLIED_FORCE)[0]
    values[1] = correction_factor*node.GetSolutionStepValue(EXTERNAL_APPLIED_FORCE)[1]
    values[2] = correction_factor*node.GetSolutionStepValue(EXTERNAL_APPLIED_FORCE)[2]
 
    node.SetSolutionStepValue(EXTERNAL_APPLIED_FORCE,values)
  
  print("End Applying Imposed Forces")
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  #Dimensions
  
  h   = 0.3
  d   = 0.15
  eps = 2

  surface = 2*(3.141592*d*d*0.25)+(3.141592*d*h)
  
  #PRESSURE ON TOP
  #for node in TOP:
        
        #r = node.GetSolutionStepValue(RADIUS,0)
        #x = node.X
        #y = node.Y
        #z = node.Z
   
   
   
   
   
   
   
   
   
        #values = Array3()
        #values[0] = 0.0
        #values[1] = 0.0
        #values[2] = 0.0
          
        #cross_section = 3.141592*r*r
        
        #vect = zeros(3, double) 

        #if ( (x*x+z*z)>=((d/2-eps*r)*(d/2-eps*r)) ): 
          
          #element.SetValue(SKIN_SPHERE,1)
          #skin_list.append(element)
          #total_cross_section = total_cross_section + cross_section 
          
          ##vector normal al centre:
          #vect_moduli = sqrt(x*x+z*z)
          ##print(vect_moduli)
          #if(vect_moduli>0.0):
            #vect[0]=-x/vect_moduli
            #vect[1]=0
            #vect[2]=-z/vect_moduli
          
          #values[0]=cross_section*Pressure*vect[0]
          #values[1]= 0.0
          #values[2]=cross_section*Pressure*vect[2]
          
          
        #if ( (y<=eps*r ) or (y>=(h-eps*r)) ): 

            #element.SetValue(SKIN_SPHERE,1)
            ##vector normal al centre:      
            #values[0]=0.0
            #values[2]=0.0
            #if ( y>h/2 ):
                #values[1]=-cross_section*Pressure
            #else:
                #values[1]= cross_section*Pressure
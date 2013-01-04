#PRESSURE CALCULATION

def ApplyPressure(model_part,solver,SKIN,BOT,TOP,LAT,XLAT,XBOT,XBOTCORNER,XTOP,XTOPCORNER): 

  print("")
  print("Applying Pressure")
  print("")
  
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
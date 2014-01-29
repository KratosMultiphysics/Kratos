from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *
# from KratosMultiphysics.MetisApplication import *
# from KratosMultiphysics.mpi import *

import math
from numpy import *

bond_00_05 = list(); bond_05_10 = list(); bond_10_15 = list(); bond_15_20 = list() ;bond_20_25 = list(); bond_25_30 = list(); bond_30_35 = list(); bond_35_40 = list(); bond_40_45 = list(); bond_45_50 = list(); bond_50_55 = list(); bond_55_60 = list(); bond_60_65 = list(); bond_65_70 = list(); bond_70_75 = list(); bond_75_80 = list(); bond_80_85 = list(); bond_85_90 = list();

sizes = [];
for i in xrange(0,18):
    sizes.append(0.0);
    
# PRESSURE CALCULATION

def OrientationStudy(contact_model_part):
  
  #OrientationChart = open("OrientationChart", 'w')
    
  #OrientationChart.write(str(j)+"    "+str(alpha_deg)+'\n')
    
  for element in contact_model_part.Elements:

    u1 = element.GetNode(1).X - element.GetNode(0).X
    u2 = element.GetNode(1).Y - element.GetNode(0).Y
    u3 = element.GetNode(1).Z - element.GetNode(0).Z
     
    alpha = abs(math.asin(abs(u2)/math.sqrt((u1*u1)+(u2*u2)+(u3*u3))))
    
    alpha_deg = alpha/math.pi*180

    if(alpha_deg >= 0.0 and alpha_deg < 5.0):
      bond_00_05.append(element)

    if(alpha_deg >= 5.0 and alpha_deg < 10.0):
      bond_05_10.append(element)
      
    if(alpha_deg >= 10.0 and alpha_deg < 15.0):
      bond_10_15.append(element)
      
    if(alpha_deg >= 15.0 and alpha_deg < 20.0):
      bond_15_20.append(element)
      
    if(alpha_deg >= 20.0 and alpha_deg < 25.0):
      bond_20_25.append(element)
      
    if(alpha_deg >= 25.0 and alpha_deg < 30.0):
      bond_25_30.append(element)
      
    if(alpha_deg >= 30.0 and alpha_deg < 35.0):
      bond_30_35.append(element)
    
    if(alpha_deg >= 35.0 and alpha_deg < 40.0):
      bond_35_40.append(element)
    
    if(alpha_deg >= 40.0 and alpha_deg < 45.0):
      bond_40_45.append(element)
      
    if(alpha_deg >= 45.0 and alpha_deg < 50.0):
      bond_45_50.append(element)
      
    if(alpha_deg >= 50.0 and alpha_deg < 55.0):
      bond_50_55.append(element)
      
    if(alpha_deg >= 55.0 and alpha_deg < 60.0):
      bond_55_60.append(element)
    
    if(alpha_deg >= 60.0 and alpha_deg < 65.0):
      bond_60_65.append(element)
      
    if(alpha_deg >= 65.0 and alpha_deg < 70.0):
      bond_65_70.append(element)
      
    if(alpha_deg >= 70.0 and alpha_deg < 75.0):
      bond_70_75.append(element)
          
    if(alpha_deg >= 75.0 and alpha_deg < 80.0):
      bond_75_80.append(element)
      
    if(alpha_deg >= 80.0 and alpha_deg < 85.0):
      bond_80_85.append(element)
      
    if(alpha_deg >= 85.0 and alpha_deg < 90.0):
      bond_85_90.append(element)
   
  i=0 
  for item in [bond_00_05, bond_05_10, bond_10_15, bond_15_20, bond_20_25, bond_25_30, bond_30_35, bond_35_40, bond_40_45,  bond_45_50, bond_50_55, bond_55_60, bond_60_65, bond_65_70, bond_70_75, bond_75_80, bond_80_85, bond_85_90]:
    
    sizes[i] = len(item)
    i+=1    

  
  #OrientationChart.close()
      
def ApplyPressure(Pressure, XLAT, XBOT, XTOP, XBOTCORNER, XTOPCORNER, alpha_top, alpha_bot, alpha_lat):

    for node in XLAT:

        r = node.GetSolutionStepValue(RADIUS)
        x = node.X
        y = node.Y
        z = node.Z

        values = Array3()
        values[0] = 0.0
        values[1] = 0.0
        values[2] = 0.0

        cross_section = 3.141592 * r * r

        vect = zeros(3, double)

        # vector normal al centre:
        vect_moduli = math.sqrt(x * x + z * z)

        if(vect_moduli > 0.0):
            vect[0] = -x / vect_moduli
            vect[1] = 0
            vect[2] = -z / vect_moduli

        values[0] = cross_section * alpha_lat * Pressure * vect[0]
        values[1] = 0.0
        values[2] = cross_section * alpha_lat * Pressure * vect[2]

        node.SetSolutionStepValue(EXTERNAL_APPLIED_FORCE, values)

    for node in XTOP:

        r = node.GetSolutionStepValue(RADIUS)
        x = node.X
        y = node.Y
        z = node.Z

        values = Array3()
        values[0] = 0.0
        values[1] = 0.0
        values[2] = 0.0

        cross_section = 3.141592 * r * r

        values[0] = 0.0
        values[1] = -cross_section * alpha_top * Pressure
        values[2] = 0.0

        node.SetSolutionStepValue(EXTERNAL_APPLIED_FORCE, values)

    for node in XBOT:

        r = node.GetSolutionStepValue(RADIUS)
        x = node.X
        y = node.Y
        z = node.Z

        values = Array3()
        values[0] = 0.0
        values[1] = 0.0
        values[2] = 0.0

        cross_section = 3.141592 * r * r

        values[0] = 0.0
        values[1] = cross_section * alpha_bot * Pressure
        values[2] = 0.0

        node.SetSolutionStepValue(EXTERNAL_APPLIED_FORCE, values)

    for node in XTOPCORNER:

        r = node.GetSolutionStepValue(RADIUS)
        x = node.X
        y = node.Y
        z = node.Z

        values = Array3()
        values[0] = 0.0
        values[1] = 0.0
        values[2] = 0.0

        cross_section = 3.141592 * r * r

        vect = zeros(3, double)

        # vector normal al centre:
        vect_moduli = math.sqrt(x * x + z * z)

        if(vect_moduli > 0.0):
            vect[0] = -x / vect_moduli
            vect[1] = 0
            vect[2] = -z / vect_moduli

        values[0] = cross_section * alpha_lat * Pressure * vect[0] * 0.70710678
        values[1] = -cross_section * alpha_top * Pressure * 0.70710678
        values[2] = cross_section * alpha_lat * Pressure * vect[2] * 0.70710678

        node.SetSolutionStepValue(EXTERNAL_APPLIED_FORCE, values)

    for node in XBOTCORNER:

        r = node.GetSolutionStepValue(RADIUS)
        x = node.X
        y = node.Y
        z = node.Z

        values = Array3()
        values[0] = 0.0
        values[1] = 0.0
        values[2] = 0.0

        cross_section = 3.141592 * r * r

        vect = zeros(3, double)

        # vector normal al centre:
        vect_moduli = math.sqrt(x * x + z * z)

        if(vect_moduli > 0.0):
            vect[0] = -x / vect_moduli
            vect[1] = 0
            vect[2] = -z / vect_moduli

        values[0] = cross_section * alpha_lat * Pressure * vect[0] * 0.70710678
        values[1] = cross_section * alpha_bot * Pressure * 0.70710678
        values[2] = cross_section * alpha_lat * Pressure * vect[2] * 0.70710678

        node.SetSolutionStepValue(EXTERNAL_APPLIED_FORCE, values)


def ApplyLateralPressure(Pressure, XLAT, XBOT, XTOP, XBOTCORNER, XTOPCORNER, alpha_top, alpha_bot, alpha_lat):

    for node in XLAT:

        r = node.GetSolutionStepValue(RADIUS)
        x = node.X
        y = node.Y
        z = node.Z

        values = Array3()
        values[0] = 0.0
        values[1] = 0.0
        values[2] = 0.0

        cross_section = 3.141592 * r * r

        vect = zeros(3, double)

        # vector normal al centre:
        vect_moduli = math.sqrt(x * x + z * z)

        if(vect_moduli > 0.0):
            vect[0] = -x / vect_moduli
            vect[1] = 0
            vect[2] = -z / vect_moduli

        values[0] = cross_section * alpha_lat * Pressure * vect[0]
        values[1] = 0.0
        values[2] = cross_section * alpha_lat * Pressure * vect[2]

        node.SetSolutionStepValue(EXTERNAL_APPLIED_FORCE, values)

    for node in XTOPCORNER:

        r = node.GetSolutionStepValue(RADIUS)
        x = node.X
        y = node.Y
        z = node.Z

        values = Array3()
        values[0] = 0.0
        values[1] = 0.0
        values[2] = 0.0

        cross_section = 3.141592 * r * r

        vect = zeros(3, double)

        # vector normal al centre:
        vect_moduli = math.sqrt(x * x + z * z)

        if(vect_moduli > 0.0):
            vect[0] = -x / vect_moduli
            vect[1] = 0
            vect[2] = -z / vect_moduli

        values[0] = cross_section * alpha_lat * Pressure * vect[0] * 0.70710678
        values[1] = 0.0
        values[2] = cross_section * alpha_lat * Pressure * vect[2] * 0.70710678

        node.SetSolutionStepValue(EXTERNAL_APPLIED_FORCE, values)

    for node in XBOTCORNER:

        r = node.GetSolutionStepValue(RADIUS)
        x = node.X
        y = node.Y
        z = node.Z

        values = Array3()
        values[0] = 0.0
        values[1] = 0.0
        values[2] = 0.0

        cross_section = 3.141592 * r * r

        vect = zeros(3, double)

        # vector normal al centre:
        vect_moduli = math.sqrt(x * x + z * z)

        if(vect_moduli > 0.0):
            vect[0] = -x / vect_moduli
            vect[1] = 0
            vect[2] = -z / vect_moduli

        values[0] = cross_section * alpha_lat * Pressure * vect[0] * 0.70710678
        values[1] = 0.0
        values[2] = cross_section * alpha_lat * Pressure * vect[2] * 0.70710678

        node.SetSolutionStepValue(EXTERNAL_APPLIED_FORCE, values)

def MeasureRadialStrain(LAT):
    
   mean_radial_strain = 0.0
   weight = 0.0
   
   for node in LAT:
     
     r = node.GetSolutionStepValue(RADIUS)
     x = node.X
     z = node.Z
     
     x0 = node.X0
     z0 = node.Z0
     
     dist_initial = math.sqrt(x0 * x0 + z0 * z0)
     dist_now = math.sqrt(x * x + z * z)
     
     node_radial_strain = (dist_now - dist_initial) / dist_initial
     
     mean_radial_strain += node_radial_strain*r*r
     weight += r*r
   
   if(weight == 0.0):
     print ("Error in MeasureRadialStrain. Lateral skin particles not well defined")
   else:
     
    radial_strain = mean_radial_strain/weight
   
   return radial_strain
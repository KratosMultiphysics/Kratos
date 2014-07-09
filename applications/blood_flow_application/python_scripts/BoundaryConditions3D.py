from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from KratosMultiphysics import *
from KratosMultiphysics.BloodFlowApplication import *
CheckForPreviousImport()

import math
import time
import sys
import config
import simulation_config

# File: 11/12/2013


class TransferTools:



    def __init__(self,model_part_3d):
		self.model_part_3d = model_part_3d
	
    def Initialize(self):
      self.inlets_3d = []
      self.outlets_3d = []
      self.skin_3d = []
      self.inlet_areas_3d = []
      self.outlet_areas_3d = []
      self.inlet_velocity_directions = []
      self.flow_3d_in = 0.0
      self.flow_3d_in_pressure = 0.0
      self.flow_3d_out = 0.0
      # compute normals and 3d areas
      Coupling_version="BoundaryConditions3D"
      print(Coupling_version)
 
      #ThrowErrors = False
      #check_process = TetrahedralMeshOrientationCheck(self.model_part_3d,ThrowErrors)
      #check_process.Execute()  
      
      normal = simulation_config.normal_direcction
      if (normal == "in"):
      # Normal estan apuntando hacia el interior con esta funcion hacemos que las normales apunten hacia el exterior
      # asegurar hacia donde apuntan las noramles debido al signo de la velocidad de entrada -1 en trasfer 1D3D
	NormalCalculationUtils().SwapNormals(self.model_part_3d)
	self.Velocity_Sig= -1
	print ("Switch the normal")
      else:
	self.Velocity_Sig= -1
        
      BodyNormalCalculationUtils().CalculateBodyNormals(self.model_part_3d, 3)
      # detect 3d inlets
      #self.normal_tools.CalculateBodyNormals(self.model_part, self.domain_size);
     
      Reference_Radius = 0.0
	
      for i in range(0,1):
	nfound = 0
	aux = []
	for node in self.model_part_3d.Nodes:
		if (node.GetSolutionStepValue(FLAG_VARIABLE) == i):
			aux.append(node)		
	if(len(aux) != 0):
		self.skin_3d.append(aux)
		#print (str(node.Id))
	else:
		break
      
      # detect 3d inlets
      normN = 0.0
      tmp = []
      #print ("2")
      for i in range(100, 101):
		nfound = 0
		aux = []
		directions = []
		normN_promedio = 0.0
		Norm_Mean = []
		Norm_Mean = [0.0, 0.0, 0.0]
		tmp = []
		total_nodes = 0.0 
		for node in self.model_part_3d.Nodes:
		      tmp = 0.0
		      if (node.GetSolutionStepValue(FLAG_VARIABLE) == i):
			      tmp = node.GetSolutionStepValue(NORMAL)
			      normN = math.sqrt(tmp[0] ** 2 + tmp[1] ** 2 + tmp[2] ** 2)	      
			      tmp[0] = tmp[0]/normN
			      tmp[1] = tmp[1]/normN
			      tmp[2] = tmp[2]/normN    
			      Norm_Mean[0]=Norm_Mean[0]+tmp[0]
			      Norm_Mean[1]=Norm_Mean[1]+tmp[1]
			      Norm_Mean[2]=Norm_Mean[2]+tmp[2]
			      total_nodes=total_nodes+1
			      
		Norm_Mean[0]=Norm_Mean[0]/total_nodes
		Norm_Mean[1]=Norm_Mean[1]/total_nodes
		Norm_Mean[2]=Norm_Mean[2]/total_nodes	  
		normN_promedio = math.sqrt(Norm_Mean[0] ** 2 + Norm_Mean[1] ** 2 + Norm_Mean[2] ** 2)
		Norm_Mean[0] =  Norm_Mean[0]/normN_promedio
		Norm_Mean[1] =  Norm_Mean[1]/normN_promedio
		Norm_Mean[2] =  Norm_Mean[2]/normN_promedio

		for node in self.model_part_3d.Nodes:
			tmp = 0.0
			if (node.GetSolutionStepValue(FLAG_VARIABLE) == i):
			      aux.append(node)
			      node.Fix(VELOCITY_X)
			      node.Fix(VELOCITY_Y)
			      node.Fix(VELOCITY_Z)
			      tmp = node.GetSolutionStepValue(NORMAL)
			      normN = math.sqrt(tmp[0] ** 2 + tmp[1] ** 2 + tmp[2] ** 2)
			      #print (str(tmp))
			      #print (str(normN))
			      #tmp = tmp/normN
			      #tmp[0] =  (tmp[0]*Norm_Mean[0])/normN
			      #tmp[1] =  (tmp[1]*Norm_Mean[1])/normN
			      #tmp[2] =  (tmp[2]*Norm_Mean[2])/normN	
			      tmp[0] =  Norm_Mean[0]
			      tmp[1] =  Norm_Mean[1]
			      tmp[2] =  Norm_Mean[2]
			      directions.append(tmp)
			      #print(tmp)	
		#print(node.Id)
		#print (i)
		#print ("3")
		#print "3D_inlet has inlet_nodes_3d = self.inlets_3d[i]been assigned", node

		if(len(aux) != 0):
			  self.inlets_3d.append(aux)
			  self.inlet_velocity_directions.append(directions)
		else:
			  break
      #raw_input()
         
      if(len(self.inlets_3d) == 0):
	  sys.exit("number of 3d_inlets are zero!! Please check your config.py file!")
	  # print "number of 3d is zero"
	  # err


        # detect 3d outlets
      for i in range(1001, 1199):
	  # print i
	  nfound = 0
	  aux = []
	  for node in self.model_part_3d.Nodes:
	      if (node.GetSolutionStepValue(FLAG_VARIABLE) == i):
		  aux.append(node)
		  node.Fix(PRESSURE)
		  # print "3D_outlet has been assigned(H)", node.Id
	  if(len(aux) != 0):
	      self.outlets_3d.append(aux)
	  else:
	      break

      if(len(self.outlets_3d) == 0):
	  # print "number of 1d is zero"
	  sys.exit("number of 3d_outlets are zero!! Please check your config.py file!")

      Total_area_inlet = 0.0
      for i in range(0, len(self.inlets_3d)):
	  inlet_nodes_3d = self.inlets_3d[i]
	  area3d = 0.0
	  normN_promedio = 0.0
	  Norm_Mean = []
	  Norm_Mean = [0.0, 0.0, 0.0]
	  tmp = []
	  tmp = 0.0
	  total_nodes = 0.0 
	  for node in inlet_nodes_3d:
	      tmp = node.GetSolutionStepValue(NORMAL)
	      normN = math.sqrt(tmp[0] ** 2 + tmp[1] ** 2 + tmp[2] ** 2)	      
	      tmp[0] = tmp[0]/normN
	      tmp[1] = tmp[1]/normN
	      tmp[2] = tmp[2]/normN    
	      Norm_Mean[0]=Norm_Mean[0]+tmp[0]
	      Norm_Mean[1]=Norm_Mean[1]+tmp[1]
	      Norm_Mean[2]=Norm_Mean[2]+tmp[2]
	      total_nodes=total_nodes+1	  
	  Norm_Mean[0]=Norm_Mean[0]/total_nodes
	  Norm_Mean[1]=Norm_Mean[1]/total_nodes
	  Norm_Mean[2]=Norm_Mean[2]/total_nodes	  
	  normN_promedio = math.sqrt(Norm_Mean[0] ** 2 + Norm_Mean[1] ** 2 + Norm_Mean[2] ** 2)
	  Norm_Mean[0] =  Norm_Mean[0]/normN_promedio
	  Norm_Mean[1] =  Norm_Mean[1]/normN_promedio
	  Norm_Mean[2] =  Norm_Mean[2]/normN_promedio
	  print(Norm_Mean)
	  for node in inlet_nodes_3d:
	      n = node.GetSolutionStepValue(NORMAL)
	      #a = math.sqrt(n[0] ** 2 + n[1] ** 2 + n[2] ** 2) Direct Method, up-stimation
	      a = math.sqrt((n[0]*Norm_Mean[0]) ** 2 + (n[1]*Norm_Mean[1]) ** 2 + (n[2]*Norm_Mean[2]) ** 2)
	      area3d += a
	      
	  self.inlet_areas_3d.append(area3d)	  	   
	  Total_area_inlet = Total_area_inlet + area3d
	  print("area3d_inlet", Total_area_inlet)	  
      Total_area_inlet = Total_area_inlet / len(self.inlets_3d)
            
      Total_area_outlet = 0.0      
      for i in range(0, len(self.outlets_3d)):
	  inlet_nodes_3d = self.outlets_3d[i]
	  area3d = 0.0
	  normN_promedio = 0.0
	  Norm_Mean = []
	  Norm_Mean = [0.0, 0.0, 0.0]
	  tmp = []
	  tmp = 0.0
	  total_nodes = 0.0 
	  for node in inlet_nodes_3d:
	      tmp = node.GetSolutionStepValue(NORMAL)
	      normN = math.sqrt(tmp[0] ** 2 + tmp[1] ** 2 + tmp[2] ** 2)	      
	      tmp[0] = tmp[0]/normN
	      tmp[1] = tmp[1]/normN
	      tmp[2] = tmp[2]/normN    
	      Norm_Mean[0]=Norm_Mean[0]+tmp[0]
	      Norm_Mean[1]=Norm_Mean[1]+tmp[1]
	      Norm_Mean[2]=Norm_Mean[2]+tmp[2]
	      total_nodes=total_nodes+1	  
	  Norm_Mean[0]=Norm_Mean[0]/total_nodes
	  Norm_Mean[1]=Norm_Mean[1]/total_nodes
	  Norm_Mean[2]=Norm_Mean[2]/total_nodes	  
	  normN_promedio = math.sqrt(Norm_Mean[0] ** 2 + Norm_Mean[1] ** 2 + Norm_Mean[2] ** 2)
	  Norm_Mean[0] =  Norm_Mean[0]/normN_promedio
	  Norm_Mean[1] =  Norm_Mean[1]/normN_promedio
	  Norm_Mean[2] =  Norm_Mean[2]/normN_promedio
	  #print(Norm_Mean)
	  for node in inlet_nodes_3d:
	      n = node.GetSolutionStepValue(NORMAL)
	      #a = math.sqrt(n[0] ** 2 + n[1] ** 2 + n[2] ** 2) Direct Method, up-stimation
	      a = math.sqrt((n[0]*Norm_Mean[0]) ** 2 + (n[1]*Norm_Mean[1]) ** 2 + (n[2]*Norm_Mean[2]) ** 2)
	      area3d += a
	      
	  self.outlet_areas_3d.append(area3d)
	  Total_area_outlet = Total_area_outlet + area3d
	  print("area3d_outlet", Total_area_outlet)
      Total_area_outlet = Total_area_outlet / len(self.inlets_3d)


    def CheckPressure(self,model_part_3d):
	for i in range(0,len(self.skin_3d)):
		skin_node = self.skin_3d[i]		
		for node in skin_node:
		  pressure_skin = node.GetSolutionStepValue(PRESSURE)
		  if (pressure_skin > (5*(config.systolic_pressure))):
			  max_pressure = 5*config.systolic_pressure
			  print ("Err: Please check your model. Pressure is higher than " ,str(max_pressure), " Pa. In node: ",str(node.Id), " pressure is: ", str(pressure_skin))
			  sys.exit("Proccess Kill")
			  break
			  

    def CheckConvergencePressure(self,model_part_3d,factor,ffit_3d,total_time):
	total_tol=0.0
	tol=0.0
	tol2=0.0
	print("Cheking Convergence")
	for i in range(0,len(self.skin_3d)):
		skin_node = self.skin_3d[i]		
		for node in skin_node:
			pressure_skin_current = node.GetSolutionStepValue(PRESSURE)
			pressure_skin_past=node.GetSolutionStepValue(PRESSURE,1)
			tol=tol+(pressure_skin_current-pressure_skin_past)**2
			tol2=tol2+(pressure_skin_past)**2

	total_tol=abs(math.sqrt(tol/tol2))
	if(total_tol<0.001):
		print ("Converge", str(total_tol))
		factor=factor+1
		print(factor)
		self.FitValues_3d(model_part_3d,ffit_3d,total_time)
	return[factor]
 
    def FitValues_3d(self,model_part_3d,ffit_3d,total_time):
		avg_press_in=0.0
		avg_press_out=0.0
		counter=0 
		press_3d=0.0
		for node in inlet_nodes_3d:
			press_3d += node.GetSolutionStepValue(PRESSURE)
			counter += 1.0
		avg_press_in = press_3d / counter		
		
		counter=0
		press_3d=0.0
		for node in outlet_nodes_3d:
			press_3d += node.GetSolutionStepValue(PRESSURE)
			counter += 1.0
		avg_press_out = press_3d / counter	
		
		dpi=0.0
		dpi=avg_press_in-avg_press_out
		
		flow_3d = 0.0
		for node in inlet_nodes_3d:
			normal = node.GetSolutionStepValue(NORMAL)
			vel = node.GetSolutionStepValue(VELOCITY)
			flow_3d += normal[0] * vel[0] + normal[1] * vel[1] + normal[2] * vel[2]
	
		#ToWrite = str(node_id) + "	 "
		ToWrite = str(total_time) + "	 "
		ToWrite += str(avg_press_in) + "		"
		ToWrite += str(avg_press_out) + "		"
		ToWrite += str(dpi) + "		"
		ToWrite += str(flow_3d) + "\n"
		ffit_3d[k].write(ToWrite)
			
	
    def InletConditions(self,model_part_3d,velocity,press):
	# Write Conditions Used
	# ARCHIVE TO SET :::::::::::::::::::::::::::>>>>>>>>>>>>>> VARIABLES
	# import config_full
	#print("Dyastolic_Pressure", dyastolic_pressure)
	print("Setting InletConditions")
	#print(self.Velocity_Sig)
	#print(velocity)		
	#raw_input()
	for i in range(0, len(self.inlets_3d)):
	  inlet_nodes_3d = self.inlets_3d[i]
	  #area3d = self.inlet_areas_3d[i]
	  directions = self.inlet_velocity_directions[i]
	  #radio3d = math.sqrt(area3d * 3.1416)    			    
	  k = 0
	  for node in inlet_nodes_3d: 
		orientation = directions[k]
		node.SetSolutionStepValue(VELOCITY, 0, directions[k] * (self.Velocity_Sig) *velocity)
		#print("setting velocity, in ", node.Id, "---> ", directions[k] * (self.Velocity_Sig) *velocity)
		k = k + 1
		
	#raw_input()
		
	print("OutletConditions")
	for i in range(0, len(self.outlets_3d)):
		outlet_nodes_3d=self.outlets_3d[i]
		for node in outlet_nodes_3d:
			node.SetSolutionStepValue(PRESSURE, 0, press)
	
	#for i in range(0,len(self.skin_3d)):
	#	skin_node = self.skin_3d[i]		
	#	for node in skin_node:
	#		node.SetSolutionStepValue(VELOCITY_X, 0, 0.0)
	#		node.SetSolutionStepValue(VELOCITY_Y, 0, 0.0)
	#		node.SetSolutionStepValue(VELOCITY_Z, 0, 0.0)	
		
	def OutletConditions(self,model_part_3d):
		press=0.0
		print("utletConditions")
		for i in range(0, len(self.outlets_3d)):  
			for node in outlet_nodes_3d:
				node.SetSolutionStepValue(PRESSURE, 0, press)
    



#-------------------------------------------------------------------------
# Setting Contitions 3d
#-------------------------------------------------------------------------
    def Setting3d(self, initial_pressure, blood_density, blood_viscosity):
	for node in self.model_part_3d.Nodes:
			node.Free(VELOCITY_X)
			node.Free(VELOCITY_Y)
			node.Free(VELOCITY_Z)
			node.Free(PRESSURE)
			node.SetSolutionStepValue(VISCOSITY, 0, blood_viscosity)
			node.SetSolutionStepValue(DENSITY, 0, blood_density)
			node.SetSolutionStepValue(PRESSURE, 0, initial_pressure)

	
	for cond in self.model_part_3d.Conditions:
			if(cond.Properties.Id == 1):  # sides --> note that this is done in an outer separated loop!! SURFACE
			    for node in cond.GetNodes():
						  node.Fix(VELOCITY_X)
						  node.Fix(VELOCITY_Y)
						  node.Fix(VELOCITY_Z)
						  node.SetSolutionStepValue(FLAG_VARIABLE, 0, 0.0)
	
	for cond in self.model_part_3d.Conditions:
			if(cond.Properties.Id > 1000):  # outlet
			    for node in cond.GetNodes():
						  # if(not node.IsFixed(VELOCITY_X)):
						  node.Fix(PRESSURE)
						  node.SetSolutionStepValue(FLAG_VARIABLE, 0, cond.Properties.Id)
	
	for cond in self.model_part_3d.Conditions:
		  if(cond.Properties.Id == 100):  # inlet
		      for node in cond.GetNodes():
						  node.Fix(VELOCITY_X)
						  node.Fix(VELOCITY_Y)
						  node.Fix(VELOCITY_Z)
						  node.SetSolutionStepValue(FLAG_VARIABLE, 0, 100.0)
	
	for node in self.model_part_3d.Nodes:
			y = node.GetSolutionStepValue(Y_WALL, 0)
			node.SetValue(Y_WALL, y)
			
	
	counter = 0.0
	for node in self.model_part_3d.Nodes:
			if(node.IsFixed(PRESSURE)):
	  			counter += 1.0
   

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

    def __init__(self, model_part_1d, model_part_3d):
        self.model_part_1d = model_part_1d
        self.model_part_3d = model_part_3d

    def Initialize(self):
      self.inlets_1d = []
      self.outlets_1d = []
      self.inlets_3d = []
      self.outlets_3d = []
      self.skin_3d = []
      self.inlet_areas_3d = []
      self.outlet_areas_3d = []
      self.inlet_velocity_directions = []
      self.flow_1d = 0.0
      self.flow_1d_pressure = 0.0
      self.flow_3d_in = 0.0
      self.flow_3d_in_pressure = 0.0
      self.flow_3d_out = 0.0
      self.flow_1d_out = 0.0
      self.flow_1d_in2 = 0.0
      self.Aprime = []
      self.Bprime = []
      # compute normals and 3d areas
      Coupling_version="CoupledTool1D_3Dv5.py:03_06_2014"
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
      
      # detect 1d inlets	        
      for i in range(100, 101):
	  nfound = 0
	  aux = []
	  for node in self.model_part_1d.Nodes:
	      if (node.GetSolutionStepValue(FLAG_VARIABLE) == i):
		  aux.append(node)
		  # node.Fix(NODAL_AREA)
		  # node.Fix(FLOW)
		  # print "1D_inlet has been assigned", node
		  # print self.model_part_1d[node].GetValue(RADIUS)
		  # inlet_nodes_1d[0].GetSolutionStepValue(FLOW)
		  # raw_input()
		  #print (node)
		  Node_reference = node
	  if(len(aux) != 0):
	      self.inlets_1d.append(aux)
	  else:
	      break
      
      if(len(self.inlets_1d) == 0):
	  # print "number of 1d is zero"
	  sys.exit("number of 1d_inlets are zero!! Please check your config.py file!")

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
		#print(Norm_Mean)
		#for node in self.model_part_3d.Nodes:
		      #tmp = 0.0
		      #if (node.GetSolutionStepValue(FLAG_VARIABLE) == i):
			      #aux.append(node)
			      #node.Fix(VELOCITY_X)
			      #node.Fix(VELOCITY_Y)
			      #node.Fix(VELOCITY_Z)
			      #tmp = node.GetSolutionStepValue(NORMAL)
			      #normN = math.sqrt(tmp[0] ** 2 + tmp[1] ** 2 + tmp[2] ** 2)
			      ##print (str(tmp))
			      ##print (str(normN))
			      ##tmp = tmp/normN
			      #tmp[0] =  tmp[0]/normN_promedio
			      #tmp[1] =  tmp[1]/normN_promedio
			      #tmp[2] =  tmp[2]/normN_promedio	                 
			      #directions.append(tmp)
			      #print(tmp)
		#print(".....")	      
		#normN=0.0
		#normN_promedio = 0.0
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

      if(len(self.inlets_1d) != len(self.inlets_3d)):
	  sys.exit("number of 1d and 3d inlets is different!! Please check your config.py file!")
	  # print "number of 1d and 3d inlets is different"
	  # err

        # detect 1d outlets

      for i in range(1001, 1199):
	  # print i
	  nfound = 0
	  aux = []
	  for node in self.model_part_1d.Nodes:
	      if (node.GetSolutionStepValue(FLAG_VARIABLE) == i):
		  aux.append(node)
		  # print "1D_outlet has been assigned(H)", node
		  # node.Fix(FLOW)
	  if(len(aux) != 0):
	      self.outlets_1d.append(aux)
	      # raw_input()
	  else:
	      break

      if(len(self.outlets_1d) == 0):
	  # print "number of 1d is zero"
	  sys.exit("number of 1d_outlets are zero!! Please check your config.py file!")

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
	  print(Norm_Mean)
	  for node in inlet_nodes_3d:
	      n = node.GetSolutionStepValue(NORMAL)
	      #a = math.sqrt(n[0] ** 2 + n[1] ** 2 + n[2] ** 2) Direct Method, up-stimation
	      a = math.sqrt((n[0]*Norm_Mean[0]) ** 2 + (n[1]*Norm_Mean[1]) ** 2 + (n[2]*Norm_Mean[2]) ** 2)
	      area3d += a
	      
	  self.outlet_areas_3d.append(area3d)
	  Total_area_outlet = Total_area_outlet + area3d
	  print("area3d_outlet", Total_area_outlet)
      Total_area_outlet = Total_area_outlet / len(self.inlets_3d)
    
      for prop in self.model_part_1d.Properties:
	  A_prime = 0.0
	  B_prime = 0.0
	  if (prop.Id == config.deactivate_list[0]):
	      Reference_Radius = prop.GetValue(RADIUS)
	      Reference_viscosity = simulation_config.blood_static_viscosity
	      Reference_density = simulation_config.blood_density
	      Area_Stenosis = 0.5 * Total_area_outlet
	      A_prime = (8 * Reference_viscosity) / (3.1416 * math.pow(Reference_Radius, 4))
	      B_prime = (Reference_density) * ((Total_area_outlet / Area_Stenosis) - 1) * ((Total_area_outlet/Area_Stenosis)-1)/(2*Total_area_outlet*Total_area_outlet)
	      self.Aprime.append(A_prime)
	      self.Bprime.append(B_prime)
	      # print "A_prime",A_prime
	      # print "B_prime",B_prime
	      break

      if (simulation_config.FitRadius):
	  meanArea = Total_area_inlet
	  meanRadius = math.sqrt(meanArea / math.pi)
	  RadiusFactor = meanRadius / Reference_Radius
	  print("1D Inlet Reference_Radius", Reference_Radius)
	  print("3D Inlet Reference Radius", meanRadius)
	  print("Radius Factor 3D/1D", RadiusFactor)
	  print("Setting 1D New Radius - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ")
	  for prop in self.model_part_1d.Properties:
		R = prop.GetValue(RADIUS)
		R_fit = R * RadiusFactor
		prop.SetValue(RADIUS, R_fit)
      
	  #print ("R", R)
	  #print ("R_FIT", R_fit)
	  #print ("RadiusFactor", RadiusFactor)
      # print "----------------------------------------------"
      #raw_input()

      import FitAB
      self.fitters_1d = []
      self.fitters_3d = []
      self.fitters_inlet_nodes = []
      self.fitters_outlet_nodes = []

      for i in range(0, len(self.outlets_1d)):
	  self.fitters_1d.append(FitAB.Fitter(self.outlets_1d[i][0].Id))
	  self.fitters_3d.append(FitAB.Fitter(self.outlets_1d[i][0].Id))

      self.fixed_flow_nodes = []
      for node in self.model_part_1d.Nodes:
	  if(node.IsFixed(FLOW) and node.GetSolutionStepValue(FLAG_VARIABLE) == 0):
	      self.fixed_flow_nodes.append(node)
      for i in range(0, len(self.fixed_flow_nodes)):
	  self.fitters_inlet_nodes.append(FitAB.Fitter(self.fixed_flow_nodes[i].Id))
	  #print(self.fixed_flow_nodes[i].Id)

      self.fixed_outlet_nodes = []
      for node in self.model_part_1d.Nodes:
	  if(node.IsFixed(PRESSURE)):
	      self.fixed_outlet_nodes.append(node)
      for i in range(0, len(self.fixed_outlet_nodes)):
	  self.fitters_outlet_nodes.append(FitAB.Fitter(self.fixed_outlet_nodes[i].Id))
	  #print(self.fixed_outlet_nodes[i].Id)

    def Velocity_Initial_Contitions(self):
        print("Initialize Velocity 3D contidion")
        print(self.Velocity_Sig) 
        #raw_input()
        #raw_input()
        for i in range(0, len(self.inlets_1d)):
            inlet_nodes_1d = self.inlets_1d[i]
            # print "3D-1D: inlet_nodes_1d [0].Id::::::>>>> ", inlet_nodes_1d[0].Id
            for i in range(0, len(self.outlets_1d)):
                outlet_nodes_1d = self.outlets_1d[i]
                pressinlet3D = outlet_nodes_1d[0].GetSolutionStepValue(PRESSURE)
                # print pressinlet3D
                # print "3D-1D: outlet_nodes_1d[0].Id::::::>>>> ", outlet_nodes_1d[0].Id
                for i in range(0, len(self.inlets_3d)):
                    inlet_nodes_3d = self.inlets_3d[i]
                    area3d = self.inlet_areas_3d[i]
                    print("area_inlet_3d", area3d)
                    directions = self.inlet_velocity_directions[i]
                    radio3d = math.sqrt(area3d * 3.1416)
                    vel1d = inlet_nodes_1d[0].GetSolutionStepValue(FLOW) / area3d
                    k = 0
                    #print("velocity 1D: ", str(vel1d))
                    # Impongo velocidad y la presion
                    for node in inlet_nodes_3d:
                        #n = node.GetSolutionStepValue(NORMAL)
                        #a = math.sqrt(n[0] ** 2 + n[1] ** 2 + n[2] ** 2)
                        #orientation = directions[k]
                        node.SetSolutionStepValue(VELOCITY, 0, directions[k] * (self.Velocity_Sig) *vel1d)
                        #print("Velocity", str(directions[k] * vel1d))
                        #print("k", k)
                        # node.SetSolutionStepValue(PRESSURE, 0, pressinlet3D)
                        k = k + 1
                        #raw_input()
	print ("final to set 3D inlet velocity Conditions")
	#raw_input()

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
	#for node in self.model_part_3d.Nodes:
		#print (str(node.Id))
		#print (str(node.GetSolutionStepValue(VELOCITY)))
		#print (str(node.GetSolutionStepValue(VELOCITY_X)))
		#print (str(node.GetSolutionStepValue(VELOCITY_Y)))
		#print (str(node.GetSolutionStepValue(VELOCITY_Z)))
		

    def Transfer1D_to_3D(self, dyastolic_pressure,summary_file_pressure):
	# Write Conditions Used
	# ARCHIVE TO SET :::::::::::::::::::::::::::>>>>>>>>>>>>>> VARIABLES
	# import config_full
	#print("Dyastolic_Pressure", dyastolic_pressure)
	print("Transfer Conditions")
	print(self.Velocity_Sig)	
	#raw_input()
	for i in range(0, len(self.inlets_1d)):
		    inlet_nodes_1d = self.inlets_1d[i]
		    #print ("NODO 1D-3D:: inlet_nodes_1d[0].Id::::::inlet 1D coupled with the 3D inlet>>>> ",inlet_nodes_1d[0].Id)
		    # print "--"
		    for i in range(0, len(self.inlets_3d)):
			    #print (i)
			    inlet_nodes_3d = self.inlets_3d[i]
			    area3d = self.inlet_areas_3d[i]
			    directions = self.inlet_velocity_directions[i]
			    radio3d = math.sqrt(area3d * 3.1416)
			    vel1d = inlet_nodes_1d[0].GetSolutionStepValue(FLOW) / area3d			    			    
			    #print ("Estoy fijando la velocidad en el outlet del 3D que proviene del nodo ")
			    #print (str(vel1d))
			    #vel1d = 0.01
			    # NOTA: REVISAR EL CAUDAL QUE ESTAMOS IMPONIENDO DE ENTRADA El hecho de considerar q/2 es debido a que el caudal que imponemos
			    # pessinlet1D =
			    # inlet_nodes_1d[0].GetSolutionStepValue(PRESSURE)
			    k = 0
			    # Impongo velocidad y la presion
			    for node in inlet_nodes_3d:
				    #n = node.GetSolutionStepValue(NORMAL)
				    #a = math.sqrt(n[0] ** 2 + n[1] ** 2 + n[2] ** 2)
				    orientation = directions[k]
				    #print (str(orientation))
				    #print (str(directions[k] * (self.Velocity_Sig) *vel1d))
				    node.SetSolutionStepValue(VELOCITY, 0, directions[k] * (self.Velocity_Sig) *vel1d)
				    #print("Velocity", str(directions[k] * vel1d))
				    k = k + 1			   
				    #raw_input()
		    # print "area3d",area3d ,"area1d",inlet_nodes_1d[0].GetSolutionStepValue(NODAL_AREA), "del nodo", inlet_nodes_1d[0].Id
		    # print "Estoy fijando la velocidad en el inlet del 3D que proviene del nodo " , inlet_nodes_1d[0].Id, " del 1D. La velocidad que estoy fijando es"
		    # print "vel1d",vel1d, "del nodo", inlet_nodes_1d[0].Id

	#raw_input()
        
        for i in range(0, len(self.outlets_1d)):
            outlet_nodes_1d = self.outlets_1d[i]
            # print outlet_nodes_1d[0].Id
            pressinlet3D = outlet_nodes_1d[0].GetSolutionStepValue(PRESSURE)
            # pressinlet3D = 0
            # print "NODO 1D-3D: outlet_nodes_1d[0].Id::::::outlet 1D coupled with the 3D outlet>>>> ",outlet_nodes_1d[0].Id
            # print "pressinlet3D ",pressinlet3D
            # raw_input()
            # double beta = E0*thickness0*1.77245385/(1.0-nu0*nu0);
            beta = outlet_nodes_1d[0].GetSolutionStepValue(BETA)
            # beta = ((outlet_nodes_1d[0].GetSolutionStepValue(YOUNG_MODULUS) * outlet_nodes_1d[0].GetSolutionStepValue(THICKNESS) * math.sqrt(math.pi)) / (
                # 1 - (outlet_nodes_1d[0].GetSolutionStepValue(POISSON_RATIO) * outlet_nodes_1d[0].GetSolutionStepValue(POISSON_RATIO))))
            A = outlet_nodes_1d[0].GetSolutionStepValue(NODAL_AREA)
            # print "AREA_1D", A
            A0 = outlet_nodes_1d[0].GetValue(NODAL_AREA)
            press = dyastolic_pressure + beta * \
                (math.sqrt(A) - math.sqrt(A0)) / \
                A0  # math.sqrt(A/A0)*beta - beta
            # press = dyastolic_pressure + beta * \
                #(math.sqrt(A) - math.sqrt(A0)) / \
                # A0  # math.sqrt(A/A0)*beta - beta
            #print("Press calculated", press)
            #print("Press solver", pressinlet3D)
            # for i in range(0, len(self.outlets_3d)):
            outlet_nodes_3d = self.outlets_3d[i]
            area3d = self.outlet_areas_3d[i]
            for node in outlet_nodes_3d:
                node.SetSolutionStepValue(PRESSURE, 0, press)
                # print " ", press, "en el nodo", node.Id
            #print ("Area3d",area3d ,"Area1d",outlet_nodes_1d[0].GetSolutionStepValue(NODAL_AREA))
            #print ("Estoy fijando la presion en el outlet del 3D que proviene del nodo " , outlet_nodes_1d[i].Id, " del 1D. La presion que estoy fijando es")
            print ("Pressure::",press, "  in  Node::", outlet_nodes_1d[0].Id)
	    ToWriteIn_Summary= "Instant pressure in the outlet node:  "  +  str(press) + " Pa " + "\n"       
	    summary_file_pressure.write(ToWriteIn_Summary)

	
    def Transfer3D_to_1D(self, dyastolic_pressure):
        inlet_flow = (self.inlets_1d[0])[0].GetSolutionStepValue(FLOW)
        self.outlets_1d[0][0].SetSolutionStepValue(FLOW, 0, inlet_flow)

        # outlet_area = self.outlets_3d[0][0].GetSolutionStepValue(NODAL_AREA)
        # outlet_beta = self.outlets_3d[0][0].GetSolutionStepValue(BETA)
        outlet_pressure = self.outlets_1d[0][0].GetSolutionStepValue(PRESSURE)

        inlet_beta = self.inlets_1d[0][0].GetSolutionStepValue(BETA)
        inlet_A0 = self.inlets_1d[0][0].GetValue(NODAL_AREA)
        Ainlet_to_prescribe = ((((outlet_pressure - dyastolic_pressure) * inlet_A0) / inlet_beta) + math.sqrt(inlet_A0)) ** 2  # A0*(avg_press/beta + 1)**2
        # self.outlets_1d[0][0].SetSolutionStepValue(NODAL_AREA,0,Ainlet_to_prescribe)

        # print "inlet node = ",self.inlets_1d[0][0].Id
        # print "outlet node = ",self.outlets_1d[0][0].Id


    # def Transfer3D_to_1D(self):
        # for i in range(0, len(self.inlets_3d)):
            # inlet_nodes_1d = self.inlets_1d[i]
            # inlet_nodes_3d = self.inlets_3d[i]
            # area3d = self.inlet_areas_3d[i]
            # outlet_nodes_1d = self.outlets_1d[i]
            # press_3d = 0.0
            # counter = 0.0
            # for node in inlet_nodes_3d:
                # press_3d += node.GetSolutionStepValue(PRESSURE)
                # counter += 1.0
            # avg_press = press_3d / counter
            # beta = inlet_nodes_1d[0].GetSolutionStepValue(YOUNG_MODULUS) * inlet_nodes_1d[
                # 0].GetSolutionStepValue(THICKNESS) * math.sqrt(math.pi)
            # beta = beta / \
                # (1.0 - (inlet_nodes_1d[0].GetSolutionStepValue(
                    # POISSON_RATIO) * inlet_nodes_1d[0].GetSolutionStepValue(POISSON_RATIO)))
            # print "beta", beta
            # beta = inlet_nodes_1d[0].GetSolutionStepValue(BETA)
            # print "beta", beta
            # print "inlet node is ", inlet_nodes_1d[0].Id
            # beta = inlet_nodes_1d[0].GetSolutionStepValue(BETA)
            # A0 = inlet_nodes_1d[0].GetValue(NODAL_AREA)
            # Ainlet_to_prescribe = (((avg_press * A0)/beta) + math.sqrt(A0))**2  # A0*(avg_press/beta + 1)**2
            # print "in Transfer AREA 3D_to_1D----->", Ainlet_to_prescribe, " to node" , inlet_nodes_1d[0].Id
            # print "Area del 1D Original", A0
            # print "Area del 1D Modificada", inlet_nodes_1d[0].GetSolutionStepValue(NODAL_AREA)


            # A0 = outlet_nodes_1d[0].GetValue(NODAL_AREA)
            # outlet_nodes_1d[0].SetSolutionStepValue(NODAL_AREA ,0, A0)
            # inlet_nodes_1d[0].SetSolutionStepValue(NODAL_AREA, 0, Ainlet_to_prescribe)

            # compute flow on inlet
            # flow = 0.0
            # for node in inlet_nodes_3d:
                # normal = node.GetSolutionStepValue(NORMAL)
                # vel = node.GetSolutionStepValue(VELOCITY)
                # flow += normal[0] * vel[0] + normal[
                    # 1] * vel[1] + normal[2] * vel[2]
            # self.flow_3d_in = str(flow)
            # print "TRASFER to 3D to 1D:::flow entering= ",flow
            # print "Flow de 1D inlet",  inlet_nodes_1d[0].GetSolutionStepValue(FLOW)
            # flow_aux=  flow

        # for i in range(0, len(self.outlets_3d)):
            # outlet_nodes_1d = self.outlets_1d[i]
            # outlet_nodes_3d = self.outlets_3d[i]
            # area3d = self.outlet_areas_3d[i]

            # compute flow on 3D outlet
            # flow = 0.0
            # for node in outlet_nodes_3d:
                # normal = node.GetSolutionStepValue(NORMAL)
                # vel = node.GetSolutionStepValue(VELOCITY)
                # flow += normal[0] * vel[0] + normal[
                    # 1] * vel[1] + normal[2] * vel[2]
            # assign flow to outlet
            # outlet_nodes_1d[0].SetSolutionStepValue(FLOW, 0, -flow)
            # print " velocity del 3D", vel
            # print " velocity del 1D", vel
            # print "in Transfer Flow 3D_to_1D----->", -flow, " to node" , outlet_nodes_1d[0].Id
        # raw_input()
#-------------------------------------------------------------------------
# FIT VALUES (WRITE FILE)
#-------------------------------------------------------------------------
#-------------------------------------------------------------------------
# Setting Contitions 3d
#-------------------------------------------------------------------------
    def Setting3d(self, initial_pressure, blood_density, blood_viscosity):
        # print self
            # print len(self.model_part_3d.Conditions)
      a=0
      for node in self.model_part_3d.Nodes:
	  node.Free(VELOCITY_X)
	  node.Free(VELOCITY_Y)
	  node.Free(VELOCITY_Z)
	  node.Free(PRESSURE)
	  node.SetSolutionStepValue(VISCOSITY, 0, blood_viscosity)
	  node.SetSolutionStepValue(DENSITY, 0, blood_density)
	  node.SetSolutionStepValue(PRESSURE, 0, initial_pressure)
# set inlet
     
      for cond in self.model_part_3d.Conditions:
	  if(cond.Properties.Id == 1):  # sides --> note that this is done in an outer separated loop!! SURFACE
	      for node in cond.GetNodes():
		  node.Fix(VELOCITY_X)
		  node.Fix(VELOCITY_Y)
		  node.Fix(VELOCITY_Z)
		  node.SetSolutionStepValue(FLAG_VARIABLE, 0, 0.0)
		  #print ("2")
		  #print (cond)
		  #print (cond.Properties.Id)
		  #a=a+1

	# set output (overwrites the others)
      #print(a)
      #a=0
      #raw_input()
      for cond in self.model_part_3d.Conditions:
	  if(cond.Properties.Id > 1000):  # outlet
	      for node in cond.GetNodes():
		  # if(not node.IsFixed(VELOCITY_X)):
		  node.Fix(PRESSURE)
		  node.SetSolutionStepValue(FLAG_VARIABLE, 0, cond.Properties.Id)
		  #print ("3")
		  #print (cond)
		  #print (cond.Properties.Id)
		  #a=a+1
      # for cond in self.model_part_3d.Conditions:
	    # if(cond.GetValue(IS_STRUCTURE) == 1):
		# for node in cond.GetNodes():
		    # node.Fix(VELOCITY_X)
		    # node.Fix(VELOCITY_Y)
		    # node.Fix(VELOCITY_Z)
		    # node.SetSolutionStepValue(VELOCITY_X,0,0.0)
	# this is only in this example...it should be set up by the problemtype
	# for node in self.model_part_3d.Nodes:
	    # node.SetSolutionStepValue(VISCOSITY,0,0.0035/1060.0)
	    # node.SetSolutionStepValue(DENSITY,0,1060.0)
	    # if(node.IsFixed(VELOCITY_X) == True and
	    # node.GetSolutionStepValue(VELOCITY_X) > 0.0001):
		# node.SetSolutionStepValue(FLAG_VARIABLE,0,100.0)
	    # if(node.IsFixed(PRESSURE) == True):
		# node.SetSolutionStepValue(FLAG_VARIABLE,0,101.0)
	    # node.SetSolutionStepValue(VELOCITY_X,0,0.0)
	    # copy Y_WALL
      #print(a)
      #a=0
      #raw_input()
      
      for cond in self.model_part_3d.Conditions:
	  if(cond.Properties.Id == 100):  # inlet
	      for node in cond.GetNodes():
		  node.Fix(VELOCITY_X)
		  node.Fix(VELOCITY_Y)
		  node.Fix(VELOCITY_Z)
		  node.SetSolutionStepValue(FLAG_VARIABLE, 0, 100.0)
		  #print (node.Id)
		  #print (cond.Properties.Id)
		  #a=a+1
	  # if(cond.Properties.Id > 1000): ##outlet

	      # for node in cond.GetNodes():

		  # node.Fix(PRESSURE)

		  # node.SetSolutionStepValue(FLAG_VARIABLE,0,cond.Properties.Id)

	  # if(cond.Properties.Id > 1001): ##outlet

	      # for node in cond.GetNodes():

		  # node.Fix(PRESSURE)

		  # node.SetSolutionStepValue(FLAG_VARIABLE,0,1002.0)

	# set sides (overwrites the inlet)
      #print(a)
      #a=0
      #raw_input()
      
      if (simulation_config.Use_Catheter):
	    for node in self.model_part_3d.Nodes:
		y = node.GetSolutionStepValue(Y_WALL, 0)
		node.SetValue(Y_WALL, y)

      counter = 0.0
      for node in self.model_part_3d.Nodes:
	  if(node.IsFixed(PRESSURE)):
	      counter += 1.0
	      #print ("4")
	      #print (node)

      #print (counter)
      #raw_input()      
        # print "n pressure nodes ",counter
    # def FitValues_Inlet(self):
      # for i in range(0,len(self.fixed_flow_nodes)):
        # node=self.fixed_flow_nodes[i]
        # for fitter in self.fitters_inlet_nodes:
            # fitter.AddPin(node[0].GetSolutionStepValue(PRESSURE))
            # fitter.AddQ(node[0].GetSolutionStepValue(FLOW))
      # for i in range(0,len(self.fixed_outlet_nodes)):
        # node=self.fixed_outlet_nodes[i]
        # for fitter in self.fitters_inlet_nodes:
            # fitter.AddPout(node[0].GetSolutionStepValue(PRESSURE))
    # raw_input()
    
    def FitValues_1d(self, total_time):
        for i in range(0, len(self.inlets_1d)):
            node = self.inlets_1d[i]
            # print "3D-1D: inlet_nodes_1d [0].Id::::::>>>> ", node[0].Id
            # nodewrite = node[0].Id
            # ToWrite = str(indextowrite) + " "
            # ToWrite = str(nodewrite) + " " + str(total_time) + " " + str(node[0].GetSolutionStepValue(PRESSURE)) + " "
            # ToWrite += str(node[0].GetSolutionStepValue(FLOW)) + " " + str(node[0].GetSolutionStepValue(NODAL_AREA)) + "\n"
            # ToWrite = str(node[0].GetSolutionStepValue(FLOW)) + "\n"
            # ffit[0].write(ToWrite)
            # AB_Matriz[row_AB][0]= total_time
            # AB_Matriz[row_AB][1]= node[0].GetSolutionStepValue(PRESSURE)
            # AB_Matriz[row_AB][2]= node[0].GetSolutionStepValue(FLOW)
            # AB_Matriz[row_AB][3]= node[0].GetSolutionStepValue(NODAL_AREA)
            for fitter in self.fitters_1d:
                fitter.AddPin(node[0].GetSolutionStepValue(PRESSURE))
                fitter.AddTotal_time(total_time)
            j = 0
            for j in range(0, len(self.outlets_1d)):
                # j=j+1
                # row_AB=row_AB+1
                node = self.outlets_1d[j]
                # print "3D-1D: outlet_nodes_1d[0].Id::::::>>>> ", node[0].Id
                # nodewrite = node[0].Id
                # ToWrite = str(indextowrite) + " "
                # ToWrite = str(nodewrite) + " " + str(total_time) + " " + str(node[0].GetSolutionStepValue(PRESSURE)) + " "
                # ToWrite += str(node[0].GetSolutionStepValue(FLOW)) + " " + str(node[0].GetSolutionStepValue(NODAL_AREA)) + "\n"
                # ToWrite = str(node[0].GetSolutionStepValue(FLOW)) + "\n"
                # ffit[j].write(ToWrite)
                # AB_Matriz[row_AB][0]= total_time
                # AB_Matriz[row_AB][1]= node[0].GetSolutionStepValue(PRESSURE)
                # AB_Matriz[row_AB][2]= node[0].GetSolutionStepValue(FLOW)
                # AB_Matriz[row_AB][3]= node[0].GetSolutionStepValue(NODAL_AREA)
                self.fitters_1d[j].AddPout(node[0].GetSolutionStepValue(PRESSURE))
                self.fitters_1d[j].AddQ(node[0].GetSolutionStepValue(FLOW))

    def FitValues_3d(self, total_time, ffit_test, kkkkkk):
        for i in range(0, len(self.inlets_1d)):
            inlet_nodes_1d = self.inlets_1d[i]
            inlet_nodes_3d = self.inlets_3d[i]
            press_3d = 0.0
            area3d = 0.0
            counter = 0.0
            # print "inlet", inlet_nodes_1d[0]
            for node in inlet_nodes_3d:
                press_3d += node.GetSolutionStepValue(PRESSURE)
                counter += 1.0
            avg_press = press_3d / counter
            # area3d=self.inlet_areas_3d
            for fitter in self.fitters_3d:
                fitter.AddPin(avg_press)
                fitter.AddTotal_time(total_time)
                # self.fitters_3d[i].AddPin(avg_press)
            # print "pressure", avg_press

            # print "HOLLLLLLLLLLLLLLLLLLAAAAAAAAAAAAAAAAAAAAAAAA"
        j = 0  # print inlet_nodes_1d[0]
        for j in range(0, len(self.outlets_1d)):
            outlet_nodes_1d = self.outlets_1d[j]
            outlet_nodes_3d = self.outlets_3d[j]
            area3d = self.outlet_areas_3d[j]
            # compute flow on 3D outlet
            flow_3d = 0.0
            for node in outlet_nodes_3d:
                normal = node.GetSolutionStepValue(NORMAL)
                vel = node.GetSolutionStepValue(VELOCITY)
                flow_3d += normal[0] * vel[0] + normal[1] * vel[1] + normal[2] * vel[2]
            # print "node", node
            # print "velo", vel
            # print "pressure1D",outlet_nodes_1d[0].GetSolutionStepValue(PRESSURE)
            # print "flow", outlet_nodes_1d[0].GetSolutionStepValue(FLOW)
            print("flow_3d", flow_3d)
            print("node", outlet_nodes_1d[0])
            print("flow_1d", outlet_nodes_1d[0].GetSolutionStepValue(FLOW))
            self.fitters_3d[j].AddQ(flow_3d)
            self.fitters_3d[j].AddPout(outlet_nodes_1d[0].GetSolutionStepValue(PRESSURE))
            # raw_input()
            # print "node outlet", outlet_nodes_1d[0]
            # print "flow", flow
            # print "pressure", outlet_nodes_1d[0].GetSolutionStepValue(PRESSURE)
            ToWrite = str(flow_3d) + "\n"
            ffit_test[kkkkkk].write(ToWrite)

            # for i in range(0, len(self.outlets_3d)):
            # outlet_nodes_1d = self.outlets_1d[i]
            # outlet_nodes_3d = self.outlets_3d[i]
            # area3d = self.outlet_areas_3d[i]
            # compute flow on 3D outlet
            # flow = 0.0
            # for node in outlet_nodes_3d:
                # normal = node.GetSolutionStepValue(NORMAL)
                # vel = node.GetSolutionStepValue(VELOCITY)
                # flow += normal[0] * vel[0] + normal[
                    # 1] * vel[1] + normal[2] * vel[2]
            # assign flow to outlet
            # outlet_nodes_1d[0].SetSolutionStepValue(FLOW, 0, -flow)
            # print " velocity del 3D", vel
            # print " velocity del 1D", vel
            # print "in Transfer Flow 3D_to_1D----->", -flow, " to node" , outlet_nodes_1d[0].Id
            
            
    def Save_Values_For_different_situations(self,total_flow_1d_value,total_in_pressure_1d_value,total_out_pressure_1d_value):	
	for i in range(0, len(self.inlets_1d)):
		inlet_nodes_1d = self.inlets_1d[i]		
		for i in range(0, len(self.inlets_3d)):
			#print (i)
			inlet_nodes_3d = self.inlets_3d[i]
			area3d = self.inlet_areas_3d[i]
			directions = self.inlet_velocity_directions[i]
			radio3d = math.sqrt(area3d * 3.1416)
			vel1d = inlet_nodes_1d[0].GetSolutionStepValue(FLOW) / area3d
			#self.fitters_1d.AddQ(inlet_nodes_1d[0].GetSolutionStepValue(FLOW))
			#print (inlet_nodes_1d[0].GetSolutionStepValue(FLOW))
			#print (inlet_nodes_1d[0].Id)
			total_flow_1d_value = total_flow_1d_value + inlet_nodes_1d[0].GetSolutionStepValue(FLOW)
			#print(total_flow_1d_value)
			total_in_pressure_1d_value=total_in_pressure_1d_value+inlet_nodes_1d[0].GetSolutionStepValue(PRESSURE)
			#print(total_in_pressure_1d_value)
	for j in range(0, len(self.outlets_1d)):
		outlet_nodes_1d = self.outlets_1d[j]
		total_out_pressure_1d_value=total_out_pressure_1d_value+outlet_nodes_1d[0].GetSolutionStepValue(PRESSURE)	
	return [total_flow_1d_value,total_in_pressure_1d_value,total_out_pressure_1d_value]
	
    
    def ComputeFFR_Health_Values(self,total_aortic_pressure,mean_pressure_1d_value,out_pressure_1d_value,mean_flow_1d_value,counter,summary_file):
	
	ToWriteIn_Summary= "FFR values in hypothic healthy situation" + "\n"
	ToWriteIn_Summary+= "total_aortic_pressure: " + str(total_aortic_pressure) + "\n"
	ToWriteIn_Summary+="mean_pressure_1d_value:  "+ str(mean_pressure_1d_value)	 + "\n"
	ToWriteIn_Summary+="out_pressure_1d_value "+ str(out_pressure_1d_value) + "\n"
	ToWriteIn_Summary+="mean_flow_1d_value "+ str(mean_flow_1d_value) + "\n"
	ToWriteIn_Summary+="-----------------------------------------------------------"+ "\n"
	summary_file.write(ToWriteIn_Summary)   
	
	Mean_Pressure_OUT_TOTAL=0.0
	Mean_Pressure_IN_TOTAL=0.0
	Mean_total_aortic_pressure=0.0	
	instant_pressure=0.0
	FFR_Value = 0.0
	
	mean_flow_1d_value = mean_flow_1d_value / counter	
	Mean_Pressure_IN_TOTAL = mean_pressure_1d_value / counter
	Mean_total_aortic_pressure=total_aortic_pressure/counter

	Mean_Pressure_OUT_TOTAL=out_pressure_1d_value/counter
	FFR_Value_Aortic=Mean_Pressure_OUT_TOTAL/Mean_total_aortic_pressure
	result_FFR_3D = Mean_Pressure_OUT_TOTAL/Mean_Pressure_IN_TOTAL
	#print(total_aortic_pressure)
	print("FFR values in hypothic healthy situation")
	print("Mean Pressure in the Aortic Root: ", Mean_total_aortic_pressure)
	print("Mean Pressure in the Inlet: ", Mean_Pressure_IN_TOTAL)	
	print("Mean Pressure in the Outlet: ", Mean_Pressure_OUT_TOTAL)
	print("FFR_Value taking the Aortic Root as reference:   ", FFR_Value_Aortic)
	print("FFR_Value taking the Inlet-node as reference:    ", result_FFR_3D)
	ToWriteIn_Summary= "FFR values in hypothic healthy situation" + "\n"
	ToWriteIn_Summary+= "Mean Pressure in the Aortic Root: " + str(Mean_total_aortic_pressure) + "\n"
	ToWriteIn_Summary+="Mean Pressure in the Inlet:  "+ str(Mean_Pressure_IN_TOTAL)	 + "\n"
	ToWriteIn_Summary+="Mean Pressure in the Outlet:  "+ str(Mean_Pressure_OUT_TOTAL) + "\n"
	ToWriteIn_Summary+="FFR_Value taking the Aortic Root as reference:   "+ str(FFR_Value_Aortic) + "\n"
	ToWriteIn_Summary+="FFR_Value taking the Inlet-node as reference:    "+ str(result_FFR_3D) + "\n"
	ToWriteIn_Summary+="Step:    "+ str(counter) + "\n"
	ToWriteIn_Summary+="-----------------------------------------------------------"+ "\n"	
	summary_file.write(ToWriteIn_Summary)     

    def ComputeFFR_Rest_Values(self,total_aortic_pressure,mean_flow_1d_FFR,mean_pressure_1d_FFR,out_pressure_1d_value_FFR,counter,summary_file):
		
	ll
	ToWriteIn_Summary= "FFR values in hypothic healthy situation" + "\n"
	ToWriteIn_Summary+= "total_aortic_pressure: " + str(total_aortic_pressure) + "\n"
	ToWriteIn_Summary+="mean_flow_1d_FFR:  "+ str(mean_flow_1d_FFR)	 + "\n"
	ToWriteIn_Summary+="mean_pressure_1d_FFR "+ str(mean_pressure_1d_FFR) + "\n"
	ToWriteIn_Summary+="out_pressure_1d_value_FFR "+ str(out_pressure_1d_value_FFR) + "\n"
	ToWriteIn_Summary+="-----------------------------------------------------------"+ "\n"
      	summary_file.write(ToWriteIn_Summary)   
	
	Mean_Pressure_OUT_TOTAL=0.0
	Mean_Pressure_IN_TOTAL=0.0
	Mean_total_aortic_pressure=0.0
	instant_pressure=0.0
	FFR_Value = 0.0
	
	Mean_total_aortic_pressure=total_aortic_pressure/counter
	Mean_Pressure_IN_TOTAL = mean_pressure_1d_FFR / counter
	Mean_Pressure_OUT_TOTAL=out_pressure_1d_value_FFR/counter
	
	FFR_Value_Aortic=Mean_Pressure_OUT_TOTAL/Mean_total_aortic_pressure
	result_FFR_3D = Mean_Pressure_OUT_TOTAL/Mean_Pressure_IN_TOTAL
	#print(total_aortic_pressure)
	print("FFR values in Rest Situation")
	print("Mean Pressure in the Aortic Root: ", Mean_total_aortic_pressure)
	print("Mean Pressure in the Inlet: ", Mean_Pressure_IN_TOTAL)	
	print("Mean Pressure in the Outlet: ", Mean_Pressure_OUT_TOTAL)
	print("FFR_Value taking the Aortic Root as reference:   ", FFR_Value_Aortic)
	print("FFR_Value taking the Inlet-node as reference:    ", result_FFR_3D)
	ToWriteIn_Summary= "FFR values in Rest Situation " + "\n"
	ToWriteIn_Summary+= "Mean Pressure in the Aortic Root: " + str(Mean_total_aortic_pressure) + "\n"
	ToWriteIn_Summary+="Mean Pressure in the Inlet:  "+ str(Mean_Pressure_IN_TOTAL)	 + "\n"
	ToWriteIn_Summary+="Mean Pressure in the Outlet:  "+ str(Mean_Pressure_OUT_TOTAL) + "\n"
	ToWriteIn_Summary+="FFR_Value taking the Aortic Root as reference:   "+ str(FFR_Value_Aortic) + "\n"
	ToWriteIn_Summary+="FFR_Value taking the Inlet-node as reference:    "+ str(result_FFR_3D) + "\n"
	ToWriteIn_Summary+="Step:    "+ str(counter) + "\n"
	ToWriteIn_Summary+="-----------------------------------------------------------"+ "\n"
	summary_file.write(ToWriteIn_Summary)     
	
    def ComputeFFR_Hypermia_Values(self,total_aortic_pressure,mean_pressure_1d_value,out_pressure_1d_value,mean_flow_1d_value,counter,summary_file):
	
	uu
	ToWriteIn_Summary= "FFR values in hypothic healthy situation" + "\n"
	ToWriteIn_Summary+= "total_aortic_pressure: " + str(total_aortic_pressure) + "\n"
	ToWriteIn_Summary+="mean_pressure_1d_value:  "+ str(mean_pressure_1d_value)	 + "\n"
	ToWriteIn_Summary+="out_pressure_1d_value "+ str(out_pressure_1d_value) + "\n"
	ToWriteIn_Summary+="mean_flow_1d_value "+ str(mean_flow_1d_value) + "\n"
	ToWriteIn_Summary+="-----------------------------------------------------------"+ "\n"
	summary_file.write(ToWriteIn_Summary)   
	
	
	Mean_Pressure_OUT_TOTAL=0.0
	Mean_Pressure_IN_TOTAL=0.0
	Mean_total_aortic_pressure=0.0
	instant_pressure=0.0
	FFR_Value = 0.0	
	mean_flow_1d_value = mean_flow_1d_value / counter	
	Mean_Pressure_IN_TOTAL = mean_pressure_1d_value / counter
	Mean_total_aortic_pressure=total_aortic_pressure/counter

	Mean_Pressure_OUT_TOTAL=out_pressure_1d_value/counter
	FFR_Value_Aortic=Mean_Pressure_OUT_TOTAL/Mean_total_aortic_pressure
	result_FFR_3D = Mean_Pressure_OUT_TOTAL/Mean_Pressure_IN_TOTAL
	#print(total_aortic_pressure)
	print("FFR values in hypothic Hypermia situation")
	print("Mean Pressure in the Aortic Root: ", Mean_total_aortic_pressure)
	print("Mean Pressure in the Inlet: ", Mean_Pressure_IN_TOTAL)	
	print("Mean Pressure in the Outlet: ", Mean_Pressure_OUT_TOTAL)
	print("FFR_Value taking the Aortic Root as reference:   ", FFR_Value_Aortic)
	print("FFR_Value taking the Inlet-node as reference:    ", result_FFR_3D)
	ToWriteIn_Summary= "FFR values in hypothic Hypermia situation" + "\n"
	ToWriteIn_Summary+= "Mean Pressure in the Aortic Root: " + str(Mean_total_aortic_pressure) + "\n"
	ToWriteIn_Summary+="Mean Pressure in the Inlet:  "+ str(Mean_Pressure_IN_TOTAL)	 + "\n"
	ToWriteIn_Summary+="Mean Pressure in the Outlet:  "+ str(Mean_Pressure_OUT_TOTAL) + "\n"
	ToWriteIn_Summary+="FFR_Value taking the Aortic Root as reference:   "+ str(FFR_Value_Aortic) + "\n"
	ToWriteIn_Summary+="FFR_Value taking the Inlet-node as reference:    "+ str(result_FFR_3D) + "\n"
	ToWriteIn_Summary+="Step:    "+ str(counter) + "\n"
	ToWriteIn_Summary+="-----------------------------------------------------------"+ "\n"	
	summary_file.write(ToWriteIn_Summary)    	
	
    #def Save_Flow_1d_inlet(self, total_time,mean_flow_1d_inlet):
	#print(outlet_nodes_1d[0].Id)
	
	#for i in range(0, len(self.inlets_3d)):
		  #inlet_nodes_1d = self.inlets_1d[i]
		  #inlet_nodes_3d = self.inlets_3d[i]
		  #outlet_nodes_1d = self.outlets_1d[i]
	
	#for i in range(0, len(self.inlets_1d)):
		#node = self.inlets_1d[i]
		#print ("3D-1D: inlet_nodes_1d [0].Id::::::>>>> ", node[0].Id)
		##nodewrite = node[0].Id
		#j = 0
		#for j in range(0, len(self.outlets_1d)):
			#node = self.outlets_1d[j]
			#print ("3D-1D: outlet_nodes_1d[0].Id::::::>>>> ", node[0].Id)
			## nodewrite = node[0].Id                
	#raw_input()
    
    def Fit_ABValues_1D(self):
        ffit_1d = []
        return_list_1D = []
        k = 0
        for fitter in self.fitters_1d:
            ffit_1d.append(k)
            results = str("Fitter_1D_" + str(k) + "_.txt")
            ffit_1d[k] = open(results, 'w')
            ffit_1d[k].write("Node-Time-Pressure_Inlet_1D-Pressure_outlet_1d-PressureDrop-flow_1d) \n")
            [node_id, A, B] = fitter.DoFitting_1D(k, ffit_1d)
            return_list_1D.append([node_id, A, B])
            k = k + 1
            print("node_id = ", node_id)
            print("A = ", A)
            print("B = ", B)
        return return_list_1D

    def Fit_ABValues_3D(self):
        ffit_3d = []
        ffit_3d_2 = []
        ffit_3d_3 = []
        ffit_3d_4 = []
        ffit_3d_5 = []
        ffit_3d_6 = []
        return_list_3D = []
        k = 0
        for fitter in self.fitters_3d:
            print("----------------------------------------------------------------------------------------------------------------")
            ffit_3d_2.append(k)
            results = str("Fitter_3D_2" + str(k) + "_.txt")
            ffit_3d_2[k] = open(results, 'w')
            ffit_3d_2[k].write("Node	Time	Pressure_Inlet_3D	Pressure_outlet_3d	PressureDrop	flow_3d) \n")
            [node_id, A, B] = fitter.DoFitting_3D_2(k, ffit_3d_2)
            print("node_id = ", node_id)
            print("A = ", A)
            print("B = ", B)
            A = 0.0
            B = 0.0
            raw_input()
            print("----------------------------------------------------------------------------------------------------------------")
            # ffit_3d_3.append(k)
            # results = str("Fitter_3D_3" + str(k) + "_.txt")
            # ffit_3d_3[k] = open(results, 'w')
            # ffit_3d_3[k].write("Node	Time	Pressure_Inlet_3D	Pressure_outlet_3d	PressureDrop	flow_3d) \n")
            #[node_id,A,B] = fitter.DoFitting_3D_3(k,ffit_3d_3)
            # print "node_id = ",node_id
            # print "A = ",A
            # print "B = ",B
            # A = 0.0
            # B = 0.0
            print("----------------------------------------------------------------------------------------------------------------")
            # ffit_3d_4.append(k)
            # results = str("Fitter_3D_4" + str(k) + "_.txt")
            # ffit_3d_4[k] = open(results, 'w')
            # ffit_3d_4[k].write("Node	Time	Pressure_Inlet_3D	Pressure_outlet_3d	PressureDrop	flow_3d) \n")
            #[node_id,A,B] = fitter.DoFitting_3D_4(k,ffit_3d_4)
            # print "node_id = ",node_id
            # print "A = ",A
            # print "B = ",B
            # A = 0.0
            # B = 0.0
            print("--------------------------------------   A-B USED  ---------------------------------------------------------------")
            ffit_3d.append(k)
            results = str("Fitter_3D_TO_BE_USED_" + str(k) + "_.txt")
            ffit_3d[k] = open(results, 'w')
            ffit_3d[k].write("Node	Time	Pressure_Inlet_3D	Pressure_outlet_3d	PressureDrop	flow_3d) \n")
            [node_id, A, B] = fitter.DoFitting_3D(k, ffit_3d)
            print("node_id = ", node_id)
            print("A = ", A)
            print("B = ", B)
            return_list_3D.append([node_id, A, B])
            A = 0.0
            B = 0.0
	    raw_input()
            print("----------------------------------------------------------")
            # ffit_3d_5.append(k)
            # results = str("Fitter_3D_5" + str(k) + "_.txt")
            # ffit_3d_5[k] = open(results, 'w')
            # ffit_3d_5[k].write("Node	Time	Pressure_Inlet_3D	Pressure_outlet_3d	PressureDrop	flow_3d) \n")
            #[node_id,A,B] = fitter.DoFitting_3D_5(k,ffit_3d_5,self.Aprime[0],self.Bprime[0])
            # print "node_id = ",node_id
            # print "A' = ",self.Aprime[0]
            # print "B' = ",self.Bprime[0]
            # print "A = ",A
            # print "B = ",B
            # A = 0.0
            # B = 0.0
            print("----------------------------------------------------------")
            # ffit_3d_6.append(k)
            # results = str("Fitter_3D_6" + str(k) + "_.txt")
            # ffit_3d_6[k] = open(results, 'w')
            # ffit_3d_6[k].write("Node	Time	Pressure_Inlet_3D	Pressure_outlet_3d	PressureDrop	flow_3d) \n")
            #[node_id,x_new] = fitter.DoFitting_3D_Gauss(k,ffit_3d_6)
            # print "node_id = ",node_id
            # print "X = ", x_new
            k = k + 1
           # ...para los fitters_3d
        return return_list_3D

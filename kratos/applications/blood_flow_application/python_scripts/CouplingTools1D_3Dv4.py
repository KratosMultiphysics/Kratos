from KratosMultiphysics import *
from KratosMultiphysics.BloodFlowApplication import *
CheckForPreviousImport()

import math
import time
import sys
import config

#File: 27/10/2013

class TransferTools:

    def __init__(self, model_part_1d, model_part_3d):
        self.model_part_1d = model_part_1d
        self.model_part_3d = model_part_3d    
    # 3

    def Initialize(self):
        self.inlets_1d = []
        self.outlets_1d = []
        self.inlets_3d = []
        self.outlets_3d = []
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
        # compute normals and 3d areas
        print "CoupledTool:v4_27102003"
	#raw_input()
        NormalCalculationUtils().SwapNormals(self.model_part_3d)
        BodyNormalCalculationUtils().CalculateBodyNormals(
            self.model_part_3d, 3) 
	  
        # detect 1d inlets

        for i in range(100, 101):
            nfound = 0
            aux = []
            for node in self.model_part_1d.Nodes:
                if (node.GetSolutionStepValue(FLAG_VARIABLE) == i):
                    aux.append(node)
                    # node.Fix(NODAL_AREA)
                    # node.Fix(FLOW)
                    print "1D_inlet has been assigned", node
                    # raw_input()
            if(len(aux) != 0):
                self.inlets_1d.append(aux)
            else:
                break

        if(len(self.inlets_1d) == 0):
            # print "number of 1d is zero"
            sys.exit(
                "number of 1d_inlets are zero!! Please check your config.py file!")

        # detect 3d inlets
        for i in range(100, 101):
            nfound = 0
            aux = []
            directions = []
            for node in self.model_part_3d.Nodes:
                if (node.GetSolutionStepValue(FLAG_VARIABLE) == i):
                    aux.append(node)
                    node.Fix(VELOCITY_X)
                    node.Fix(VELOCITY_Y)
                    node.Fix(VELOCITY_Z)
                    tmp = node.GetSolutionStepValue(NORMAL)
                    normN = math.sqrt(tmp[0] ** 2 + tmp[1] ** 2 + tmp[2] ** 2)
                    tmp /= normN
                    directions.append(tmp)
                    print "3D_inlet has inlet_nodes_3d = self.inlets_3d[i]been assigned", node
            if(len(aux) != 0):
                self.inlets_3d.append(aux)
                self.inlet_velocity_directions.append(directions)
            else:
                break
            # raw_input()

        if(len(self.inlets_3d) == 0):
            sys.exit(
                "number of 3d_inlets are zero!! Please check your config.py file!")
            # print "number of 3d is zero"
            # err

        if(len(self.inlets_1d) != len(self.inlets_3d)):
            sys.exit(
                "number of 1d and 3d inlets is different!! Please check your config.py file!")
            # print "number of 1d and 3d inlets is different"
            # err

        # detect 1d outlets
        for i in range(1001, 1199):
            nfound = 0
            aux = []
            for node in self.model_part_1d.Nodes:
                if (node.GetSolutionStepValue(FLAG_VARIABLE) == i):
                    aux.append(node)
                    #node.Fix(FLOW)
            if(len(aux) != 0):
                self.outlets_1d.append(aux)
                print "1D_outlet has been assigned(H)", node
                # raw_input()
            else:
                break

        if(len(self.outlets_1d) == 0):
            # print "number of 1d is zero"
            sys.exit(
                "number of 1d_outlets are zero!! Please check your config.py file!")

        # detect 3d outlets
        for i in range(1001, 1199):
            nfound = 0
            aux = []
            for node in self.model_part_3d.Nodes:
                if (node.GetSolutionStepValue(FLAG_VARIABLE) == i):
                    aux.append(node)
                    node.Fix(PRESSURE)
            if(len(aux) != 0):
                self.outlets_3d.append(aux)
            else:
                break

        if(len(self.outlets_3d) == 0):
            # print "number of 1d is zero"
            sys.exit(
                "number of 3d_outlets are zero!! Please check your config.py file!")

        for i in range(0, len(self.inlets_3d)):
            inlet_nodes_3d = self.inlets_3d[i]
            area3d = 0.0
            for node in inlet_nodes_3d:
                n = node.GetSolutionStepValue(NORMAL)
                a = math.sqrt(n[0] ** 2 + n[1] ** 2 + n[2] ** 2)
                area3d += a
            self.inlet_areas_3d.append(area3d)
            # print "area3d_inlet", area3d

        for i in range(0, len(self.outlets_3d)):
            inlet_nodes_3d = self.outlets_3d[i]
            area3d = 0.0
            for node in inlet_nodes_3d:
                n = node.GetSolutionStepValue(NORMAL)
                a = math.sqrt(n[0] ** 2 + n[1] ** 2 + n[2] ** 2)
                area3d += a
            self.outlet_areas_3d.append(area3d)
            print "area_inlet_3d", area3d
            #raw_input()

	import FitAB
	self.fitters_1d = []
	self.fitters_3d = []
	self.fitters_inlet_nodes=[]
	self.fitters_outlet_nodes=[]
	
	for i in range(0,len(self.outlets_1d)):
	  self.fitters_1d.append(FitAB.Fitter(self.outlets_1d[i][0].Id))        
	  self.fitters_3d.append(FitAB.Fitter(self.outlets_1d[i][0].Id))
	
	self.fixed_flow_nodes =[]
	for node in self.model_part_1d.Nodes:
	  if(node.IsFixed(FLOW) == True and node.GetSolutionStepValue(FLAG_VARIABLE) == 0):
	    self.fixed_flow_nodes.append(node)  	  
	for i in range(0,len(self.fixed_flow_nodes)):
	  self.fitters_inlet_nodes.append(FitAB.Fitter(self.fixed_flow_nodes[i].Id))
	  print self.fixed_flow_nodes[i].Id

	self.fixed_outlet_nodes=[]
	for node in self.model_part_1d.Nodes:
	  if(node.IsFixed(PRESSURE) == True):
	    self.fixed_outlet_nodes.append(node)	
	for i in range(0,len(self.fixed_outlet_nodes)):
	  self.fitters_outlet_nodes.append(FitAB.Fitter(self.fixed_outlet_nodes[i].Id))
	  print self.fixed_outlet_nodes[i].Id   
	
    def Initial_Contitions(self):
        #initial_pressure = 0  # TODO
        initial_pressure = config.dyastolic_pressure
        #print "Inicializo 3D"
        #print initial_pressure
        for i in range(0, len(self.inlets_1d)):
            inlet_nodes_1d = self.inlets_1d[i]
            print "3D-1D: inlet_nodes_1d [0].Id::::::>>>> ", inlet_nodes_1d[0].Id
            for i in range(0, len(self.outlets_1d)):
                outlet_nodes_1d = self.outlets_1d[i]
                pressinlet3D = outlet_nodes_1d[0].GetSolutionStepValue(PRESSURE)
                #print pressinlet3D
                print "3D-1D: outlet_nodes_1d[0].Id::::::>>>> ", outlet_nodes_1d[0].Id
                for i in range(0, len(self.inlets_3d)):
                    inlet_nodes_3d = self.inlets_3d[i]
                    area3d = self.inlet_areas_3d[i]
                    print "area_inlet_3d", area3d
                    directions = self.inlet_velocity_directions[i]
                    radio3d = math.sqrt(area3d * 3.1416)
                    vel1d = inlet_nodes_1d[0].GetSolutionStepValue(FLOW) / area3d
                    k = 0
                    print vel1d
                    # Impongo velocidad y la presion
                    for node in inlet_nodes_3d:
                        n = node.GetSolutionStepValue(NORMAL)
                        a = math.sqrt(n[0] ** 2 + n[1] ** 2 + n[2] ** 2)
                        orientation = directions[k]
                        node.SetSolutionStepValue(VELOCITY, 0, directions[k] * vel1d)
                        # node.SetSolutionStepValue(PRESSURE, 0, pressinlet3D)
                        k = k + 1
            # raw_input()

    def Transfer1D_to_3D(self):
        # ARCHIVE TO SET :::::::::::::::::::::::::::>>>>>>>>>>>>>> VARIABLES
        # import config_full
        initial_pressure = config.dyastolic_pressure
        #initial_pressure = 0
	print initial_pressure
        print "Transfer1D_to_3D"
        for i in range(0, len(self.inlets_1d)):
            inlet_nodes_1d = self.inlets_1d[i]
            #print "NODO 1D-3D:: inlet_nodes_1d[0].Id::::::inlet 1D coupled with the 3D inlet>>>> ",inlet_nodes_1d[0].Id
            # print "--"
            # raw_input()
            for i in range(0, len(self.inlets_3d)):
                inlet_nodes_3d = self.inlets_3d[i]
                area3d = self.inlet_areas_3d[i]
                directions = self.inlet_velocity_directions[i]
                radio3d = math.sqrt(area3d * 3.1416)
                vel1d = inlet_nodes_1d[0].GetSolutionStepValue(FLOW) / area3d
                # NOTA: REVISAR EL CAUDAL QUE ESTAMOS IMPONIENDO DE ENTRADA El hecho de considerar q/2 es debido a que el caudal que imponemos
                # pessinlet1D =
                # inlet_nodes_1d[0].GetSolutionStepValue(PRESSURE)
                k = 0
                # Impongo velocidad y la presion
                for node in inlet_nodes_3d:
                    n = node.GetSolutionStepValue(NORMAL)
                    a = math.sqrt(n[0] ** 2 + n[1] ** 2 + n[2] ** 2)
                    orientation = directions[k]
                    node.SetSolutionStepValue(VELOCITY, 0, directions[k] * vel1d)
                    k = k + 1
            #print "area3d",area3d ,"area1d",inlet_nodes_1d[0].GetSolutionStepValue(NODAL_AREA), "del nodo", inlet_nodes_1d[0].Id
            #print "Estoy fijando la velocidad en el inlet del 3D que proviene del nodo " , inlet_nodes_1d[0].Id, " del 1D. La velocidad que estoy fijando es"
            #print "vel1d",vel1d, "del nodo", inlet_nodes_1d[0].Id
            

        for i in range(0, len(self.outlets_1d)):
            outlet_nodes_1d = self.outlets_1d[i]
            pressinlet3D = outlet_nodes_1d[0].GetSolutionStepValue(PRESSURE)
            # pressinlet3D = 0
            #print "NODO 1D-3D: outlet_nodes_1d[0].Id::::::outlet 1D coupled with the 3D outlet>>>> ",outlet_nodes_1d[0].Id
            #print "pressinlet3D ",pressinlet3D
            #raw_input()
            # double beta = E0*thickness0*1.77245385/(1.0-nu0*nu0);
            beta = outlet_nodes_1d[0].GetSolutionStepValue(BETA)
            #beta = ((outlet_nodes_1d[0].GetSolutionStepValue(YOUNG_MODULUS) * outlet_nodes_1d[0].GetSolutionStepValue(THICKNESS) * math.sqrt(math.pi)) / (
                #1 - (outlet_nodes_1d[0].GetSolutionStepValue(POISSON_RATIO) * outlet_nodes_1d[0].GetSolutionStepValue(POISSON_RATIO))))
            A = outlet_nodes_1d[0].GetSolutionStepValue(NODAL_AREA)
            #print "AREA_1D", A
            A0 = outlet_nodes_1d[0].GetValue(NODAL_AREA)
            press = initial_pressure + beta * \
                (math.sqrt(A) - math.sqrt(A0)) / \
                A0  # math.sqrt(A/A0)*beta - beta
            #press = initial_pressure + beta * \
                #(math.sqrt(A) - math.sqrt(A0)) / \
                #A0  # math.sqrt(A/A0)*beta - beta
            print "prress calculated",press 
            print "press solver", pressinlet3D
            for i in range(0, len(self.outlets_3d)):
                outlet_nodes_3d = self.outlets_3d[i]
                area3d = self.outlet_areas_3d[i]
                for node in outlet_nodes_3d:
                    node.SetSolutionStepValue(PRESSURE, 0, press)     
            #print "area3d",area3d ,"area1d",outlet_nodes_1d[0].GetSolutionStepValue(NODAL_AREA)
            #print "Estoy fijando la presion en el outlet del 3D que proviene del nodo " , outlet_nodes_1d[0].Id, " del 1D. La presion que estoy fijando es"
            #print "pression",press, "del nodo", outlet_nodes_1d[0].Id

    def Transfer3D_to_1D(self):
        initial_pressure = config.dyastolic_pressure
        inlet_flow = (self.inlets_1d[0])[0].GetSolutionStepValue(FLOW)
        self.outlets_1d[0][0].SetSolutionStepValue(FLOW,0,inlet_flow)

        #outlet_area = self.outlets_3d[0][0].GetSolutionStepValue(NODAL_AREA)
        #outlet_beta = self.outlets_3d[0][0].GetSolutionStepValue(BETA)
        outlet_pressure = self.outlets_1d[0][0].GetSolutionStepValue(PRESSURE)

        inlet_beta = self.inlets_1d[0][0].GetSolutionStepValue(BETA)
        inlet_A0  = self.inlets_1d[0][0].GetValue(NODAL_AREA)
        Ainlet_to_prescribe = ((((outlet_pressure-initial_pressure) * inlet_A0)/inlet_beta) + math.sqrt(inlet_A0))**2  # A0*(avg_press/beta + 1)**2
        #self.outlets_1d[0][0].SetSolutionStepValue(NODAL_AREA,0,Ainlet_to_prescribe)
        
        print "inlet node = ",self.inlets_1d[0][0].Id
        print "outlet node = ",self.outlets_1d[0][0].Id
		
		
    #def Transfer3D_to_1D(self):
        #for i in range(0, len(self.inlets_3d)):
            #inlet_nodes_1d = self.inlets_1d[i]
            #inlet_nodes_3d = self.inlets_3d[i]
            #area3d = self.inlet_areas_3d[i]
            #outlet_nodes_1d = self.outlets_1d[i]
            #press_3d = 0.0
            #counter = 0.0
            #for node in inlet_nodes_3d:
                #press_3d += node.GetSolutionStepValue(PRESSURE)
                #counter += 1.0
            #avg_press = press_3d / counter
            ##beta = inlet_nodes_1d[0].GetSolutionStepValue(YOUNG_MODULUS) * inlet_nodes_1d[
                ##0].GetSolutionStepValue(THICKNESS) * math.sqrt(math.pi)
            ##beta = beta / \
                ##(1.0 - (inlet_nodes_1d[0].GetSolutionStepValue(
                    ##POISSON_RATIO) * inlet_nodes_1d[0].GetSolutionStepValue(POISSON_RATIO)))
            ##print "beta", beta
            #beta = inlet_nodes_1d[0].GetSolutionStepValue(BETA)
            #print "beta", beta
            #print "inlet node is ", inlet_nodes_1d[0].Id
            #beta = inlet_nodes_1d[0].GetSolutionStepValue(BETA)
            #A0 = inlet_nodes_1d[0].GetValue(NODAL_AREA)
            #Ainlet_to_prescribe = (((avg_press * A0)/beta) + math.sqrt(A0))**2  # A0*(avg_press/beta + 1)**2
            #print "in Transfer AREA 3D_to_1D----->", Ainlet_to_prescribe, " to node" , inlet_nodes_1d[0].Id
            #print "Area del 1D Original", A0 
            #print "Area del 1D Modificada", inlet_nodes_1d[0].GetSolutionStepValue(NODAL_AREA)
            
           
            ## A0 = outlet_nodes_1d[0].GetValue(NODAL_AREA)
            ## outlet_nodes_1d[0].SetSolutionStepValue(NODAL_AREA ,0, A0)
            #inlet_nodes_1d[0].SetSolutionStepValue(NODAL_AREA, 0, Ainlet_to_prescribe)

            ## compute flow on inlet
            #flow = 0.0
            #for node in inlet_nodes_3d:
                #normal = node.GetSolutionStepValue(NORMAL)
                #vel = node.GetSolutionStepValue(VELOCITY)
                #flow += normal[0] * vel[0] + normal[
                    #1] * vel[1] + normal[2] * vel[2]
            ## self.flow_3d_in = str(flow)
            #print "TRASFER to 3D to 1D:::flow entering= ",flow
            #print "Flow de 1D inlet",  inlet_nodes_1d[0].GetSolutionStepValue(FLOW)
            ## flow_aux=  flow

        #for i in range(0, len(self.outlets_3d)):
            #outlet_nodes_1d = self.outlets_1d[i]
            #outlet_nodes_3d = self.outlets_3d[i]
            #area3d = self.outlet_areas_3d[i]

            ## compute flow on 3D outlet
            #flow = 0.0
            #for node in outlet_nodes_3d:
                #normal = node.GetSolutionStepValue(NORMAL)
                #vel = node.GetSolutionStepValue(VELOCITY)
                #flow += normal[0] * vel[0] + normal[
                    #1] * vel[1] + normal[2] * vel[2]

            ## assign flow to outlet
            #outlet_nodes_1d[0].SetSolutionStepValue(FLOW, 0, -flow)
            #print " velocity del 3D", vel
            #print " velocity del 1D", vel
            #print "in Transfer Flow 3D_to_1D----->", -flow, " to node" , outlet_nodes_1d[0].Id

        # raw_input()
#-------------------------------------------------------------------------

# FIT VALUES (WRITE FILE)

#-------------------------------------------------------------------------
        

#-------------------------------------------------------------------------

# Setting Contitions 3d

#-------------------------------------------------------------------------

    def Setting3d(self):
        # print self
     # print len(self.model_part_3d.Conditions)
        for node in self.model_part_3d.Nodes:
            node.Free(VELOCITY_X)
            node.Free(VELOCITY_Y)
            node.Free(VELOCITY_Z)
            node.Free(PRESSURE)
            node.SetSolutionStepValue(VISCOSITY, 0, 0.0035 / 1060.0)
            node.SetSolutionStepValue(DENSITY, 0, 1060.0)
            node.SetSolutionStepValue(PRESSURE,0,8000)
        # set inlet
        for cond in self.model_part_3d.Conditions:
            if(cond.Properties.Id == 100):  # inlet
                for node in cond.GetNodes():
                    node.Fix(VELOCITY_X)
                    node.Fix(VELOCITY_Y)
                    node.Fix(VELOCITY_Z)
                    node.SetSolutionStepValue(FLAG_VARIABLE, 0, 100.0)
            # if(cond.Properties.Id > 1000): ##outlet

                # for node in cond.GetNodes():

                    # node.Fix(PRESSURE)

                    # node.SetSolutionStepValue(FLAG_VARIABLE,0,cond.Properties.Id)

            # if(cond.Properties.Id > 1001): ##outlet

                # for node in cond.GetNodes():

                    # node.Fix(PRESSURE)

                    # node.SetSolutionStepValue(FLAG_VARIABLE,0,1002.0)

        # set sides (overwrites the inlet)
        for cond in self.model_part_3d.Conditions:
            if(cond.Properties.Id == 1):  # sides --> note that this is done in an outer separated loop!!
                for node in cond.GetNodes():
                    node.Fix(VELOCITY_X)
                    node.Fix(VELOCITY_Y)
                    node.Fix(VELOCITY_Z)
                    node.SetSolutionStepValue(FLAG_VARIABLE, 0, 0.0)

        # set output (overwrites the others)

        for cond in self.model_part_3d.Conditions:
            if(cond.Properties.Id > 1000):  # outlet
                for node in cond.GetNodes():
                    # if(not node.IsFixed(VELOCITY_X)):
                    node.Fix(PRESSURE)
                    node.SetSolutionStepValue(FLAG_VARIABLE, 0, cond.Properties.Id)

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
        for node in self.model_part_3d.Nodes:
            y = node.GetSolutionStepValue(Y_WALL, 0)
            node.SetValue(Y_WALL, y)

        counter = 0.0

        for node in self.model_part_3d.Nodes:
            if(node.IsFixed(PRESSURE)):
                counter += 1.0


        # print "n pressure nodes ",counter
       
    def FitValues_Inlet(self):
      for i in range(0,len(self.fixed_flow_nodes)):
	node=self.fixed_flow_nodes[i]
	for fitter in self.fitters_inlet_nodes:
	  fitter.AddPin(node[0].GetSolutionStepValue(PRESSURE))
	  fitter.AddQ(node[0].GetSolutionStepValue(FLOW))	  
	  
      for i in range(0,len(self.fixed_outlet_nodes)):
	node=self.fixed_outlet_nodes[i]
	for fitter in self.fitters_inlet_nodes:
	  fitter.AddPout(node[0].GetSolutionStepValue(PRESSURE))
    #raw_input()
	  
    def FitValues_1d(self):     
      for i in range(0, len(self.inlets_1d)):
	node = self.inlets_1d[i]
	##print "3D-1D: inlet_nodes_1d [0].Id::::::>>>> ", node[0].Id
	#nodewrite = node[0].Id
	##ToWrite = str(indextowrite) + " "
	##ToWrite = str(nodewrite) + " " + str(total_time) + " " + str(node[0].GetSolutionStepValue(PRESSURE)) + " "
	##ToWrite += str(node[0].GetSolutionStepValue(FLOW)) + " " + str(node[0].GetSolutionStepValue(NODAL_AREA)) + "\n"
	#ToWrite = str(node[0].GetSolutionStepValue(FLOW)) + "\n"
	#ffit[0].write(ToWrite)
	#AB_Matriz[row_AB][0]= total_time
	#AB_Matriz[row_AB][1]= node[0].GetSolutionStepValue(PRESSURE)
	#AB_Matriz[row_AB][2]= node[0].GetSolutionStepValue(FLOW)
	#AB_Matriz[row_AB][3]= node[0].GetSolutionStepValue(NODAL_AREA)
	for fitter in self.fitters_1d:
	  fitter.AddPin(node[0].GetSolutionStepValue(PRESSURE))
	j=0
	for j in range(0, len(self.outlets_1d)):
	  #j=j+1
	  #row_AB=row_AB+1
	  node = self.outlets_1d[j]	 
	  ##print "3D-1D: outlet_nodes_1d[0].Id::::::>>>> ", node[0].Id        
	  #nodewrite = node[0].Id
	  ##ToWrite = str(indextowrite) + " "
	  ##ToWrite = str(nodewrite) + " " + str(total_time) + " " + str(node[0].GetSolutionStepValue(PRESSURE)) + " "
	  ##ToWrite += str(node[0].GetSolutionStepValue(FLOW)) + " " + str(node[0].GetSolutionStepValue(NODAL_AREA)) + "\n"
	  #ToWrite = str(node[0].GetSolutionStepValue(FLOW)) + "\n"
	  #ffit[j].write(ToWrite)
	  #AB_Matriz[row_AB][0]= total_time
	  #AB_Matriz[row_AB][1]= node[0].GetSolutionStepValue(PRESSURE)
	  #AB_Matriz[row_AB][2]= node[0].GetSolutionStepValue(FLOW)
	  #AB_Matriz[row_AB][3]= node[0].GetSolutionStepValue(NODAL_AREA)
	  self.fitters_1d[j].AddPout(node[0].GetSolutionStepValue(PRESSURE))
	  self.fitters_1d[j].AddQ(node[0].GetSolutionStepValue(FLOW))
		        
    def FitValues_3d(self):
      for i in range(0, len(self.inlets_1d)):
	inlet_nodes_1d = self.inlets_1d[i]
	inlet_nodes_3d = self.inlets_3d[i]
	press_3d = 0.0
	flow_3d = 0.0
	area3d=0.0
	counter = 0.0
	print "inlet", inlet_nodes_1d[0]
	for node in inlet_nodes_3d:
	  press_3d += node.GetSolutionStepValue(PRESSURE)
	  counter += 1.0
	  avg_press = press_3d / counter	
	#area3d=self.inlet_areas_3d
	for fitter in self.fitters_3d:
	  fitter.AddPin(avg_press)
	  #self.fitters_3d[i].AddPin(avg_press)   	
	#print "pressure", avg_press
	j=0
	for j in range(0, len(self.outlets_1d)):
	  outlet_nodes_1d = self.outlets_1d[j]
	  outlet_nodes_3d = self.outlets_3d[j]
	  area3d = self.outlet_areas_3d[j]
	  # compute flow on 3D outlet
	  flow = 0.0
	  for node in outlet_nodes_3d:
	      normal = node.GetSolutionStepValue(NORMAL)
	      vel = node.GetSolutionStepValue(VELOCITY)
	      flow += normal[0] * vel[0] + normal[1] * vel[1] + normal[2] * vel[2]         
	  self.fitters_3d[j].AddQ(flow)  
	  self.fitters_3d[j].AddPout(outlet_nodes_1d[0].GetSolutionStepValue(PRESSURE))  
	  #print "node outlet", outlet_nodes_1d[0]
	  #print "flow", flow
	  #print "pressure", outlet_nodes_1d[0].GetSolutionStepValue(PRESSURE)
      
    def Fit_ABValues_1D(self):
      ffit_1d=[]
      k=0
      for fitter in self.fitters_1d:
	  ffit_1d.append(k)
	  results = str("Fitter_1D_" + str(k) + "_.txt")
	  ffit_1d[k] = open(results, 'w')
	  ffit_1d[k].write("NODE-PRESSURE_OUTLET_1D   -    PRESSURE_INLET_1D        -      FLOW_1D) \n")
	  [node_id,A,B] = fitter.DoFitting_1D(k,ffit_1d)
	  k=k+1
	  print "node_id = ",node_id
	  print "A = ",A
	  print "B = ",B

    def Fit_ABValues_3D(self): 
      ffit_3d=[]
      k=0
      for fitter in self.fitters_3d:
	  ffit_3d.append(k)
	  results = str("Fitter_3D_" + str(k) + "_.txt")
	  ffit_3d[k] = open(results, 'w')		
	  ffit_3d[k].write("NODE-PRESSURE_OUTLET_3D   -    PRESSURE_INLET_3D        -      FLOW_3D) \n")
	  [node_id,A,B] = fitter.DoFitting_3D(k,ffit_3d)
	  k=k+1
	  print "node_id = ",node_id
	  print "A = ",A
	  print "B = ",B
	 # ...para los fitters_3d
      
      
      #i=0
      #j=0

      #GradientPress=[]
      #dd=len(range(row_AB))
      #kk=dd/3      
      #for i in range(kk):
	#GradientPress.append([0]*1)   #Variables to save
      
      #A_Values=[]
      #for i in range(2): #nodes_number
	#A_Values.append([0]*2)   #Variables to save
      
      #B_values=[]
      #for i in range(2): #nodes_number
	#B_values.append([0]*1)   #Variables to save
      
      #print len(range(row_AB))
      ##AB_Matriz[][0]=time
      ##AB_Matriz[][1]=Pressure
      ##AB_Matriz[][2]=flow
      ##AB_Matriz[][3]=NodalArea
      #j=0
      #i=0
      #while i < len(range(row_AB)):
	##print i	
	##print "pressure", AB_Matriz[i][1]
	##print "pressure", AB_Matriz[i+1][1]
	##print "pressure", AB_Matriz[i+2][1]	
	#GradientPress[j]=AB_Matriz[i+1][1]-AB_Matriz[i][1]
	#print GradientPress
	##FlowFlow=AB_Matriz[i+1][2]*AB_Matriz[i+1][2]
	##A=(GradientPress-FlowFlow)/AB_Matriz[i+1][2]
	##B=(GradientPress-AB_Matriz[i+1][2])/FlowFlow
	##ToWrite = str(GradientPress) + " " + str(FlowFlow) + "\n" 
	##ToWrite = str(A) + " " + str(B) + "\n"
	##ffit_A[0].write(ToWrite)
	#i=i+3
	#j=j+1

      #TotalGradientPress=0
      #for i in range(kk):
	#TotalGradientPress=TotalGradientPress+GradientPress[i]
      
      #a11=0
      #a21=0
      #a12=0
      #a22=0
      #b1=0
      #b2=0
      #for i in range(0,len(row_AB)):
	#Q2 = AB_Matriz[i][2]*AB_Matriz[i][2]
	#a11=a11+(Q2)
	#a12=a12+(Q2*AB_Matriz[i][2])
	#a22=a22+(Q2*Q2)
	#b1=b1+(AB_Matriz[i][2]*GradientPress[i])
	#b2=b2+(Q2*GradientPress[i])
      #a21 = a12
		
      #A_Values[0][0]= a11 
      ## q2 
      #A_Values[0][1]= a12 
      ## q3
      #A_Values[1][0]= a12 
      ## q3
      #A_Values[1][1]= a22 
      ## q4
      #B_values[0]=b1 
      ##PQ
      #B_values[1]=b2 
      ##PQ2
            
      #print a11
      #print a12 
      #print a21
      #print a22
      #print b1
      #print b2
      
      #ToWrite = str(a11) + " " + str(a12) + " " + str(a21) + " " + str(a22) + " " 
      #ToWrite += str(b1) + " " + str(b2) + " " + str(TotalGradientPress) + "\n"
      #ffit_A[0].write(ToWrite)
      
      ##transposed = []
      ##for i in range(2):
	##transposed.append([row[i] for row in A_Values])
    
      #B=(b2-a21/a11*b1)/(a22-a12*a21/a11)
      #A=b1/a11-a12/a11*B
      #ToWrite = str(A) + " " + str(B) + "\n"
      #ffit_A[0].write(ToWrite)
      
  
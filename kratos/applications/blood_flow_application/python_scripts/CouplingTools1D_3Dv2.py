from KratosMultiphysics import *
from KratosMultiphysics.BloodFlowApplication import *
CheckForPreviousImport()

import math

class TransferTools:
  
    def __init__( self, model_part_1d, model_part_3d ):
	self.model_part_1d = model_part_1d
	self.model_part_3d = model_part_3d
    
    
    
    ################################################################3
    def Initialize(self):
	self.inlets_1d = []
	self.outlets_1d = []
	self.inlets_3d = []
	self.outlets_3d = []
	self.inlet_areas_3d = []
	self.outlet_areas_3d = []
	
	##detect 1d inlets
	for i in range(100,101):
	  nfound = 0
	  aux = []
	  for node in self.model_part_1d.Nodes:
	      if (node.GetSolutionStepValue(FLAG_VARIABLE) == i):
		  aux.append(node)
		  node.Fix(NODAL_AREA)
		  #node.Fix(FLOW)
	  if(len(aux) != 0):
	      self.inlets_1d.append(aux)
	  else:
	      break
	    
	##detect 3d inlets
	for i in range(100,101):
	    nfound = 0
	    aux = []
	    for node in self.model_part_3d.Nodes:
		if (node.GetSolutionStepValue(FLAG_VARIABLE) == i):
		    aux.append(node)
		    node.Fix(VELOCITY_X)
		    node.Fix(VELOCITY_Y)
		    node.Fix(VELOCITY_Z)
	    if(len(aux) != 0):
		self.inlets_3d.append(aux)
	    else:
		break	
		
	if(len(self.inlets_1d) != len(self.inlets_3d)):
	    print "number of 1d and 3d inlets is different"
	    err
	  
	  
	  
	##detect 1d outlets
	for i in range(101,199):
	  nfound = 0
	  aux = []
	  for node in self.model_part_1d.Nodes:
	      if (node.GetSolutionStepValue(FLAG_VARIABLE) == i):
		  aux.append(node)
		  node.Fix(FLOW)
	  if(len(aux) != 0):
	      self.outlets_1d.append(aux)
	  else:
	      break
	    
	for i in range(101,199):
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
	    
	#try len(self.inlets_1d) == len(self.inlets_3d):
	  #print "found ",len(inlets_1d)," domains"
	#except RuntimeError:
	  #print "number of 1d and 3d inlets is different"    
	  
	  
	#compute normals and 3d areas
	BodyNormalCalculationUtils().CalculateBodyNormals(self.model_part_3d,3)
	for i in range(0,len(self.inlets_3d)):
	    inlet_nodes_3d = self.inlets_3d[i]
	    area3d = 0.0
	    for node in inlet_nodes_3d:
	      n = node.GetSolutionStepValue(NORMAL)
	      a = math.sqrt(n[0]**2 +  n[1]**2 + n[2]**2 )
	      area3d += a  
	    self.inlet_areas_3d.append(area3d)
	    
	for i in range(0,len(self.outlets_3d)):
	    inlet_nodes_3d = self.outlets_3d[i]
	    area3d = 0.0
	    for node in inlet_nodes_3d:
	      n = node.GetSolutionStepValue(NORMAL)
	      a = math.sqrt(n[0]**2 +  n[1]**2 + n[2]**2 )
	      area3d += a  
	    self.outlet_areas_3d.append(area3d)	
	    
	print "inlet areas = ",self.inlet_areas_3d
	print "outlet areas = ",self.outlet_areas_3d
	
    
    
    
    
    def Transfer1D_to_3D( self  ):      
	for i in range(0,len(self.inlets_3d)):  
	  inlet_nodes_1d = self.inlets_1d[i]
	  inlet_nodes_3d = self.inlets_3d[i]
	  area3d = self.inlet_areas_3d[i]
	  
	  vel1d = inlet_nodes_1d[0].GetSolutionStepValue(FLOW) / area3d
	  for node in inlet_nodes_3d:
	    node.SetSolutionStepValue( VELOCITY_X, 0, vel1d)
	    
	for i in range(0,len(self.outlets_3d)):  
	  outlet_nodes_1d = self.outlets_1d[i]
	  outlet_nodes_3d = self.outlets_3d[i]
	  beta = outlet_nodes_1d[0].GetSolutionStepValue(YOUNG_MODULUS) * outlet_nodes_1d[0].GetSolutionStepValue(THICKNESS) * math.sqrt(math.pi)
	  A= outlet_nodes_1d[0].GetSolutionStepValue(NODAL_AREA)
	  A0 = outlet_nodes_1d[0].GetValue(NODAL_AREA)
	  press = math.sqrt(A/A0)*beta - beta          
	  for node in outlet_nodes_3d:
	      node.SetSolutionStepValue( PRESSURE, 0, press)
      
      
    def Transfer3D_to_1D( self ):
      for i in range(0,len(self.inlets_3d)):  
	  inlet_nodes_1d = self.inlets_1d[i]
	  inlet_nodes_3d = self.inlets_3d[i]
	  area3d = self.inlet_areas_3d[i]
	  
	  press_3d = 0.0
	  counter = 0.0
	  for node in inlet_nodes_3d:
	    press_3d += node.GetSolutionStepValue( PRESSURE)
	    counter += 1.0
	  avg_press = press_3d/counter
	  
	  #TODO: make it to read from the input
	  #print "inlet_nodes_1d[0].GetSolutionStepValue(YOUNG_MODULUS)",inlet_nodes_1d[0].GetSolutionStepValue(YOUNG_MODULUS)
	  #print "inlet_nodes_1d[0].GetSolutionStepValue(THICKNESS)",inlet_nodes_1d[0].GetSolutionStepValue(THICKNESS)
	  #print inlet_nodes_1d[0].Id
	  #beta = inlet_nodes_1d[0].GetSolutionStepValue(YOUNG_MODULUS) * inlet_nodes_1d[0].GetSolutionStepValue(THICKNESS) * math.sqrt(math.pi)
	  #A = area3d*(avg_press/beta + 1)**2
	  #inlet_nodes_1d[0].SetSolutionStepValue(NODAL_AREA ,0, A) 

	  
      for i in range(0,len(self.outlets_3d)):  
	  outlet_nodes_1d = self.outlets_1d[i]
	  outlet_nodes_3d = self.outlets_3d[i]
	  area3d = self.outlet_areas_3d[i]
	    
	  #assign flow to the outlet
	  flow = 0.0
	  for node in outlet_nodes_3d:
	    normal = node.GetSolutionStepValue(NORMAL);
	    vel = node.GetSolutionStepValue(VELOCITY)
	    flow += normal[0]*vel[0] + normal[1]*vel[1] + normal[2]*vel[2]
	  print "flow = ",flow
	    
	  outlet_nodes_1d[0].SetSolutionStepValue(FLOW, 0, flow)
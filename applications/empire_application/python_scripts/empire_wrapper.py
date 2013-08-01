# import Kratos libraries
from KratosMultiphysics import *
from KratosMultiphysics.EmpireApplication import *
import os
from ctypes import *
CheckForPreviousImport()

class EmpireWrapper:

    def __init__(self,fluid_model_part):
        self.libempire_api 	  = cdll.LoadLibrary(os.environ['EMPIRE_API_LIBSO_ON_MACHINE'])
	self.interface_model_part = ModelPart("InterfaceModelPart");   
	self.ale_wrapper_process  = ALEWrapperProcess(fluid_model_part,self.interface_model_part)	

    def initialize(self):
	self.ale_wrapper_process.ExtractInterface()
	print "extracted interface"

    def recvDisplacement(self): 
	size = self.interface_model_part.GetNodes().Size() * 3
	c_displacements = (c_double * size)(0)	
	self.libempire_api.EMPIRE_API_recvDataField("defaultField", size, c_displacements)
	
	disp = Vector(3)
	i = 0
	for nodes in self.interface_model_part.GetNodes():
		disp = [c_displacements[3*i+0],c_displacements[3*i+1],c_displacements[3*i+2]]	
		nodes.SetSolutionStepValue(DISPLACEMENT,disp)
		print nodes.GetSolutionStepValue(DISPLACEMENT)
		i = i + 1

    def sendForces(self):
	forces = []
	self.ale_wrapper_process.ExtractForcesFromModelPart( forces )
	print forces

	c_forces = (c_double * len(forces))(*forces)
	c_size = len(c_forces)
	self.libempire_api.EMPIRE_API_sendDataField("defaultField", c_size , c_forces)

    def sendMesh(self):
	numNodes 	= []
	numElems 	= []
	nodes 		= []
	nodeIDs 	= []
	numNodesPerElem = []
	elems 		= []

	self.ale_wrapper_process.ExtractMeshInfo(numNodes,numElems,nodes,nodeIDs,numNodesPerElem,elems)

	c_numNodes = (c_int * len(numNodes))(*numNodes)
	c_numElems = (c_int * len(numElems))(*numElems)
	c_nodes = (c_double * len(nodes))(*nodes)
	c_nodeIDs = (c_int * len(nodeIDs))(*nodeIDs)
	c_numNodesPerElem = (c_int * len(numNodesPerElem))(*numNodesPerElem)
	c_elems = (c_int * len(elems))(*elems)	

	self.libempire_api.EMPIRE_API_sendMesh("defaultMesh",c_numNodes[0], c_numElems[0], c_nodes, c_nodeIDs, c_numNodesPerElem, c_elems)

    def recvConvergenceSignal(self):
	isConvergent = self.libempire_api.EMPIRE_API_recvConvergenceSignal()
	return isConvergent

# import libraries
from KratosMultiphysics import *
from KratosMultiphysics.EmpireApplication import *
import os
from ctypes import *
CheckForPreviousImport()

class EmpireWrapper:

    # -------------------------------------------------------------------------------------------------
    def __init__(self,model_part,libempire_api):
	self.libempire_api        = libempire_api
	self.model_part           = model_part 
	self.interface_model_part = ModelPart("InterfaceModelPart");  
	self.ale_wrapper_process  = ALEWrapperProcess(self.model_part,self.interface_model_part)
	########################
	self.utilities 		  = VariableUtils()
    # -------------------------------------------------------------------------------------------------	

    # -------------------------------------------------------------------------------------------------
    # extract interface nodes and conditions from the fluid model part which have a flag "IS_INTERFACE"
    def extractInterface(self):
	self.ale_wrapper_process.ExtractInterface()
    # -------------------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------------------
    # convert char* to python-compatible string
    def getUserDefinedText(self,stringText):
	empireFunc = self.libempire_api.EMPIRE_API_getUserDefinedText
	empireFunc.restype = c_char_p
	text = empireFunc(stringText)
	return text
    # -------------------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------------------
    def recvDisplacement(self): 
	# initialize vector storing the displacements
	size = self.interface_model_part.GetNodes().Size() * 3
	c_displacements = (c_double * size)(0)

	# receive displacements from empire	
	self.libempire_api.EMPIRE_API_recvDataField("defaultField", size, c_displacements)

	disp = Vector(3)
	i = 0
	# assign displacements to nodes of interface for current time step
	for nodes in (self.interface_model_part).Nodes:
		disp[0] = c_displacements[3*i+0]
		disp[1] = c_displacements[3*i+1]
		disp[2] = c_displacements[3*i+2]

		if(abs(disp[0]) < 1e-9):
			disp[0] = 0
		if(abs(disp[1]) < 1e-9):
			disp[1] = 0
		if(abs(disp[2]) < 1e-9):
			disp[2] = 0

		nodes.SetSolutionStepValue(DISPLACEMENT,0,disp)
		i = i + 1
    # -------------------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------------------
    def sendPressure(self):
	# extract pressure from interface nodes as result of the CFD simulation
	pressure = []
	self.ale_wrapper_process.ExtractPressureFromModelPart( pressure )

	# convert list containg the forces to ctypes
	c_pressure = (c_double * len(pressure))(*pressure)
	c_size = len(c_pressure)

	# send pressure to empire
	self.libempire_api.EMPIRE_API_sendDataField("defaultField", c_size , c_pressure)
    # -------------------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------------------
    def sendMesh(self):
	numNodes 	= []
	numElems 	= []
	nodes 		= []
	nodeIDs 	= []
	numNodesPerElem = []
	elems 		= []

	# extract interface mesh information
	self.ale_wrapper_process.ExtractMeshInfo(numNodes,numElems,nodes,nodeIDs,numNodesPerElem,elems)

	# convert python lists to ctypes required for empire-function call
	c_numNodes = (c_int * len(numNodes))(*numNodes)
	c_numElems = (c_int * len(numElems))(*numElems)
	c_nodes = (c_double * len(nodes))(*nodes)
	c_nodeIDs = (c_int * len(nodeIDs))(*nodeIDs)
	c_numNodesPerElem = (c_int * len(numNodesPerElem))(*numNodesPerElem)
	c_elems = (c_int * len(elems))(*elems)

	# send mesh information to empire
	self.libempire_api.EMPIRE_API_sendMesh("defaultMesh",c_numNodes[0], c_numElems[0], c_nodes, c_nodeIDs, c_numNodesPerElem, c_elems)
    # -------------------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------------------
    def recvConvergenceSignal(self):
	# receive convergence information from empire
	isConvergent = self.libempire_api.EMPIRE_API_recvConvergenceSignal()
	return isConvergent
    # -------------------------------------------------------------------------------------------------

#######################################################################################################
    # Functions required for Kratos-Kratos-FSI (whereas Kratos solves also the structural side)
    # -------------------------------------------------------------------------------------------------
    def recvPressure(self):
	# initialize pressure variable
	size = self.interface_model_part.GetNodes().Size() * 1
	c_pressure = (c_double * size)(0)

	# receive pressure from empire	
	self.libempire_api.EMPIRE_API_recvDataField("defaultField", size, c_pressure)

	pressure = Vector(1)
	i = 0
	# assign pressure to nodes of interface for current time step
	for nodes in (self.interface_model_part).Nodes:
		pressure[0] = c_pressure[i]
		nodes.SetSolutionStepValue(POSITIVE_FACE_PRESSURE,0,pressure[0]);
		i = i + 1
    # -------------------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------------------	
    def sendDisplacements(self):
	# extract displacements from interface nodes as result of the CFD simulation
	displacements = []
	self.ale_wrapper_process.ExtractDisplacementsFromModelPart( displacements )

	# convert list containg the displacements to ctypes
	c_displacements = (c_double * len(displacements))(*displacements)
	c_size = len(c_displacements)

	# send displacements to empire
	self.libempire_api.EMPIRE_API_sendDataField("defaultField", c_size , c_displacements)
    # -------------------------------------------------------------------------------------------------
#######################################################################################################

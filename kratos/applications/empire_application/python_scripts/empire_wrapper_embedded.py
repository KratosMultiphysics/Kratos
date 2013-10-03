# import libraries
from ctypes import *
import os
from KratosMultiphysics import *
from KratosMultiphysics.EmpireApplication import *
from KratosMultiphysics.IncompressibleFluidApplication import *

CheckForPreviousImport()

class EmpireWrapper:

    # -------------------------------------------------------------------------------------------------
    def __init__(self,model_part,libempire_api):
	self.libempire_api        = libempire_api
	self.model_part           = model_part 
	self.interface_model_part = ModelPart("InterfaceModelPart");  
        self.wrapper_process      = WrapperProcess(self.model_part,self.interface_model_part)
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
        self.wrapper_process.ExtractPressureFromModelPart( pressure )

	# convert list containg the forces to ctypes
	c_pressure = (c_double * len(pressure))(*pressure)
	c_size = len(c_pressure)

	# send pressure to empire
	self.libempire_api.EMPIRE_API_sendDataField("defaultField", c_size , c_pressure)
    # -------------------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------------------
    def sendForces(self):
	# extract pressure from interface nodes as result of the CFD simulation
	forces = []
        self.wrapper_process.ExtractForcesFromModelPart( forces )

	# convert list containg the forces to ctypes
	c_forces = (c_double * len(forces))(*forces)
	c_size = len(c_forces)

	# send pressure to empire
	self.libempire_api.EMPIRE_API_sendDataField("defaultField", c_size , c_forces)
    # -------------------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------------------
    def sendInteraceMesh(self):
	numNodes        = []
	numElems        = []
	nodeCoordinates = []
	nodeIDs         = []
	numNodesPerElem = []
	elements        = []
	sizes		= []

        # extract interface nodes and conditions from the fluid model part which have a flag "IS_INTERFACE"
        self.wrapper_process.ExtractInterface()

	# Extract interface mesh info from above extracted interface
	self.wrapper_process.ExtractMeshInfo(numNodes,numElems,nodeCoordinates,nodeIDs,numNodesPerElem,elements)
	
	print "------------------------------------------------------------------------"
	print "Size of the nodeIDS-Array to be sent to the fluid: ", numNodes[0]
	print "Size of the nodeCoordinates-Array to be sent to the fluid: ", numNodes[0]*3
	print "Size of the elements-Array to be sent to the fluid: ", numElems[0]*3
	print "------------------------------------------------------------------------"

	# receive number of nodes of structure mesh (EMPIRE communication requrires the following data formatting)
	sizes.append(numNodes[0])
	sizes.append(numElems[0])

	c_sizes = (c_double * 2)(*sizes)
	self.libempire_api.EMPIRE_API_sendSignal_double("Interface_SizeInfo", 2, c_sizes)

	c_nodeIDs         = (c_double * numNodes[0])(*nodeIDs)
	c_nodeCoordinates = (c_double * (numNodes[0]*3))(*nodeCoordinates)
	c_elements        = (c_double * (numElems[0]*3))(*elements)

	# send mesh
	self.libempire_api.EMPIRE_API_sendSignal_double("Interface_NodeIDs", numNodes[0] , c_nodeIDs)
	self.libempire_api.EMPIRE_API_sendSignal_double("Interface_NodeCoordinates", numNodes[0]*3 , c_nodeCoordinates)
	self.libempire_api.EMPIRE_API_sendSignal_double("Interface_Elements", (numElems[0]*3) , c_elements)
    # -------------------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------------------
    def recvInterfaceMesh(self):
	nodeIDs		= []        
	nodeCoordinates = []
        elements        = []
	numNodes        = []
	numElems        = []
	sizes		= []
	c_sizes		= []

        # receive number of nodes of structure mesh
        c_sizes = (c_double * 2)(*sizes)
        self.libempire_api.EMPIRE_API_recvSignal_double("Interface_SizeInfo", 2, c_sizes)

	# Save size information (EMPIRE communication requrires the following data formatting)
	numNodes.append(int(c_sizes[0]))
        numElems.append(int(c_sizes[1]))

	c_nodeIDs         = (c_double * numNodes[0])(*nodeIDs)
        c_nodeCoordinates = (c_double * (numNodes[0]*3))(*nodeCoordinates)
        c_elements        = (c_double * (numElems[0]*3))(*elements)

        # receive mesh data from empire
	# Note: EMPIREs send and recvSignal only processes doubles, i.e. also the ids will be seen as doubles --> conversion to int later
	self.libempire_api.EMPIRE_API_recvSignal_double("Interface_NodeIDs", numNodes[0], c_nodeIDs)        
	self.libempire_api.EMPIRE_API_recvSignal_double("Interface_NodeCoordinates", (numNodes[0]*3), c_nodeCoordinates)
        self.libempire_api.EMPIRE_API_recvSignal_double("Interface_Elements", (numElems[0]*3), c_elements)

        # create python list to work with in c++
        for i in range(0,numNodes[0]):
                nodeIDs.append(int(c_nodeIDs[i]))
        for i in range(0,numNodes[0]*3):
                nodeCoordinates.append(c_nodeCoordinates[i])
	for i in range(0,(numElems[0]*3)):
		elements.append(int(c_elements[i]))

        # create interface part with the received mesh information
        self.wrapper_process.CreateEmbeddedInterfacePart(nodeIDs, nodeCoordinates, elements)
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

	i = 0
	# assign pressure to nodes of interface for current time step
	for nodes in (self.interface_model_part).Nodes:
		pressure = c_pressure[i]
		nodes.SetSolutionStepValue(POSITIVE_FACE_PRESSURE,0,pressure);
		i = i + 1
    # -------------------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------------------
    def recvForces(self):
	# initialize pressure variable
	size = self.interface_model_part.GetNodes().Size() * 3
	c_forces = (c_double * size)(0)

	# receive pressure from empire	
	self.libempire_api.EMPIRE_API_recvDataField("defaultField", size, c_forces)

	pressure = Vector(3)
	i = 0
	# assign pressure to nodes of interface for current time step
	for nodes in (self.interface_model_part).Nodes:
		forces[0] = c_pressure[3*i+0]
		forces[1] = c_pressure[3*i+1]
		forces[2] = c_pressure[3*i+2]
		nodes.SetSolutionStepValue(FORCE,0,forces);
		i = i + 1

	# Necessary to add forces to RS of structure
	AssignNoSlipCondition().AssignNoSlipCondition2D(self.model_part)
    # -------------------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------------------	
    def sendDisplacements(self):
	# extract displacements from interface nodes as result of the CFD simulation
	displacements = []
        self.wrapper_process.ExtractDisplacementsFromModelPart( displacements )

	# convert list containg the displacements to ctypes
	c_displacements = (c_double * len(displacements))(*displacements)
	c_size = len(c_displacements)

	# send displacements to empire
	self.libempire_api.EMPIRE_API_sendDataField("defaultField", c_size , c_displacements)
    # -------------------------------------------------------------------------------------------------
#######################################################################################################

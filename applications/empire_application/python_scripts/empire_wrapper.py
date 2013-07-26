# import Kratos libraries
from KratosMultiphysics import *
from KratosMultiphysics.EmpireApplication import *
import os
from ctypes import *
CheckForPreviousImport()

# guelle

class EmpireWrapper:

    def __init__(self,fluid_model_part):
        self.interfaceDisplacementU = (c_double * 7)(0)
        self.isConvergent = 0
	self.size = 7
        self.libempire_api = cdll.LoadLibrary(os.environ['EMPIRE_API_LIBSO_ON_MACHINE'])
	self.ale_wrapper_process = ALEWrapperProcess(fluid_model_part)	

    def recvDisplacement(self): 
        self.libempire_api.EMPIRE_API_recvSignal_double("interfaceDisplacementU", self.size, self.interfaceDisplacementU)

        print "Received Displacement:"
        print 'Received interfaceDisplacementU: {}'.format(self.interfaceDisplacementU[3])

    def sendForces(self):
        #self.interfaceForceY[0] = 200
	forces_X = []
	forces_Y = []
	forces_Z = []
	self.ale_wrapper_process.ExtractForcesFromModelPart( forces_X , forces_Y , forces_Z )

	interfaceForceY = (c_double * len(forces_X))(*forces_X)
        self.libempire_api.EMPIRE_API_sendSignal_double("interfaceForceY", self.size, interfaceForceY)
        #print 'Sent interfaceForceY: {}'.format(interfaceForceY[0])

    def sendMesh(self):

 	#* \param[in] numNodes number of nodes
 	#* \param[in] numElems number of elements
 	#* \param[in] nodes coordinates of all nodes
 	#* \param[in] nodeIDs IDs of all nodes
 	#* \param[in] numNodesPerElem number of nodes per element
 	#* \param[in] elems connectivity table of all elements
	print "send Mesh"
	#FillArrays(fluid_model_part,structure_model_part,numNodes,numElems,nodes,nodeIDs,numNodesPerElem,elems)

	#self.libempire_api.EMPIRE_API_sendMesh(int numNodes, int numElems, double *nodes, int *nodeIDs,int *numNodesPerElem, int *elems)

    def recvConvergenceSignal(self):
	self.isConvergent = self.libempire_api.EMPIRE_API_recvConvergenceSignal()
	return self.isConvergent

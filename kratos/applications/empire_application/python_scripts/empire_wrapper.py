# import Kratos libraries
from KratosMultiphysics import *
from KratosEmpireApplication import *
import os
from ctypes import *
CheckForPreviousImport()

class EmpireWrapper:

    def __init__(self,fluid_model_part):
	self.interfaceDisplacementU = (c_double * 1)(0)
	self.interfaceForceY 	    = (c_double * 1)(0)
	self.numNodes		    = 
	self.size                   = 1
	self.isConvergent  	    = 0
	self.libempire_api 	    = cdll.LoadLibrary(os.environ['EMPIRE_API_LIBSO_ON_MACHINE'])

    def recvDisplacement(self): 
        self.libempire_api.EMPIRE_API_recvSignal_double("interfaceDisplacementU", self.size, self.interfaceDisplacementU)

        print "Received Displacement:"
        print 'Received interfaceDisplacementU: {}'.format(self.interfaceDisplacementU[0])

    def sendForces(self):
	self.interfaceForceY[0] = 200
        self.libempire_api.EMPIRE_API_sendSignal_double("interfaceForceY", self.size, self.interfaceForceY)
        print 'Sent interfaceForceY: {}'.format(self.interfaceForceY[0])

    def sendMesh(self):

 	#* \param[in] numNodes number of nodes
 	#* \param[in] numElems number of elements
 	#* \param[in] nodes coordinates of all nodes
 	#* \param[in] nodeIDs IDs of all nodes
 	#* \param[in] numNodesPerElem number of nodes per element
 	#* \param[in] elems connectivity table of all elements

	#FillArrays(fluid_model_part,structure_model_part,numNodes,numElems,nodes,nodeIDs,numNodesPerElem,elems)

	self.libempire_api.EMPIRE_API_sendMesh(int numNodes, int numElems, double *nodes, int *nodeIDs,int *numNodesPerElem, int *elems)

    def recvConvergenceSignal(self):
	self.isConvergent = self.libempire_api.EMPIRE_API_recvConvergenceSignal()
	return self.isConvergent

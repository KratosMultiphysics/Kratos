from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# import libraries
from ctypes import *
import os
from KratosMultiphysics import *
from KratosMultiphysics.EmpireApplication import *
from KratosMultiphysics.IncompressibleFluidApplication import *

CheckForPreviousImport()


class EmpireWrapper:

    # -------------------------------------------------------------------------------------------------
    def __init__(self, model_part, libempire_api):
        self.libempire_api = libempire_api
        self.model_part = model_part
        self.interface_model_part = ModelPart("InterfaceModelPart")
        self.interface_model_part.AddNodalSolutionStepVariable(VELOCITY)
        self.interface_model_part.AddNodalSolutionStepVariable(DISPLACEMENT)
        self.interface_model_part.AddNodalSolutionStepVariable(POSITIVE_FACE_PRESSURE)
        self.interface_model_part.AddNodalSolutionStepVariable(NEGATIVE_FACE_PRESSURE)
        self.interface_model_part.SetBufferSize(3)
        self.wrapper_process = WrapperProcess(self.model_part, self.interface_model_part, 3)

        # specifically needed for exchange of mesh information (constant mesh information)
        self.nodeIDs = []
        self.elements = []
        self.interfaceX0 = []
        self.interfaceY0 = []
        self.interfaceZ0 = []
    # -------------------------------------------------------------------------------------------------

    def getInterfacePart(self):
        return self.interface_model_part

    # -------------------------------------------------------------------------------------------------
    # convert char* to python-compatible string
    def getUserDefinedText(self, stringText):
        empireFunc = self.libempire_api.EMPIRE_API_getUserDefinedText
        empireFunc.restype = c_char_p
        text = empireFunc(stringText)
        return text
    # -------------------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------------------
    def recvDisplacements(self):
        # initialize vector storing the displacements
        size = self.interface_model_part.GetNodes().Size() * 3
        c_displacements = (c_double * size)(0)

        # receive displacements from empire
        self.libempire_api.EMPIRE_API_recvSignal_double("displacements", size, c_displacements)

        # extract delta_t from fluid model part to calculate the updated velocity
        delta_T = self.model_part.ProcessInfo[DELTA_TIME]

        # initialize vectors and lists
        disp = Vector(3)
        disp_old = Vector(3)
        vel = Vector(3)
        vel_old = Vector(3)
        nodeCoordinates = []
        displacement = []
        velocity = []

        # update coordinates of interface of current time step
        i = 0
        for nodes in (self.interface_model_part).Nodes:
            disp[0] = c_displacements[3 * i]
            disp[1] = c_displacements[3 * i + 1]
            disp[2] = c_displacements[3 * i + 2]

            displacement.append(disp[0])
            displacement.append(disp[1])
            displacement.append(disp[2])

            if(abs(disp[0]) < 1e-9):
                disp[0] = 0
            if(abs(disp[1]) < 1e-9):
                disp[1] = 0
            if(abs(disp[2]) < 1e-9):
                disp[2] = 0

            # Compute coordinates of new interface_model_part
            nodeCoordinates.append(self.interfaceX0[i] + disp[0])
            nodeCoordinates.append(self.interfaceY0[i] + disp[1])
            nodeCoordinates.append(self.interfaceZ0[i] + disp[2])

            # Set velocity of new interface_model_part
            disp_old = nodes.GetSolutionStepValue(DISPLACEMENT, 1)
            vel[0] = (disp[0] - disp_old[0]) / delta_T
            vel[1] = (disp[1] - disp_old[1]) / delta_T
            vel[2] = (disp[2] - disp_old[2]) / delta_T

            velocity.append(vel[0])
            velocity.append(vel[1])
            velocity.append(vel[2])

            i = i + 1

        # create interface part with the received mesh information
        self.wrapper_process.CreateEmbeddedInterfacePart(self.nodeIDs, nodeCoordinates, self.elements, displacement, velocity)

        # assign displacements and velocities to nodes of updated interface
        # i = 0
        # for nodes in (self.interface_model_part).Nodes:
            # Set displacement of updated interface_model_part
            # disp[0] = c_displacements[3*i]
            # disp[1] = c_displacements[3*i+1]
            # disp[2] = c_displacements[3*i+2]
            # if(abs(disp[0]) < 1e-9):
                # disp[0] = 0
            # if(abs(disp[1]) < 1e-9):
                # disp[1] = 0
            # if(abs(disp[2]) < 1e-9):
                # disp[2] = 0
            # nodes.SetSolutionStepValue(DISPLACEMENT,0,disp)
            # Set velocity of new interface_model_part
            # disp_old = nodes.GetSolutionStepValue(DISPLACEMENT,1)
            # vel[0] = (disp[0]-disp_old[0]) / delta_T
            # vel[1] = (disp[1]-disp_old[1]) / delta_T
            # vel[2] = (disp[2]-disp_old[2]) / delta_T
            # nodes.SetSolutionStepValue(VELOCITY,0,vel)
            # i = i + 1
    # -------------------------------------------------------------------------------------------------
    # -------------------------------------------------------------------------------------------------
    def sendPressure(self):
        # extract pressure from interface nodes as result of the CFD simulation
        pressure = []
        self.wrapper_process.ExtractPressureFromEmbeddedModelPart(pressure)

        # List pressure has the format: Node1.POSITIVE_FACE_PRESSURE, Node1.NEGATIVE_FACE_PRESSURE,
        #				Node2.POSITIVE_FACE_PRESSURE, Node2.NEGATIVE_FACE_PRESSURE,
        # such that array has size node_size * 2

        # convert list containg the forces to ctypes
        c_pressure = (c_double * len(pressure))(*pressure)
        c_size = len(c_pressure)

        # send pressure to empire
        self.libempire_api.EMPIRE_API_sendSignal_double("pressure", c_size, c_pressure)
    # -------------------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------------------
    def sendForces(self):
        # extract pressure from interface nodes as result of the CFD simulation
        forces = []
        self.wrapper_process.ExtractForcesFromModelPart(forces)

        # convert list containg the forces to ctypes
        c_forces = (c_double * len(forces))(*forces)
        c_size = len(c_forces)

        # send pressure to empire
        self.libempire_api.EMPIRE_API_sendSignal_double("force", c_size, c_forces)
    # -------------------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------------------
    def sendInterfaceMesh(self):
        numNodes = []
        numElems = []
        nodeCoordinates = []
        nodeIDs = []
        numNodesPerElem = []
        elements = []
        sizes = []

        # extract interface nodes and conditions from the fluid model part which have a flag "IS_INTERFACE"
        self.wrapper_process.ExtractInterface()

        # Extract interface mesh info from above extracted interface
        self.wrapper_process.ExtractMeshInfo(numNodes, numElems, nodeCoordinates, nodeIDs, numNodesPerElem, elements)

        print("------------------------------------------------------------------------")
        print("Size of the nodeIDS-Array to be sent to the fluid: ", numNodes[0])
        print("Size of the nodeCoordinates-Array to be sent to the fluid: ", numNodes[0] * 3)
        print("Size of the elements-Array to be sent to the fluid: ", numElems[0] * 3)
        print("------------------------------------------------------------------------")

        # receive number of nodes of structure mesh (EMPIRE communication requrires the following data formatting)
        sizes.append(numNodes[0])
        sizes.append(numElems[0])

        c_sizes = (c_double * 2)(*sizes)
        self.libempire_api.EMPIRE_API_sendSignal_double("Interface_SizeInfo", 2, c_sizes)

        c_nodeIDs = (c_double * numNodes[0])(*nodeIDs)
        c_nodeCoordinates = (c_double * (numNodes[0] * 3))(*nodeCoordinates)
        c_elements = (c_double * (numElems[0] * 3))(*elements)

        # send mesh
        self.libempire_api.EMPIRE_API_sendSignal_double("Interface_NodeIDs", numNodes[0], c_nodeIDs)
        self.libempire_api.EMPIRE_API_sendSignal_double("Interface_NodeCoordinates", numNodes[0] * 3, c_nodeCoordinates)
        self.libempire_api.EMPIRE_API_sendSignal_double("Interface_Elements", (numElems[0] * 3), c_elements)
    # -------------------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------------------
    def recvInterfaceMesh(self):
        nodeCoordinates = []
        numNodes = []
        numElems = []
        sizes = []
        c_sizes = []
        zero_list = []

        # receive number of nodes of structure mesh
        c_sizes = (c_double * 2)(*sizes)
        self.libempire_api.EMPIRE_API_recvSignal_double("Interface_SizeInfo", 2, c_sizes)

        # Save size information (EMPIRE communication requrires the following data formatting)
        numNodes.append(int(c_sizes[0]))
        numElems.append(int(c_sizes[1]))

        c_nodeIDs = (c_double * numNodes[0])(*self.nodeIDs)
        c_nodeCoordinates = (c_double * (numNodes[0] * 3))(*nodeCoordinates)
        c_elements = (c_double * (numElems[0] * 3))(*self.elements)

        # receive mesh data from empire
        # Note: EMPIREs send and recvSignal only processes doubles, i.e. also the ids will be seen as doubles --> conversion to int later
        self.libempire_api.EMPIRE_API_recvSignal_double("Interface_NodeIDs", numNodes[0], c_nodeIDs)
        self.libempire_api.EMPIRE_API_recvSignal_double("Interface_NodeCoordinates", (numNodes[0] * 3), c_nodeCoordinates)
        self.libempire_api.EMPIRE_API_recvSignal_double("Interface_Elements", (numElems[0] * 3), c_elements)

        # create python list to work with in c++
        for i in range(0, numNodes[0]):
            self.nodeIDs.append(int(c_nodeIDs[i]))
            self.interfaceX0.append(c_nodeCoordinates[3 * i])
            self.interfaceY0.append(c_nodeCoordinates[3 * i + 1])
            self.interfaceZ0.append(c_nodeCoordinates[3 * i + 2])
        for i in range(0, numNodes[0] * 3):
            nodeCoordinates.append(c_nodeCoordinates[i])
            zero_list.append(0)
        for i in range(0, (numElems[0] * 3)):
            self.elements.append(int(c_elements[i]))

        # create interface part with the received mesh information
        self.wrapper_process.CreateEmbeddedInterfacePart(self.nodeIDs, nodeCoordinates, self.elements, zero_list, zero_list)
    # -------------------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------------------
    def recvConvergenceSignal(self):
        # receive convergence information from empire
        isConvergent = self.libempire_api.EMPIRE_API_recvConvergenceSignal()
        return isConvergent
    # -------------------------------------------------------------------------------------------------

#
    # Functions required for Kratos-Kratos-FSI (whereas Kratos solves also the structural side)
    # -------------------------------------------------------------------------------------------------
    def recvPressure(self):
        # initialize pressure variable
        size = self.interface_model_part.GetNodes().Size() * 1
        c_pressure = (c_double * size)(0)

        # receive pressure from empire
        self.libempire_api.EMPIRE_API_recvSignal_double("pressure", size, c_pressure)

        i = 0
        # assign pressure to nodes of interface for current time step
        for nodes in (self.interface_model_part).Nodes:
            pressure = c_pressure[i]
            nodes.SetSolutionStepValue(NEGATIVE_FACE_PRESSURE, 0, -pressure);
            # nodes.SetSolutionStepValue(POSITIVE_FACE_PRESSURE,0,pressure);
            i = i + 1
    # -------------------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------------------
    def recvForces(self):
        # initialize pressure variable
        size = self.interface_model_part.GetNodes().Size() * 3
        c_forces = (c_double * size)(0)

        # receive pressure from empire
        self.libempire_api.EMPIRE_API_recvSignal_double("force", size, c_forces)

        forces = Vector(3)
        i = 0
        # assign pressure to nodes of interface for current time step
        for nodes in (self.interface_model_part).Nodes:
            forces[0] = c_forces[3 * i]
            forces[1] = c_forces[3 * i + 1]
            forces[2] = c_forces[3 * i + 2]
            nodes.SetSolutionStepValue(FORCE, 0, forces);
            i = i + 1

        # Necessary to add forces to RS of structure
        AssignNoSlipCondition().AssignNoSlipCondition2D(self.model_part)
    # -------------------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------------------
    def sendDisplacements(self):
        # extract displacements from interface nodes as result of the CFD simulation
        displacements = []
        self.wrapper_process.ExtractDisplacementsFromModelPart(displacements)

        # convert list containg the displacements to ctypes
        c_displacements = (c_double * len(displacements))(*displacements)
        c_size = len(c_displacements)

        # send displacements to empire
        self.libempire_api.EMPIRE_API_sendSignal_double("displacements", c_size, c_displacements)
    # -------------------------------------------------------------------------------------------------

#

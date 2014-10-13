from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# import libraries
from KratosMultiphysics import *
from KratosMultiphysics.EmpireApplication import *
from KratosMultiphysics.IncompressibleFluidApplication import *
from ctypes import *
import os

CheckForPreviousImport()


class EmpireWrapper:

    # -------------------------------------------------------------------------------------------------
    def __init__(self, model_parts_list, libempire_api, dimension=3):
        self.libempire_api = libempire_api
        self.model_parts_list = model_parts_list
        N = len(self.model_parts_list)
        self.interface_model_parts_list = [ModelPart("InterfaceModelPart_"+str(i)) for i in range(N)]
        self.wrapper_processes_list = [WrapperProcess(self.model_parts_list[i], self.interface_model_parts_list[i], dimension) for i in range(N)]
    # -------------------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------------------
    # extract interface nodes and conditions from the fluid model part which have a flag "IS_INTERFACE"
    def extractInterface(self):
        for wrapper_process in self.wrapper_processes_list:
            wrapper_process.ExtractInterface()
    # -------------------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------------------
    # convert char* to python-compatible string
    def getUserDefinedText(self, stringText):
        empireFunc = self.libempire_api.EMPIRE_API_getUserDefinedText
        empireFunc.restype = c_char_p
        text = empireFunc(stringText)
        return text
    # -------------------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------------------
    def recvDisplacement(self):
        for interface_model_part in self.interface_model_parts_list:
            # initialize vector storing the displacements
            size = interface_model_part.GetNodes().Size() * 3
            c_displacements = (c_double * size)(0)
            
            # receive displacements from empire
            self.libempire_api.EMPIRE_API_recvDataField("defaultField", size, c_displacements)
            
            disp = Vector(3)
            i = 0
            # assign displacements to nodes of interface for current time step
            for node in (interface_model_part).Nodes:
                disp[0] = c_displacements[3 * i + 0]
                disp[1] = c_displacements[3 * i + 1]
                disp[2] = c_displacements[3 * i + 2]
                
                node.SetSolutionStepValue(DISPLACEMENT, 0, disp)
                
                i = i + 1
    # -------------------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------------------
    def sendPressure(self):
        for wrapper_process in self.wrapper_processes_list:
            # extract pressure from interface nodes as result of the CFD simulation
            pressure = []
            wrapper_process.ExtractPressureFromModelPart(pressure)
            
            # convert list containg the forces to ctypes
            c_pressure = (c_double * len(pressure))(*pressure)
            c_size = len(c_pressure)

            # send pressure to empire
            self.libempire_api.EMPIRE_API_sendDataField("defaultField", c_size, c_pressure)
    # -------------------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------------------
    def sendForces(self):
        for wrapper_process in self.wrapper_processes_list:
            # extract pressure from interface nodes as result of the CFD simulation
            forces = []
            wrapper_process.ExtractForcesFromModelPart(forces)

            # convert list containg the forces to ctypes
            c_forces = (c_double * len(forces))(*forces)
            c_size = len(c_forces)

            # send pressure to empire
            self.libempire_api.EMPIRE_API_sendDataField("defaultField", c_size, c_forces)
    # -------------------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------------------
    def sendMesh(self):
        for wrapper_process in self.wrapper_processes_list:
            numNodes = []
            numElems = []
            nodes = []
            nodeIDs = []
            numNodesPerElem = []
            elems = []

            # extract interface mesh information
            wrapper_process.ExtractMeshInfo(numNodes, numElems, nodes, nodeIDs, numNodesPerElem, elems)
            
            # Turek - Shell Element BEGIN
            #for j in range(1,len(nodes),3):
            #    nodes[j] = 0.20
            # Turek - Shell Element END

            # convert python lists to ctypes required for empire-function call
            c_numNodes = (c_int * len(numNodes))(*numNodes)
            c_numElems = (c_int * len(numElems))(*numElems)
            c_nodes = (c_double * len(nodes))(*nodes)
            c_nodeIDs = (c_int * len(nodeIDs))(*nodeIDs)
            c_numNodesPerElem = (c_int * len(numNodesPerElem))(*numNodesPerElem)
            c_elems = (c_int * len(elems))(*elems)

            # send mesh information to empire
            self.libempire_api.EMPIRE_API_sendMesh("defaultMesh", c_numNodes[0], c_numElems[0], c_nodes, c_nodeIDs, c_numNodesPerElem, c_elems)
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
        for interface_model_part in self.interface_model_parts_list:
            # initialize pressure variable
            size = interface_model_part.GetNodes().Size() * 1
            c_pressure = (c_double * size)(0)

            # receive pressure from empire
            self.libempire_api.EMPIRE_API_recvDataField("defaultField", size, c_pressure)
            
            pressure = Vector(1)
            i = 0
            # assign pressure to nodes of interface for current time step
            for nodes in (interface_model_part).Nodes:
                pressure[0] = c_pressure[i]
                # nodes.SetSolutionStepValue(POSITIVE_FACE_PRESSURE,0,pressure[0]);
                nodes.SetSolutionStepValue(NEGATIVE_FACE_PRESSURE, 0, pressure[0])
                i = i + 1
    # -------------------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------------------
    def recvForces(self):
        for interface_model_part in self.interface_model_parts_list:
            # initialize pressure variable
            size = interface_model_part.GetNodes().Size() * 3
            c_forces = (c_double * size)(0)
            
            # receive pressure from empire
            self.libempire_api.EMPIRE_API_recvDataField("defaultField", size, c_forces)
            
            forces = Vector(3)
            i = 0
            # assign pressure to nodes of interface for current time step
            for nodes in (interface_model_part).Nodes:
                forces[0] = c_forces[3 * i + 0]
                forces[1] = c_forces[3 * i + 1]
                forces[2] = c_forces[3 * i + 2]
                nodes.SetSolutionStepValue(FORCE, 0, forces)
                i = i + 1

        # Necessary to add forces to RS of structure
        for model_part in self.model_parts_list:
            AssignNoSlipCondition().AssignNoSlipCondition2D(model_part)
    # -------------------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------------------
    def sendDisplacements(self):
        for wrapper_process in self.wrapper_processes_list:
            # extract displacements from interface nodes as result of the CFD simulation
            displacements = []
            wrapper_process.ExtractDisplacementsFromModelPart(displacements)

            # convert list containg the displacements to ctypes
            c_displacements = (c_double * len(displacements))(*displacements)
            c_size = len(c_displacements)

            # send displacements to empire
            self.libempire_api.EMPIRE_API_sendDataField("defaultField", c_size, c_displacements)
    # -------------------------------------------------------------------------------------------------
#

from __future__ import print_function, absolute_import, division # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# import libraries
from KratosMultiphysics import *
from KratosMultiphysics.ExternalSolversApplication import *
# from KratosMultiphysics.IncompressibleFluidApplication import *
from KratosMultiphysics.SolidMechanicsApplication import *
from KratosMultiphysics.EmpireApplication import *
from ctypes import *
import os
from numpy import linalg as nla

CheckForPreviousImport()



class EmpireWrapper:
    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    def __init__(self, model_part, interface_mesh_number_list, interface_condition_number_list, libempire_api, is_2d = False):
        self.model_part = model_part
        self.interface_mesh_number_list = interface_mesh_number_list
        self.interface_condition_number_list = interface_condition_number_list
        self.libempire_api = libempire_api
        self.is_2d = is_2d
        
        # define interface model part...
        # ... via interface mesh number list...
        conditions = self.model_part.Conditions
        for interface_mesh_number in interface_mesh_number_list:
            interface_node_number_list = self.model_part.GetNodes(interface_mesh_number)
            for condition in conditions:
                condition_nodes = condition.GetNodes()
                is_interface_condition = True
                for condition_node in condition_nodes:
                    if condition_node.Id in interface_node_number_list:
                        is_interface_condition = is_interface_condition * True
                    else:
                        is_interface_condition = is_interface_condition * False
                if is_interface_condition:
                    condition.Set(INTERFACE, True)
                    for condition_node in condition_nodes:
                        condition_node.Set(INTERFACE, True)
        # ... via interface condition number list...
        for interface_condition_number in self.interface_condition_number_list:
            interface_condition = self.model_part.Conditions[interface_condition_number]
            interface_condition.Set(INTERFACE, True)
            for interface_node in interface_condition.GetNodes():
                interface_node.Set(INTERFACE,True)
            # ...
        self.interface_model_part = ModelPart("InterfaceModelPart")
        
        self.interface_model_part.AddNodalSolutionStepVariable(DISPLACEMENT)
        self.interface_model_part.AddNodalSolutionStepVariable(REACTION)
        self.interface_model_part.AddNodalSolutionStepVariable(POINT_LOAD)
        self.wrapper_process = WrapperProcess(self.model_part, self.interface_model_part, 3)
        print("defined interface model part and initialized EMPIRE wrapper")
    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    def connect(self, input_file_name):
        self.libempire_api.EMPIRE_API_Connect(input_file_name)
        print("connected to Emperor using input file *", input_file_name, "*")
    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    def sendMesh(self, extend_in_x = 0.0, extend_in_y = 0.0, extend_in_z = 1.0):
                
        # extract interface = extract conditions and nodes from model part if flag "INTERFACE" is true
        self.wrapper_process.ExtractInterface()
        self.number_of_interface_nodes = self.interface_model_part.GetNodes().Size()
        self.interface_nodes = (self.interface_model_part).Nodes
        
        # extract interface mesh information
        numNodes = []
        numElems = []
        nodeCoors = []
        nodeIDs = []
        numNodesPerElem = []
        elemTable = []
        self.wrapper_process.ExtractMeshInfo(numNodes, numElems, nodeCoors, nodeIDs, numNodesPerElem, elemTable)
        
        # extend mesh to 3d
        if self.is_2d:
            self.number_of_interface_nodes = self.number_of_interface_nodes * 2
            old_numNodes = numNodes[0]
            numNodes[0] = numNodes[0] * 2
            old_numElems = numElems[0]
            numElems[0] = numElems[0]
            old_nodeCoors = nodeCoors
            nodeCoors = nodeCoors + nodeCoors
            for i in range(old_numNodes):
                nodeCoors[                 3*i + 0] = old_nodeCoors[3*i + 0]
                nodeCoors[                 3*i + 1] = old_nodeCoors[3*i + 1]
                nodeCoors[                 3*i + 2] = old_nodeCoors[3*i + 2]
                nodeCoors[3*old_numNodes + 3*i + 0] = old_nodeCoors[3*i + 0] + extend_in_x
                nodeCoors[3*old_numNodes + 3*i + 1] = old_nodeCoors[3*i + 1] + extend_in_y
                nodeCoors[3*old_numNodes + 3*i + 2] = old_nodeCoors[3*i + 2] + extend_in_z
            max_nodeIDs = max(nodeIDs)
            nodeIDs = nodeIDs + nodeIDs
            for i in range(old_numNodes):
                nodeIDs[old_numNodes + i] = nodeIDs[old_numNodes + i] + max_nodeIDs
            numNodesPerElem = numNodesPerElem
            for i in range(old_numElems):
                numNodesPerElem[i] = 4
            old_elemTable = elemTable
            elemTable = elemTable + elemTable
            for i in range(old_numElems):
                elemTable[4*i + 0] = old_elemTable[2*i + 0]
                elemTable[4*i + 1] = old_elemTable[2*i + 1]
                elemTable[4*i + 2] = old_elemTable[2*i + 1] + max_nodeIDs
                elemTable[4*i + 3] = old_elemTable[2*i + 0] + max_nodeIDs
            print("extended mesh to 3d")
        
        # # for debug
        # print("extracted mesh info")
        # print(" - numNodes = ", numNodes)
        # print(" - numElems = ", numElems)
        # print(" - nodeCoors = ", nodeCoors)
        # print(" - nodeIDs = ", nodeIDs)
        # print(" - numNodesPerElem = ", numNodesPerElem)
        # print(" - elemTable = ",elemTable)
        
        # convert lists to ctypes
        c_nodeCoors = (c_double * (numNodes[0] * 3))(*nodeCoors)
        c_nodeIDs = (c_int * numNodes[0])(*nodeIDs)
        c_numNodesPerElem = (c_int * numElems[0])(*numNodesPerElem)
        c_elemTable = (c_int * len(elemTable))(*elemTable)
        
        # send mesh to Emperor
        self.libempire_api.EMPIRE_API_sendMesh("defaultMesh", numNodes[0], numElems[0], c_nodeCoors, c_nodeIDs, c_numNodesPerElem, c_elemTable)
        print("sent mesh to Emperor")
    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    def receiveDisplacements(self):
        # initialize displacement variable
        c_displacements = (c_double * (self.number_of_interface_nodes * 3))(0.0)
        
        # receive displacements from Emperor
        self.libempire_api.EMPIRE_API_recvDataField("defaultField", self.number_of_interface_nodes * 3, c_displacements)
        
        # reduce to 2d
        if self.is_2d:
            for i in range(int(self.number_of_interface_nodes * 3 / 2)):
                c_displacements[i] = (c_displacements[i] + c_displacements[int(self.number_of_interface_nodes * 3 / 2) + i]) / 2.0
        
        # # for debug
        # abs_displacement = 0.0
        # max_abs_displacement = 0.0
        # for displacement in c_displacements:
        #     abs_displacement = abs(displacement)
        #     if abs_displacement > max_abs_displacement:
        #         max_abs_displacement = abs_displacement
        # print("absolute maximum of received displacements = ", max_abs_displacement)
        
        # assign displacements to interface nodes
        displacement = Vector(3)
        i = 0
        for interface_node in self.interface_nodes:
            displacement[0] = c_displacements[3*i + 0]
            displacement[1] = c_displacements[3*i + 1]
            displacement[2] = c_displacements[3*i + 2]
            interface_node.SetSolutionStepValue(DISPLACEMENT, 0, displacement)
            i = i + 1
        
        # # for debug
        # abs_displacement = 0.0
        # max_abs_displacement = 0.0
        # for interface_node in self.interface_nodes:
        #     abs_displacement = abs(interface_node.GetSolutionStepValue(DISPLACEMENT, 0)[0])
        #     if abs_displacement > max_abs_displacement:
        #         max_abs_disp = abs_displacement
        #     abs_displacement = abs(interface_node.GetSolutionStepValue(DISPLACEMENT, 0)[1])
        #     if abs_displacement > max_abs_displacement:
        #         max_abs_disp = abs_displacement
        #     abs_displacement = abs(interface_node.GetSolutionStepValue(DISPLACEMENT, 0)[2])
        #     if abs_displacement > max_abs_displacement:
        #         max_abs_displacement = abs_displacement
        # print("absolute maximum of set displacements = ", max_abs_displacement)
        print("> > > received displacements from Emperor")
    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    def receiveForces(self):
        # initialize force variable
        c_forces = (c_double * (self.number_of_interface_nodes * 3))(0.0)
        
        # receive forces from Emperor
        self.libempire_api.EMPIRE_API_recvDataField("defaultField", self.number_of_interface_nodes * 3, c_forces)
        
        # reduce to 2d
        if self.is_2d:
            for i in range(int(self.number_of_interface_nodes * 3 / 2)):
                c_forces[i] = c_forces[i] + c_forces[int(self.number_of_interface_nodes * 3 / 2) + i]
        
        # # for debug
        int_force = Vector(3)
        int_force[0] = 0.0
        int_force[1] = 0.0
        int_force[2] = 0.0
        for i in range(self.number_of_interface_nodes):
            int_force[0] = int_force[0] + c_forces[3*i + 0]
            int_force[1] = int_force[1] + c_forces[3*i + 1]
            int_force[2] = int_force[2] + c_forces[3*i + 2]
        # if is_2d:
        #     int_force[0] = int_force[0] / 2.0
        #     int_force[1] = int_force[1] / 2.0
        #     int_force[2] = int_force[2] / 2.0
        print("integral of received forces = [", int_force[0], ", ", int_force[1], ", ", int_force[2], "]")
        
        # assign forces to interface nodes
        force = Vector(3)
        i = 0
        for interface_node in self.interface_nodes:
            force[0] = c_forces[3*i + 0]
            force[1] = c_forces[3*i + 1]
            force[2] = c_forces[3*i + 2]
            interface_node.SetSolutionStepValue(POINT_LOAD, 0, force)
            i = i + 1
        
        # # for debug
        int_force[0] = 0.0
        int_force[1] = 0.0
        int_force[2] = 0.0
        i = 0
        for node in (self.model_part).Nodes:
            int_force[0] = int_force[0] + node.GetSolutionStepValue(POINT_LOAD, 0)[0]
            int_force[1] = int_force[1] + node.GetSolutionStepValue(POINT_LOAD, 0)[1]
            int_force[2] = int_force[2] + node.GetSolutionStepValue(POINT_LOAD, 0)[2]
        print("integral of set point loads = [", int_force[0], ", ", int_force[1], ", ", int_force[2], "]")
        print("> > > received forces from Emperor")
    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    def receiveControlInput(self, control_input_node_number_list, direction = 1):
        # initialize control input variable
        c_control_input = (c_double * 1)(0.0)
        
        # receive control input from Emperor
        self.libempire_api.EMPIRE_API_recvSignal_double(b"controlInput", 1, c_control_input)
        
        # assign control input to nodes of control_input_node_number_list
        control_input = Vector(3)
        control_input[0] = 0.0
        control_input[1] = 0.0
        control_input[2] = 0.0
        control_input[direction] = c_control_input[0]
        for node_number in control_input_node_number_list:
            # Why not... "self.model_part.Nodes[nodeNum].SetSolutionStepValue(DISPLACEMENT, 0, control_input)"... but...
            self.model_part.Nodes[node_number].SetSolutionStepValue(DISPLACEMENT, 1, control_input) # ... ???
        
        # # for debug
        # print("received control input = ", c_control_input[0])
        # print("set Dirichlet BC = [", control_input[0] , ", ", control_input[1] , ", ", control_input[2] , "]")
        print("> > > received control input from Emperor")
    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    def sendForces(self):
        # # for debug
        # int_force = Vector(3)
        # int_force[0] = 0.0
        # int_force[1] = 0.0
        # int_force[2] = 0.0
        # for interface_node in self.interface_nodes:
        #     int_force[0] = int_force[0] + interface_node.GetSolutionStepValue(REACTION, 0)[0]
        #     int_force[1] = int_force[1] + interface_node.GetSolutionStepValue(REACTION, 0)[1]
        #     int_force[2] = int_force[2] + interface_node.GetSolutionStepValue(REACTION, 0)[2]
        # print("integral of calculated reactions = [", int_force[0], ", ", int_force[1], ", ", int_force[2], "]")
        
        # extract forces as simulation result from interface nodes
        forces = []
        self.wrapper_process.ExtractForcesFromModelPart(forces)
        
        # extend to 3d
        if self.is_2d:
            forces = forces + forces
            for i in range(self.number_of_interface_nodes * 3):
                forces[i] = forces[i] / 2.0
        
        # convert list to ctypes
        c_forces = (c_double * (self.number_of_interface_nodes * 3))(*forces)
        
        # # for debug
        # int_force[0] = 0.0
        # int_force[1] = 0.0
        # int_force[2] = 0.0
        # for i in range(self.number_of_interface_nodes):
        #     int_force[0] = int_force[0] + c_forces[3*i + 0]
        #     int_force[1] = int_force[1] + c_forces[3*i + 1]
        #     int_force[2] = int_force[2] + c_forces[3*i + 2]
        # if is_2d:
        #     int_force[0] = int_force[0] / 2.0
        #     int_force[1] = int_force[1] / 2.0
        #     int_force[2] = int_force[2] / 2.0
        # print("integral of sent forces = [", int_force[0], ", ", int_force[1], ", ", int_force[2], "]")
        
        # send forces to Emperor
        self.libempire_api.EMPIRE_API_sendDataField("defaultField", self.number_of_interface_nodes * 3, c_forces)
        print("< < < sent forces to Emperor")
    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    def sendDisplacements(self):
        # extract displacements as simulation result from interface nodes
        displacements = []
        self.wrapper_process.ExtractDisplacementsFromModelPart(displacements)
        
        # extend to 3d
        if self.is_2d:
            displacements = displacements + displacements
        
        # # for debug
        # abs_displacement = 0.0
        # max_abs_displacement = 0.0
        # for displacement in displacements:
        #     abs_displacement = abs(displacement)
        #     if abs_displacement > max_abs_displacement:
        #         max_abs_displacement = abs_displacement
        # print("absolute maximum of sent displacements = ", max_abs_displacement)
        
        # convert list to ctypes
        c_displacements = (c_double * (self.number_of_interface_nodes * 3))(*displacements)
        
        # send displacements to Emperor
        self.libempire_api.EMPIRE_API_sendDataField("defaultField", self.number_of_interface_nodes * 3, c_displacements)
        print("< < < sent displacements to Emperor")
    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    def sendMeasuredOutput(self, measured_output_node_number_list, direction = 1):
        # initialize variables
        measured_output = []
        size = 0
        displacement = Vector(3)
        displacement[0] = 0.0
        displacement[1] = 0.0
        displacement[2] = 0.0
        
        # extract measured output as simulation result from measured_output_node_number_list nodes
        for measured_output_node_number in measured_output_node_number_list:
            displacement = self.model_part.Nodes[measured_output_node_number].GetSolutionStepValue(DISPLACEMENT, 0)
            measured_output.append(displacement[direction])
            size = size + 1
        
        # convert list to ctypes
        c_measured_output = (c_double * size)(*measured_output)
        
        # send measured output to Emperor
        self.libempire_api.EMPIRE_API_sendSignal_double(b"measured_output", size, c_measured_output)
        print("< < < sent measured output to Emperor")
    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    def receiveConvergenceSignal(self):
        # receive convergence signal from Emperor
        is_converged = self.libempire_api.EMPIRE_API_recvConvergenceSignal()
        if is_converged:
            print("+ + + received convergence signal *true* from Emperor")
        else:
            print("- - - received convergence signal *false* from Emperor")
        return is_converged
    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    def disconnect(self):
        self.libempire_api.EMPIRE_API_Disconnect()
        print("disconnected from Emperor")
    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # convert char* to python-compatible string
    def getUserDefinedText(self, string_text):
        empire_func = self.libempire_api.EMPIRE_API_getUserDefinedText
        empire_func.restype = c_char_p
        text = empire_func(string_text)
        print("got user defined text *", text, "* from input file")
        return text
    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

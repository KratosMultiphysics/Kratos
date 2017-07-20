from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# import libraries
from KratosMultiphysics import *
from KratosMultiphysics.EmpireApplication import *
import empire_tools
import numpy as np
import ctypes as ctp
import os

CheckForPreviousImport()


class EmpireWrapper:
    # Wrapper for the EMPIRE API (/EMPIRE-Core/EMPIRE_API/src/include/EMPIRE_API.h)
    # -------------------------------------------------------------------------------------------------
    def __init__(self, model_part, libempire_api):
        self.libempire_api = libempire_api
        self.model_part = model_part
        self.tools = empire_tools.EmpireTools(self.model_part)
    # -------------------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------------------
    @staticmethod
    def Connect(libempire_api, xml_input_file):
        ''' Establishes the necessary connection with the Emperor '''
        libempire_api.EMPIRE_API_Connect(xml_input_file.encode())
    # -------------------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------------------
    @staticmethod
    def Disconnect(libempire_api):
        ''' Performs disconnection and finalization operations to the Emperor '''
        libempire_api.EMPIRE_API_Disconnect()
    # -------------------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------------------
    def GetUserDefinedText(self):
        ''' Get user defined text by the element name in the XML input file
        \param[in] elementName name of the XML element
        \return user defined text
        
        char *EMPIRE_API_getUserDefinedText(char *elementName); '''
        
        pass
    # -------------------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------------------
    def SendMesh(self):
        ''' Send the mesh to the server
        \param[in] name name of the mesh
        \param[in] numNodes number of nodes
        \param[in] numElems number of elements
        \param[in] nodes coordinates of all nodes
        \param[in] nodeIDs IDs of all nodes
        \param[in] numNodesPerElem number of nodes per element
        \param[in] elems connectivity table of all elements
         
        void EMPIRE_API_sendMesh(char *name, int numNodes, int numElems, double *nodes, int *nodeIDs,
                int *numNodesPerElem, int *elems);
        '''

        # extract interface mesh information
        numNodes = [];          numElems = []
        nodeCoors = [];         nodeIDs = []
        numNodesPerElem = [];   elemTable = []
        self.tools.ExtractInterface(numNodes, numElems, nodeCoors, nodeIDs, numNodesPerElem, elemTable)
        

        # convert python lists to ctypes, required for empire-function call
        c_numNodes = (ctp.c_int * len(numNodes))(*numNodes) # required?
        c_numElems = (ctp.c_int * len(numElems))(*numElems) # required?
        c_nodeCoors = (ctp.c_double * len(nodeCoors))(*nodeCoors)
        c_nodeIDs = (ctp.c_int * len(nodeIDs))(*nodeIDs)
        c_numNodesPerElem = (ctp.c_int * len(numNodesPerElem))(*numNodesPerElem)
        c_elemTable = (ctp.c_int * len(elemTable))(*elemTable)

        # send mesh to Emperor
        self.libempire_api.EMPIRE_API_sendMesh("defaultMesh", 
                                               c_numNodes[0], c_numElems[0], 
                                               c_nodeCoors, c_nodeIDs, 
                                               c_numNodesPerElem, c_elemTable)
    # -------------------------------------------------------------------------------------------------
    
    # -------------------------------------------------------------------------------------------------
    def ReceiveMesh(self):
        ''' Recieve mesh from the server
        \param[in] name name of the mesh
        \param[in] numNodes number of nodes
        \param[in] numElems number of elements
        \param[in] nodes coordinates of all nodes
        \param[in] nodeIDs IDs of all nodes
        \param[in] numNodesPerElem number of nodes per element
        \param[in] elems connectivity table of all elements
        
        void EMPIRE_API_recvMesh(char *name, int *numNodes, int *numElems, double **nodes, int **nodeIDs,
                int **numNodesPerElem, int **elem); '''

        # numNodes = []
        # numElems = []
        # nodeCoors = []
        # nodeIDs = []
        # numNodesPerElem = []
        # elemTable = []

        # c_numNodes = (ctp.c_int * len(numNodes))(*numNodes)
        # c_numElems = (ctp.c_int * len(numElems))(*numElems)
        # c_nodeCoors = (ctp.c_double * len(nodeCoors))(*nodeCoors)
        # c_nodeIDs = (ctp.c_int * len(nodeIDs))(*nodeIDs)
        # c_numNodesPerElem = (ctp.c_int * len(numNodesPerElem))(*numNodesPerElem)
        # c_elemTable = (ctp.c_int * len(elemTable))(*elemTable)

        # c_numNodes = ctp.POINTER(ctp.c_int)
        # c_numElems = ctp.POINTER(ctp.c_int)

        # ***********
        # numNodes = ctp.c_int(0)
        # numElems = ctp.c_int(0)
        # nodeCoors = ctp.c_double(0)
        # nodeIDs = ctp.c_int(0)
        # numNodesPerElem = ctp.c_int(0)
        # elemTable = ctp.c_int(0)

        c_numNodes = ctp.pointer(ctp.c_int(0))
        c_numElems = ctp.pointer(ctp.c_int(0))
        c_nodeCoors = ctp.pointer(ctp.pointer(ctp.c_double(0)))
        c_nodeIDs = ctp.pointer(ctp.pointer(ctp.c_int(0)))
        c_numNodesPerElem = ctp.pointer(ctp.pointer(ctp.c_int(0)))
        c_elemTable = ctp.pointer(ctp.pointer(ctp.c_int(0)))
        
        self.libempire_api.EMPIRE_API_recvMesh("defaultMesh", 
                                                c_numNodes, c_numElems, 
                                                c_nodeCoors, c_nodeIDs, 
                                                c_numNodesPerElem, c_elemTable)

        numNodes = c_numNodes.contents.value
        numElems = c_numElems.contents.value
        nodeCoors = c_nodeCoors.contents
        nodeIDs = c_nodeIDs.contents
        numNodesPerElem = c_numNodesPerElem.contents
        elemTable = c_elemTable.contents
        
        self.tools.ConstructMesh(numNodes, numElems, nodeCoors, nodeIDs, numNodesPerElem, elemTable)
    # -------------------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------------------
    def SendDataField(self, kratos_variable):
        ''' Send data field to the server
        \param[in] name name of the field
        \param[in] sizeOfArray size of the array (data field)
        \param[in] dataField the data field to be sent
        
        void EMPIRE_API_sendDataField(char *name, int sizeOfArray, double *dataField); '''

        # extract data field from nodes
        data_field = []
        self.tools.GetDataField(kratos_variable, data_field)

        # convert list containg the displacements to ctypes
        c_data_field = (ctp.c_double * len(data_field))(*data_field)
        c_size = len(c_data_field)

        # send displacements to empire
        self.libempire_api.EMPIRE_API_sendDataField("defaultField", c_size, c_data_field)
    # -------------------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------------------
    def ReceiveDataField(self, kratos_variable):   
        ''' Receive data field from the server
        \param[in] name name of the field
        \param[in] sizeOfArray size of the array (data field)
        \param[out] dataField the data field to be received
        
        void EMPIRE_API_recvDataField(char *name, int sizeOfArray, double *dataField); '''

        # Determine Size of Variable
        try:
            first_node = next(iter(self.model_part.Nodes))
            value = first_node.GetSolutionStepValue(kratos_variable)
            if (isinstance(value, float) or isinstance(value, int)): # Variable is a scalar
                size_of_variable = 1
            else:
                size_of_variable = len(first_node.GetSolutionStepValue(kratos_variable))
        except StopIteration:
            raise TypeError("size_of_variable could not be determined")

        # initialize vector storing the values
        size_data_field = self.model_part.NumberOfNodes() * size_of_variable
        c_data_field = (ctp.c_double * size_data_field)(0)
        
        # receive displacements from empire
        self.libempire_api.EMPIRE_API_recvDataField("defaultField", size_data_field, c_data_field)
        
        # if( (type(self.var) != KratosMultiphysics.VectorVariable) and (type(self.var) != KratosMultiphysics.Array1DVariable3) ):
        #     raise Exception("Variable type is incorrect. Must be a vector or a array_1d vector.")

        # if(type(self.var) != KratosMultiphysics.Array1DComponentVariable and type(self.var) != KratosMultiphysics.DoubleVariable and type(self.var) != KratosMultiphysics.VectorVariable):

        self.tools.SetDataField(kratos_variable, c_data_field, size_of_variable)
    # -------------------------------------------------------------------------------------------------
    
    
    # -------------------------------------------------------------------------------------------------
    def SendSignal(self):
        ''' Send signal to the server
        \param[in] name name of the signal
        \param[in] sizeOfArray size of the array (signal)
        \param[in] signal the signal
        
        void EMPIRE_API_sendSignal_double(char *name, int sizeOfArray, double *signal); '''

        pass
    # -------------------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------------------
    def ReceiveSignal(self):
        ''' Receive signal from the server
        \param[in] name name of the signal
        \param[in] sizeOfArray size of the array (signal)
        \param[in] signal the signal
        
        void EMPIRE_API_recvSignal_double(char *name, int sizeOfArray, double *signal); '''
        
        pass
    # -------------------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------------------
    def SendConvergenceSignal(self):
        ''' Send the convergence signal of an loop
        \param[in] signal 1 means convergence, 0 means non-convergence
        
        void EMPIRE_API_sendConvergenceSignal(int signal); '''

        pass
    # -------------------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------------------
    def ReceiveConvergenceSignal(self):
        '''Receive the convergence signal of an loop
        \return 1 means convergence, 0 means non-convergence
        
        int EMPIRE_API_recvConvergenceSignal(); '''

        pass
    # -------------------------------------------------------------------------------------------------







        # # receive dataFields from Emperor
        # self.client.recvDataField(dataFieldName_inEMPIRE)
        # numpy_dataField_vec=  self.client.dataField(dataFieldName_inEMPIRE)
        # # reduce to 2d
        # if self.is_2d:
        #     for i in range(int(self.number_of_interface_nodes * 3 / 2)):
        #         numpy_dataField_vec.array[i] = (numpy_dataField_vec.array[i] + numpy_dataField_vec.array[int(self.number_of_interface_nodes * 3 / 2) + i]) / 2.0
        
        # # # for debug
        # # abs_dataField = 0.0
        # # max_abs_dataField = 0.0
        # # for dataField in numpy_dataField_vec:
        # #     abs_dataField = abs(dataField)
        # #     if abs_dataField > max_abs_dataField:
        # #         max_abs_dataField = abs_dataField
        # # print("absolute maximum of received dataFields = ", max_abs_dataField)
        
        # # assign dataFields to interface nodes
        # dataField = Vector(3)
        # i = 0
        # for interface_node in self.model_part.Nodes:
        # #for interface_node in self.interface_nodes:
        #     dataField[0] = numpy_dataField_vec.array[3*i + 0]
        #     dataField[1] = numpy_dataField_vec.array[3*i + 1]
        #     dataField[2] = numpy_dataField_vec.array[3*i + 2]
        #     if(dataFieldName_inKRATOS == "FORCE"):
        #         interface_node.SetSolutionStepValue(POINT_LOAD, 0, dataField)
        #     if(dataFieldName_inKRATOS == "DISPLACEMENT"):
        #         interface_node.SetSolutionStepValue(DISPLACEMENT, 0, dataField)

        #     i = i + 1
        
        # # for debug
        # abs_dataField = 0.0
        # max_abs_dataField = 0.0
        # for interface_node in self.interface_nodes:
        #     abs_dataField = abs(interface_node.GetSolutionStepValue(dataField, 0)[0])
        #     if abs_dataField > max_abs_dataField:
        #         max_abs_disp = abs_dataField
        #     abs_dataField = abs(interface_node.GetSolutionStepValue(dataField, 0)[1])
        #     if abs_dataField > max_abs_dataField:
        #         max_abs_disp = abs_dataField
        #     abs_dataField = abs(interface_node.GetSolutionStepValue(DISPLACEMENT, 0)[2])
        #     if abs_dataField > max_abs_dataField:
        #         max_abs_dataField = abs_dataField
        # print("absolute maximum of set dataFields = ", max_abs_dataField)
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
        is_converged = self.client.recvConvergence()
        return is_converged
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






































    # -------------------------------------------------------------------------------------------------
    # extract interface nodes and conditions from the fluid model part which have a flag "IS_INTERFACE"
    def extractInterface(self):
        for wrapper_process in self.tools_processes_list:
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
        for wrapper_process in self.tools_processes_list:
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
        for wrapper_process in self.tools_processes_list:
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
        for wrapper_process in self.tools_processes_list:
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
        for wrapper_process in self.tools_processes_list:
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

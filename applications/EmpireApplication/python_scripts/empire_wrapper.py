from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# Other imports
import ctypes as ctp


class EmpireWrapper:
    # Source of Implementation: https://code.activestate.com/recipes/52558/
    # storage for the instance reference
    __instance = None

    def __init__(self, echo_level=0):
        """ Create singleton instance """
        # Check whether we already have an instance
        if EmpireWrapper.__instance is None:
            # Create and remember instance
            EmpireWrapper.__instance = EmpireWrapper.__EmpireWrapper(echo_level)

    def __getattr__(self, attr):
        """ Delegate access to implementation """
        return getattr(self.__instance, attr)

    def __setattr__(self, attr, value):
        """ Delegate access to implementation """
        return setattr(self.__instance, attr, value)

    class __EmpireWrapper:
        # Wrapper for the EMPIRE API (/EMPIRE-Core/EMPIRE_API/src/include/EMPIRE_API.h)
        # Implemented as Singleton, bcs otherwise the EMPIRE library can be imported several times
        ##### Constructor #####
        # -------------------------------------------------------------------------------------------------
        def __init__(self, echo_level):
            self.model_parts = {}
            self._load_empire_library()
            self.echo_level = echo_level
        # -------------------------------------------------------------------------------------------------

        ##### Public Functions #####
        # -------------------------------------------------------------------------------------------------
        def Connect(self, xml_input_file):
            ''' Establishes the necessary connection with the Emperor '''
            info_msg  = 'Attempting to connect with xml-file "'
            info_msg += xml_input_file + '"'
            KratosMultiphysics.Logger.PrintInfo('EMPIRE-Wrapper', info_msg)
            self.libempire_api.EMPIRE_API_Connect(xml_input_file.encode())
            KratosMultiphysics.Logger.PrintInfo('EMPIRE-Wrapper', 'Successfully connected')
        # -------------------------------------------------------------------------------------------------

        # -------------------------------------------------------------------------------------------------
        def Disconnect(self):
            ''' Performs disconnection and finalization operations to the Emperor '''
            KratosMultiphysics.Logger.PrintInfo('EMPIRE-Wrapper', 'Attempting to disconnect')
            self.libempire_api.EMPIRE_API_Disconnect()
            KratosMultiphysics.Logger.PrintInfo('EMPIRE-Wrapper', 'Successfully disconnected')
        # -------------------------------------------------------------------------------------------------

        # -------------------------------------------------------------------------------------------------
        def SendMesh(self, mesh_name, model_part):
            ''' Send the mesh to the server
            \param[in] name name of the mesh
            \param[in] numNodes number of nodes
            \param[in] numElems number of elements
            \param[in] nodes coordinates of all nodes
            \param[in] nodeIDs IDs of all nodes
            \param[in] numNodesPerElem number of nodes per element
            \param[in] elems connectivity table of all elements

            void EMPIRE_API_sendMesh(char *name, int numNodes, int numElems, double *nodes, int *nodeIDs,
                    int *numNodesPerElem, int *elems); '''
            # mesh_name: name of mesh in the emperor input

            # Save the ModelPart for data-field exchange later
            self._save_model_part(mesh_name, model_part)

            # extract interface mesh information
            numNodes = [];          numElems = []
            nodeCoors = [];         nodeIDs = []
            numNodesPerElem = [];   elemTable = []
            self._get_mesh(model_part, numNodes, numElems, nodeCoors, nodeIDs, numNodesPerElem, elemTable)

            # convert python lists to ctypes, required for empire-function call
            c_numNodes = (ctp.c_int * len(numNodes))(*numNodes)
            c_numElems = (ctp.c_int * len(numElems))(*numElems)
            c_nodeCoors = (ctp.c_double * len(nodeCoors))(*nodeCoors)
            c_nodeIDs = (ctp.c_int * len(nodeIDs))(*nodeIDs)
            c_numNodesPerElem = (ctp.c_int * len(numNodesPerElem))(*numNodesPerElem)
            c_elemTable = (ctp.c_int * len(elemTable))(*elemTable)

            if self.echo_level > 0:
                KratosMultiphysics.Logger.PrintInfo('EMPIRE-Wrapper', 'Attempting to send mesh "' + mesh_name + '"')
            # send mesh to Emperor
            self.libempire_api.EMPIRE_API_sendMesh("mesh_name.encode()",
                                                c_numNodes[0], c_numElems[0],
                                                c_nodeCoors, c_nodeIDs,
                                                c_numNodesPerElem, c_elemTable)
            if self.echo_level > 0:
                KratosMultiphysics.Logger.PrintInfo('EMPIRE-Wrapper', 'Sucessfully sent mesh')
        # -------------------------------------------------------------------------------------------------

        # -------------------------------------------------------------------------------------------------
        def ReceiveMesh(self, mesh_name, model_part):
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
            # mesh_name: name of mesh in the emperor input

            # Save the ModelPart for data-field exchange later
            self._save_model_part(mesh_name, model_part)

            c_numNodes = ctp.pointer(ctp.c_int(0))
            c_numElems = ctp.pointer(ctp.c_int(0))
            c_nodeCoors = ctp.pointer(ctp.pointer(ctp.c_double(0)))
            c_nodeIDs = ctp.pointer(ctp.pointer(ctp.c_int(0)))
            c_numNodesPerElem = ctp.pointer(ctp.pointer(ctp.c_int(0)))
            c_elemTable = ctp.pointer(ctp.pointer(ctp.c_int(0)))

            if self.echo_level > 0:
                KratosMultiphysics.Logger.PrintInfo('EMPIRE-Wrapper', 'Attempting to receive mesh "' + mesh_name + '"')
            self.libempire_api.EMPIRE_API_recvMesh(mesh_name.encode(),
                                                   c_numNodes, c_numElems,
                                                   c_nodeCoors, c_nodeIDs,
                                                   c_numNodesPerElem, c_elemTable)
            if self.echo_level > 0:
                KratosMultiphysics.Logger.PrintInfo('EMPIRE-Wrapper', 'Sucessfully received mesh')

            numNodes = c_numNodes.contents.value
            numElems = c_numElems.contents.value
            nodeCoors = c_nodeCoors.contents
            nodeIDs = c_nodeIDs.contents
            numNodesPerElem = c_numNodesPerElem.contents
            elemTable = c_elemTable.contents

            self._set_mesh(model_part, numNodes, numElems, nodeCoors, nodeIDs, numNodesPerElem, elemTable)
        # -------------------------------------------------------------------------------------------------

        # -------------------------------------------------------------------------------------------------
        def SendDataField(self, mesh_name, data_field_name, kratos_variables):
            ''' Send data field to the server
            \param[in] name name of the field
            \param[in] sizeOfArray size of the array (data field)
            \param[in] dataField the data field to be sent

            void EMPIRE_API_sendDataField(char *name, int sizeOfArray, double *dataField); '''
            # mesh_name: name of mesh in the emperor input
            # data_field_name: name of dataField in the emperor input

            # get ModelPart
            model_part = self.model_parts[mesh_name]

            if not type(kratos_variables) == list:
                kratos_variables = [kratos_variables]

            # extract data field from nodes
            data_field = self._get_data_field(model_part, kratos_variables)

            # convert list containg the data field to ctypes
            c_data_field = (ctp.c_double * len(data_field))(*data_field)
            c_size = len(c_data_field)

            if self.echo_level > 1:
                info_msg  = 'Attempting to send data-field "' + data_field_name
                info_msg += '" for mesh "' + mesh_name + '"\n                '
                info_msg += 'With Kratos-Vars: ' + ", ".join([var.Name() for var in kratos_variables])
                KratosMultiphysics.Logger.PrintInfo('EMPIRE-Wrapper', info_msg)
            # send data field to EMPIRE
            self.libempire_api.EMPIRE_API_sendDataField("data_field_name.encode()", c_size, c_data_field)
            if self.echo_level > 1:
                KratosMultiphysics.Logger.PrintInfo('EMPIRE-Wrapper', 'Successfully sent data-field')
        # -------------------------------------------------------------------------------------------------

        # -------------------------------------------------------------------------------------------------
        def ReceiveDataField(self, mesh_name, data_field_name, kratos_variables):
            ''' Receive data field from the server
            \param[in] name name of the field
            \param[in] sizeOfArray size of the array (data field)
            \param[out] dataField the data field to be received

            void EMPIRE_API_recvDataField(char *name, int sizeOfArray, double *dataField); '''
            # mesh_name: name of mesh in the emperor input
            # data_field_name: name of dataField in the emperor input

            # get ModelPart
            model_part = self.model_parts[mesh_name]

            if not type(kratos_variables) == list:
                kratos_variables = [kratos_variables]

            # Determine Sizes of Variables
            sizes_of_variables = self._sizes_of_variables(model_part, kratos_variables)
            self._check_size_of_variables(sizes_of_variables)

            # initialize vector storing the values
            size_data_field = model_part.NumberOfNodes() * sum(sizes_of_variables)
            c_size_data_field = ctp.c_int(size_data_field)
            c_data_field = (ctp.c_double * size_data_field)(0)

            if self.echo_level > 1:
                info_msg  = 'Attempting to receive data-field "' + data_field_name
                info_msg += '" for mesh "' + mesh_name + '"\n                '
                info_msg += 'With Kratos-Vars: ' + ", ".join([var.Name() for var in kratos_variables])
                KratosMultiphysics.Logger.PrintInfo('EMPIRE-Wrapper', info_msg)
            # receive data field from empire
            self.libempire_api.EMPIRE_API_recvDataField(data_field_name.encode(), c_size_data_field, c_data_field)
            if self.echo_level > 1:
                KratosMultiphysics.Logger.PrintInfo('EMPIRE-Wrapper', 'Successfully received data-field')

            self._set_data_field(model_part, kratos_variables, c_data_field, sizes_of_variables)
        # -------------------------------------------------------------------------------------------------

        # -------------------------------------------------------------------------------------------------
        def SendArray(self, array_name, array_to_send):
            ''' Send signal to the server
            \param[in] name name of the signal
            \param[in] sizeOfArray size of the array (signal)
            \param[in] signal the signal

            void EMPIRE_API_sendSignal_double(char *name, int sizeOfArray, double *signal); '''
            # array_name: name of signal in the emperor input

            # convert array to ctypes
            c_signal = (ctp.c_double * len(array_to_send))(*array_to_send)
            c_size = len(c_signal)

            if self.echo_level > 1:
                KratosMultiphysics.Logger.PrintInfo('EMPIRE-Wrapper', 'Attempting to send array "' + array_name + '"')
            self.libempire_api.EMPIRE_API_sendSignal_double(array_name.encode(), c_size, c_signal)
            if self.echo_level > 1:
                KratosMultiphysics.Logger.PrintInfo('EMPIRE-Wrapper', 'Successfully sent array')
        # -------------------------------------------------------------------------------------------------

        # -------------------------------------------------------------------------------------------------
        def ReceiveArray(self, array_name, array_size):
            ''' Receive signal from the server
            \param[in] name name of the signal
            \param[in] sizeOfArray size of the array (signal)
            \param[in] signal the signal

            void EMPIRE_API_recvSignal_double(char *name, int sizeOfArray, double *signal); '''
            # array_name: name of signal in the emperor input

            # initialize vector storing the values
            c_signal = (ctp.c_double * array_size)(0)

            if self.echo_level > 1:
                KratosMultiphysics.Logger.PrintInfo('EMPIRE-Wrapper', 'Attempting to receive array "' + array_name + '"')
            self.libempire_api.EMPIRE_API_recvSignal_double(array_name.encode(), array_size, c_signal)
            if self.echo_level > 1:
                KratosMultiphysics.Logger.PrintInfo('EMPIRE-Wrapper', 'Successfully received array')

            return self._convert_to_list(array_size, c_signal)
        # -------------------------------------------------------------------------------------------------

        # -------------------------------------------------------------------------------------------------
        def SendConvergenceSignal(self, signal):
            ''' Send the convergence signal of an loop
            \param[in] signal 1 means convergence, 0 means non-convergence
            void EMPIRE_API_sendConvergenceSignal(int signal); '''

            if self.echo_level > 1:
                KratosMultiphysics.Logger.PrintInfo('EMPIRE-Wrapper', 'Attempting to send convergence-signal')
            self.libempire_api.EMPIRE_API_sendConvergenceSignal(signal)
            if self.echo_level > 1:
                KratosMultiphysics.Logger.PrintInfo('EMPIRE-Wrapper', 'Successfully sent convergence-signal')
        # -------------------------------------------------------------------------------------------------

        # -------------------------------------------------------------------------------------------------
        def ReceiveConvergenceSignal(self):
            '''Receive the convergence signal of an loop
            \return 1 means convergence, 0 means non-convergence

            int EMPIRE_API_recvConvergenceSignal(); '''

            if self.echo_level > 1:
                KratosMultiphysics.Logger.PrintInfo('EMPIRE-Wrapper', 'Attempting to receive convergence-signal')
            signal = self.libempire_api.EMPIRE_API_recvConvergenceSignal()
            if self.echo_level > 1:
                KratosMultiphysics.Logger.PrintInfo('EMPIRE-Wrapper', 'Successfully received convergence-signal')
            return signal
        # -------------------------------------------------------------------------------------------------

        ##### Private Functions #####
        # -------------------------------------------------------------------------------------------------
        def _load_empire_library(self):
            if hasattr(self, 'libempire_api'): # the library has been loaded already
                raise ImportError("The EMPIRE library must be loaded only once!")

            import os

            if "EMPIRE_API_LIBSO_ON_MACHINE" not in os.environ:
                raise ImportError("The EMPIRE environment is not set!")

            KratosMultiphysics.Logger.PrintInfo('EMPIRE-Wrapper', 'Attempting to load empire-api')
            try: # OpenMPI
                self.libempire_api = ctp.CDLL(os.environ['EMPIRE_API_LIBSO_ON_MACHINE'], ctp.RTLD_GLOBAL)
                KratosMultiphysics.Logger.PrintInfo('EMPIRE-Wrapper', 'Using standard OpenMPI')
            except: # Intel MPI or OpenMPI compiled with "–disable-dlopen" option
                self.libempire_api = ctp.cdll.LoadLibrary(os.environ['EMPIRE_API_LIBSO_ON_MACHINE'])
                KratosMultiphysics.Logger.PrintInfo('EMPIRE-Wrapper', 'Using Intel MPI or OpenMPI compiled with "–disable-dlopen" option')
            KratosMultiphysics.Logger.PrintInfo('EMPIRE-Wrapper', 'Successfully loaded empire-api')
        # -------------------------------------------------------------------------------------------------

        # -------------------------------------------------------------------------------------------------
        def _get_mesh(self, model_part, num_nodes, num_elements, node_coords, node_IDs, num_nodes_per_element, element_table):
            num_nodes.append(model_part.NumberOfNodes())
            num_elements.append(model_part.NumberOfElements())

            for node in model_part.Nodes:
                node_coords.append(node.X)
                node_coords.append(node.Y)
                node_coords.append(node.Z)
                node_IDs.append(node.Id)

            for elem in model_part.Elements:
                num_nodes_per_element.append(len(elem.GetNodes()))
                for node in elem.GetNodes():
                    element_table.append(node.Id)
        # -------------------------------------------------------------------------------------------------

        # -------------------------------------------------------------------------------------------------
        def _set_mesh(self, model_part, num_nodes, num_elements, node_coords, node_IDs, num_nodes_per_element, element_table):
            # This function requires an empty ModelPart
            # It constructs Nodes and Elements from what was received from EMPIRE

            # Some checks to validate input:
            if model_part.NumberOfNodes() != 0:
                raise Exception("ModelPart is not empty, it has some Nodes!")
            if model_part.NumberOfElements() != 0:
                raise Exception("ModelPart is not empty, it has some Elements!")
            if model_part.NumberOfConditions() != 0:
                raise Exception("ModelPart is not empty, it has some Conditions!")

            # Create Nodes
            for i in range(num_nodes):
                model_part.CreateNewNode(node_IDs[i], node_coords[3*i+0], node_coords[3*i+1], node_coords[3*i+2]) # Id, X, Y, Z

            # Create dummy Property for Element
            prop = model_part.GetProperties()[1]

            element_table_counter = 0
            # Create Elements
            for i in range(num_elements):
                num_nodes_element = num_nodes_per_element[i]
                if num_nodes_element == 2:
                    name_element = "Element2D2N"
                elif num_nodes_element == 3:
                    name_element = "Element2D3N"
                elif num_nodes_element == 4:
                    name_element = "Element2D4N"
                else:
                    raise Exception("Wrong number of nodes for creating the element")

                element_nodes = []
                for j in range(num_nodes_element):
                    element_nodes.append(int(element_table[element_table_counter]))
                    element_table_counter += 1

                model_part.CreateNewElement(name_element, i+1, element_nodes, prop)
        # -------------------------------------------------------------------------------------------------

        # -------------------------------------------------------------------------------------------------
        def _get_data_field(self, model_part, kratos_variables):
            for var in kratos_variables:
                if not model_part.HasNodalSolutionStepVariable(var):
                    err_msg  = 'ModelPart "' + model_part.Name + '" does not have'
                    err_msg += ' "' + var.Name() + '" as SolutionStepVariable!'
                    raise Exception(err_msg)

            sizes_of_variables = self._sizes_of_variables(model_part, kratos_variables)
            self._check_size_of_variables(sizes_of_variables)

            num_nodes = model_part.NumberOfNodes()
            sum_sizes = sum(sizes_of_variables)

            data_field = [0.0] * (num_nodes * sum_sizes) # preallocate

            for node_i, node in enumerate(model_part.Nodes):
                size_index = 0
                for size_of_variable, variable in zip(sizes_of_variables, kratos_variables):
                    data_value = node.GetSolutionStepValue(variable)

                    if size_of_variable == 1:
                        data_field[node_i * sum_sizes + size_index] = data_value
                        size_index += 1
                    else:
                        for k in range(size_of_variable):
                            data_field[node_i * sum_sizes + size_index] = data_value[k]
                            size_index += 1

            return data_field
        # -------------------------------------------------------------------------------------------------

        # -------------------------------------------------------------------------------------------------
        def _set_data_field(self, model_part, kratos_variables, data_field, size_of_variables):
            # check if size of data field is correct
            if len(data_field) != model_part.NumberOfNodes() * sum(size_of_variables):
                raise Exception("received data field has wrong size!")

            for var in kratos_variables:
                if not model_part.HasNodalSolutionStepVariable(var):
                    err_msg  = 'ModelPart "' + model_part.Name + '" does not have'
                    err_msg += ' "' + var.Name() + '" as SolutionStepVariable!'
                    raise Exception(err_msg)

            # Preallocate values
            values = []
            for size_of_variable in size_of_variables:
                if size_of_variable > 1:
                    values.append(KratosMultiphysics.Vector(size_of_variable))
                else:
                    values.append(0.0)

            sum_sizes = sum(size_of_variables)

            # assign values to nodes of interface for current time step
            for node_i, node in enumerate(model_part.Nodes):
                size_index = 0
                for size_of_variable, variable, value in zip(size_of_variables, kratos_variables, values):
                    if size_of_variable == 1:
                        value = data_field[sum_sizes * node_i + size_index]
                        size_index += 1
                    else:
                        for k in range(size_of_variable):
                            value[k] = data_field[sum_sizes * node_i + size_index]
                            size_index += 1

                    node.SetSolutionStepValue(variable, 0, value)
        # -------------------------------------------------------------------------------------------------

        # -------------------------------------------------------------------------------------------------
        def _sizes_of_variables(self, model_part, kratos_variables):
            # this function is very general, even though EMPIRE works with Scalar and Vector quantities only!
            sizes_of_variables = []
            first_node = next(iter(model_part.Nodes))
            for variable in kratos_variables:
                try:
                    value = first_node.GetSolutionStepValue(variable)
                    if (isinstance(value, float) or isinstance(value, int)): # Variable is a scalar
                        size_of_variable = 1
                    else:
                        size_of_variable = len(value)
                    sizes_of_variables.append(size_of_variable)
                except StopIteration:
                    raise TypeError("Size of Variable \"" + variable + "\" could not be determined")

            return sizes_of_variables
        # -------------------------------------------------------------------------------------------------

        # -------------------------------------------------------------------------------------------------
        def _save_model_part(self, mesh_name, model_part):
            # Save the model_part for data-field exchange later
            if mesh_name in self.model_parts:
                raise ValueError("Mesh exsts already")
            else:
                self.model_parts.update({mesh_name : model_part})
        # -------------------------------------------------------------------------------------------------

        # -------------------------------------------------------------------------------------------------
        def _convert_to_list(self, array_size, c_signal):
            converted_list = [0.0] * array_size # preallocate

            for i in range(array_size):
                converted_list[i] = c_signal[i]

            return converted_list
        # -------------------------------------------------------------------------------------------------

        # -------------------------------------------------------------------------------------------------
        def _check_size_of_variables(self, size_of_variables):
            sum_sizes = sum(size_of_variables)

            possible_sizes = [ 1,  # Scalar
                               3,  # Vector
                               6 ] # doubleVector

            if sum_sizes not in possible_sizes:
                err_msg =  "Wrong size of variables: " + str(sum_sizes) + " !\n"
                err_msg += "Curently only Scalar, Vector and doubleVector are supported by Empire!"
                raise Exception(err_msg)
        # -------------------------------------------------------------------------------------------------

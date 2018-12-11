from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.MappingApplication as KratosMapping

# Importing the base class
from co_simulation_base_io import CoSimulationBaseIO

# Other imports
import numpy as np
from co_simulation_tools import csprint, bold, green, red

def Create(solvers, solver_name, level):
    return KratosIO(solvers, solver_name, level)

class KratosIO(CoSimulationBaseIO):
    def __init__(self, solvers, solver_name, level):
        super(KratosIO, self).__init__(solvers, solver_name, level)

        self.mappers = {}
        self.mapper_flags = {
            "add_values" : KratosMapping.Mapper.ADD_VALUES,
            "swap_sign" : KratosMapping.Mapper.SWAP_SIGN,
            "conservative" : KratosMapping.Mapper.USE_TRANSPOSE
        }

        self.kratos_vars = {} # dict storing name-KratosVars,
        # hopefully faster than accessing KratosComponents all the time


    def ImportData(self, data_settings, from_client):
        data_format = data_settings["data_format"]
        data_name = data_settings["data_name"]

        if data_format == "numpy_array":
            # TODO check if var in ModelPart!
            # In this case the from_client is the solver itself
            data_definition = from_client.GetDataDefinition(data_name)
            geometry_name = data_definition["geometry_name"]
            var_name = data_definition["data_identifier"]
            data_array = data_settings["data_array"]

            model_part = from_client.model[geometry_name]
            kratos_var = self.__GetKratosVariable(var_name)

            if type(kratos_var) == KratosMultiphysics.DoubleVariable or type(kratos_var) == KratosMultiphysics.Array1DComponentVariable:
                SetData(model_part, kratos_var, data_array)
            elif type(kratos_var) == KratosMultiphysics.Array1DVariable3:
                domain_size = model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]
                if not domain_size in [1,2,3]:
                    raise Exception("DOMAIN_SIZE has to be 1, 2 or 3!")
                num_nodes = NumberOfNodes(model_part)
                if data_array.size != num_nodes*domain_size:
                    raise Exception("Size of data does not match number of nodes x domain size!")
                ext = ["_X", "_Y", "_Z"]
                for i in range(domain_size):
                    component_var = self.__GetKratosVariable(kratos_var.Name()+ext[i])
                    range_begin = i*num_nodes
                    range_end = (i+1)*num_nodes
                    SetData(model_part, component_var, data_array[range_begin:range_end])
            else:
                err_msg  = 'Type of variable "' + kratos_var.Name() + '" is not valid\n'
                err_msg += 'It can only be double, component or array3d!'
                raise Exception(err_msg)

        elif data_format == "kratos_modelpart":
            to_client = self.solvers[self.solver_name]
            self.__Map(from_client, to_client, data_settings)

        elif data_format == "scalar_value":
            # TODO check if var in ModelPart!
            # In this case the from_client is the solver itself
            data_definition = from_client.GetDataDefinition(data_name)
            geometry_name = data_definition["geometry_name"]
            var_name = data_definition["data_identifier"]
            value = data_settings["scalar_value"]

            model_part = from_client.model[geometry_name]
            kratos_var = self.__GetKratosVariable(var_name)

            if type(kratos_var) == KratosMultiphysics.DoubleVariable or type(kratos_var) == KratosMultiphysics.Array1DComponentVariable:
                for node in Nodes(model_part):
                    node.SetSolutionStepValue(kratos_var, value)
                # KratosMultiphysics.VariableUtils().SetScalarVar(kratos_var, value, Nodes(model_part))
            else:
                err_msg  = 'Type of variable "' + kratos_var.Name() + '" is not valid\n'
                err_msg += 'It can only be double, component!'
                raise Exception(err_msg)

        else:
            raise Exception("The requested data_format is not implemented in KratosIO!")


    def ExportData(self, data_settings, to_client):
        data_format = data_settings["data_format"]
        data_name = data_settings["data_name"]

        if data_format == "numpy_array":
            # TODO check if var in ModelPart!
            # In this case the to_client is the solver itself
            data_definition = to_client.GetDataDefinition(data_name)
            geometry_name = data_definition["geometry_name"]
            var_name = data_definition["data_identifier"]
            data_array = data_settings["data_array"]
            buffer_index = data_settings["buffer_index"]

            model_part = to_client.model[geometry_name]
            kratos_var = self.__GetKratosVariable(var_name)

            if type(kratos_var) == KratosMultiphysics.DoubleVariable or type(kratos_var) == KratosMultiphysics.Array1DComponentVariable:
                required_size = NumberOfNodes(model_part)
                if not data_array.size == required_size:
                    # data_array = np.resize(data_array, (1,required_size))
                    data_array.resize(required_size, refcheck=False)
                ExtractData(model_part, kratos_var, data_array, buffer_index)
            elif type(kratos_var) == KratosMultiphysics.Array1DVariable3:
                domain_size = model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]
                if not domain_size in [1,2,3]:
                    raise Exception("DOMAIN_SIZE has to be 1, 2 or 3!")

                num_nodes = NumberOfNodes(model_part)
                required_size = num_nodes * domain_size
                if not data_array.size == required_size:
                    # data_array = np.resize(data_array, (1,required_size))
                    data_array.resize(required_size, refcheck=False)

                ext = ["_X", "_Y", "_Z"]
                for i in range(domain_size):
                    component_var = self.__GetKratosVariable(kratos_var.Name()+ext[i])
                    range_begin = i*num_nodes
                    range_end = (i+1)*num_nodes
                    ExtractData(model_part, component_var, data_array[range_begin:range_end], buffer_index)
            else:
                err_msg  = 'Type of variable "' + kratos_var.Name() + '" is not valid\n'
                err_msg += 'It can only be double, component or array3d!'
                raise Exception(err_msg)

        elif data_format == "kratos_modelpart":
            from_client = self.solvers[self.solver_name]
            self.__Map(from_client, to_client, data_settings)

        else:
            raise Exception("The requested data_format is not implemented in KratosIO!")

    def ExportMesh(self, data_settings, to_client):

        data_format = data_settings["data_format"]
        data_name = data_settings["data_name"]

        if data_format == "numpy_array_mesh":
            # In this case the to_client is the solver itself
            data_definition = to_client.GetDataDefinition(data_name)
            geometry_name = data_definition["geometry_name"]
            model_part = to_client.model[geometry_name]

            num_nodes = NumberOfNodes(model_part)
            num_elements = NumberOfElements(model_part)

            node_coords = np.zeros(num_nodes * 3)
            node_IDs = np.zeros(num_nodes, dtype=int)

            for i, node in enumerate(Nodes(model_part)):
                node_coords[i*3]    = node.X
                node_coords[i*3+1]  = node.Y
                node_coords[i*3+2]  = node.Z
                node_IDs[i] = node.Id

            num_nodes_per_element = []
            element_table = []

            for i, elem in enumerate(Elements(model_part)):
                num_nodes_per_element.append(len(elem.GetNodes()))
                for node in elem.GetNodes():
                    element_table.append(node.Id)

            num_nodes_per_element = np.array(num_nodes_per_element)
            element_table = np.array(element_table)

            numpy_mesh = {
                "num_nodes" : num_nodes,
                "num_elements" : num_elements,
                "node_coords" : node_coords,
                "node_IDs" : node_IDs,
                "num_nodes_per_element" : num_nodes_per_element,
                "element_table" : element_table
            }

            data_settings["mesh"] = numpy_mesh

        else:
            raise Exception("The requested data_format is not implemented in KratosIO!")

    def __GetMapper(self, from_client, to_client, data_settings):
        data_name = data_settings["data_name"]
        data_definition_from = from_client.GetDataDefinition(data_name)
        data_definition_to   = to_client.GetDataDefinition(data_name)

        geometry_name_from = data_definition_from["geometry_name"]
        geometry_name_to = data_definition_to["geometry_name"]

        mapper = None
        is_inverse_mapper = False
        if geometry_name_from in self.mappers: # a "Map"-Mapper exists
            if geometry_name_to in self.mappers[geometry_name_from]:
                mapper = self.mappers[geometry_name_from][geometry_name_to]

        if mapper == None and geometry_name_to in self.mappers: # an "InverseMap"-Mapper exists
            if geometry_name_from in self.mappers[geometry_name_to]:
                mapper = self.mappers[geometry_name_to][geometry_name_from]
                is_inverse_mapper = True

        if mapper == None: # no mapper for this geometry-pair exists, initializing a "Map"-Mapper
            if not geometry_name_from in self.mappers:
                self.mappers[geometry_name_from] = {}


            client_mesh_from = from_client.model[geometry_name_from]
            client_mesh_to = to_client.model[geometry_name_to]

            mapper_settings = KratosMultiphysics.Parameters("""{
                "mapper_type" : ""
            }""")
            mapper_settings["mapper_type"].SetString(data_settings["io_settings"]["mapper_type"])
            mapper = KratosMapping.MapperFactory.CreateMapper(client_mesh_from, client_mesh_to, mapper_settings)

            self.mappers[geometry_name_from][geometry_name_to] = mapper

            # Printing information related to mapping
            if self.echo_level > 2:
                info_msg  = bold("Mapper created") + ' for solver "' + self.solver_name + '": from "'
                info_msg += from_client._Name() + ':' + geometry_name_from + '" to "'
                info_msg += to_client._Name() + ':' + geometry_name_to + '"'
                csprint(self.lvl, info_msg)

        return mapper, is_inverse_mapper

    def __Map(self, from_client, to_client, data_settings):
            mapper, is_inverse_mapper = self.__GetMapper(from_client, to_client, data_settings)

            data_name = data_settings["data_name"]

            data_definition_from = from_client.GetDataDefinition(data_name)
            data_definition_to = to_client.GetDataDefinition(data_name)

            if is_inverse_mapper:
                var_origin = self.__GetKratosVariable(data_definition_to["data_identifier"])
                var_dest = self.__GetKratosVariable(data_definition_from["data_identifier"])
            else:
                var_origin = self.__GetKratosVariable(data_definition_from["data_identifier"])
                var_dest = self.__GetKratosVariable(data_definition_to["data_identifier"])

            var_origin_for_mapping = var_origin
            var_dest_for_mapping = var_dest

            mapper_flags = KratosMultiphysics.Flags()
            if "mapper_args" in data_settings["io_settings"]:
                for flag_name in data_settings["io_settings"]["mapper_args"]:
                    mapper_flags |= self.mapper_flags[flag_name]

            if "type_of_quantity" in data_definition_from:
                if data_definition_from["type_of_quantity"] == "nodal_point":
                    redistribution_tolerance = 1e-8
                    redistribution_max_iters = 50
                    geometry_name = data_definition_from["geometry_name"]
                    # Convert the nodal point quantities to distributed quantities before mapping
                    errr #Talk to Philipp befoe using this => VAUX_EQ_TRACTION has to be added to the ModelPart!
                    KratosMultiphysics.VariableRedistributionUtility.DistributePointValues(
                        from_client.model[geometry_name],
                        var_origin,
                        KratosMultiphysics.VAUX_EQ_TRACTION,
                        redistribution_tolerance,
                        redistribution_max_iters)
                    var_origin_for_mapping = KratosMultiphysics.VAUX_EQ_TRACTION
                    if self.echo_level > -1:
                        info_msg  = bold("Distributing Point Values of ")
                        info_msg += bold("Variable: ") + var_origin.Name()
                        info_msg += bold(" On: ") + geometry_name
                        csprint(self.lvl, info_msg)

            distribute_on_dest=False
            if "type_of_quantity" in data_definition_to:
                if data_definition_to["type_of_quantity"] == "nodal_point":
                    var_dest_for_mapping = KratosMultiphysics.VAUX_EQ_TRACTION
                    distribute_on_dest = True

            if is_inverse_mapper:
                mapper.InverseMap(var_origin_for_mapping, var_dest_for_mapping, mapper_flags)
            else:
                mapper.Map(var_origin_for_mapping, var_dest_for_mapping, mapper_flags)

            if distribute_on_dest:
                geometry_name = data_definition_to["geometry_name"]
                # Convert the transferred traction loads to point loads
                KratosMultiphysics.VariableRedistributionUtility.ConvertDistributedValuesToPoint(
                    to_client.model[geometry_name],
                    KratosMultiphysics.VAUX_EQ_TRACTION,
                    var_dest)
                if self.echo_level > -1:
                    info_msg  = bold("Converting Distributed-Values to Point-Values ")
                    info_msg += bold("Variable: ") + var_dest.Name()
                    info_msg += bold(" On: ") + geometry_name
                    csprint(self.lvl, info_msg)

            if self.echo_level > 3:
                pre_string = ""
                if is_inverse_mapper:
                    pre_string = "Inverse-"
                info_msg  = bold(pre_string+"Mapping with: ")
                info_msg += bold("Origin_Variable: ") + var_origin.Name() + " | "
                info_msg += bold("Destination_Variable: ") + var_dest.Name()
                if "mapper_args" in data_settings["io_settings"]:
                    info_msg += " | " + bold("Mapper-Flags: ") + ", ".join(data_settings["io_settings"]["mapper_args"])
                csprint(self.lvl, info_msg)

    def __GetKratosVariable(self, var_name):
        # TODO properly check which one is faster
        if not var_name in self.kratos_vars:
            self.kratos_vars[var_name] = KratosMultiphysics.KratosGlobals.GetVariable(var_name)

        return self.kratos_vars[var_name]

        # return KratosMultiphysics.KratosGlobals.GetVariable(var_name)


def ExtractData(model_part, kratos_var, data_array, buffer_index):
    for i, node in enumerate(Nodes(model_part)):
        data_array[i] = node.GetSolutionStepValue(kratos_var, buffer_index)

def SetData(model_part, kratos_var, data_array):
    num_nodes = NumberOfNodes(model_part)
    if data_array.size != num_nodes:
        raise Exception("Size of data does not match number of nodes!")
    for i, node in enumerate(Nodes(model_part)):
        node.SetSolutionStepValue(kratos_var, data_array[i])


def Nodes(model_part):
    # Wrapper to avoid long call
    return model_part.GetCommunicator().LocalMesh().Nodes

def NumberOfNodes(model_part):
    return len(Nodes(model_part)) # Mesh does currently not expose NumberOfNodes!

def Elements(model_part):
    # Wrapper to avoid long call
    return model_part.GetCommunicator().LocalMesh().Elements

def NumberOfElements(model_part):
    return len(Elements(model_part)) # Mesh does currently not expose NumberOfElements!

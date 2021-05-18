# Importing the Kratos Library
import KratosMultiphysics as KM

# Mapping imports
from KratosMultiphysics.MappingApplication import Mapper
from KratosMultiphysics.MappingApplication.python_mapper import PythonMapper

# other imports
import os
import ctypes as ctp

def Create(model_part_origin, model_part_destination, mapper_settings):
    return EmpireMortarMapper(model_part_origin, model_part_destination, mapper_settings)


class EmpireMortarMapper(PythonMapper):
    """Wrapper for the Mortar mapper of EMPIRE

    Usage:
    Empire needs to be compiled separately
    It can then be used by either of the two ways:
    - (default) use "startEMPIRE" to bring "EMPIRE_MAPPER_LIBSO_ON_MACHINE" to the environment
    - use "path_mapper_lib" to sprecify the path to "libEMPIRE_MapperLib.so" (by default located in "EMPIRE-Core/lib/")
    """

    mapper_count = 0
    instances = 0
    mapper_lib = None

    def __init__(self, model_part_origin, model_part_destination, mapper_settings):
        super().__init__(model_part_origin, model_part_destination, mapper_settings)

        if model_part_origin.IsDistributed() or model_part_destination.IsDistributed():
            raise Exception('{} does not support mapping with distributed ModelParts!'.format(self._ClassName()))

        if EmpireMortarMapper.mapper_lib:
            KM.Logger.PrintInfo("EmpireMortarMapper", "Mapper lib is already loaded")
        else:
            KM.Logger.PrintInfo("EmpireMortarMapper", "Attempting to load mapper lib")
            EmpireMortarMapper.mapper_lib = LoadEmpireMapperLib(self.mapper_settings["path_mapper_lib"].GetString())

        self.mapper_name = "EmpireMortarMapper_"+str(EmpireMortarMapper.mapper_count)

        self.mesh_name_origin = model_part_origin.FullName()+"_o_"+str(EmpireMortarMapper.mapper_count)
        self.mesh_name_destination = model_part_destination.FullName()+"_d_"+str(EmpireMortarMapper.mapper_count)

        self.__CreateEmpireFEMesh(self.model_part_origin, self.mesh_name_origin)
        self.__CreateEmpireFEMesh(self.model_part_destination, self.mesh_name_destination)

        self.__CreateMapper()

        EmpireMortarMapper.mapper_count += 1 # required for identification purposes
        EmpireMortarMapper.instances += 1
        self.__inverse_mapper = None

    def __del__(self):
        if EmpireMortarMapper.mapper_lib.hasMapper(ConvertToChar(self.mapper_name)):
            EmpireMortarMapper.mapper_lib.deleteMapper(ConvertToChar(self.mapper_name))

        if EmpireMortarMapper.mapper_lib.hasMesh(ConvertToChar(self.mesh_name_origin)):
            EmpireMortarMapper.mapper_lib.deleteMesh(ConvertToChar(self.mesh_name_origin))
        if EmpireMortarMapper.mapper_lib.hasMesh(ConvertToChar(self.mesh_name_destination)):
            EmpireMortarMapper.mapper_lib.deleteMesh(ConvertToChar(self.mesh_name_destination))

        EmpireMortarMapper.instances -= 1
        if EmpireMortarMapper.instances == 0: # last mapper was destoyed
            if self.echo_level > 1:
                KM.Logger.PrintInfo('EmpireMortarMapper', 'Destroying last instance, deleting all meshes & mappers')
            #  delete everything to make sure nothing is left
            EmpireMortarMapper.mapper_lib.deleteAllMappers()
            EmpireMortarMapper.mapper_lib.deleteAllMeshes()

    def UpdateInterface(self):
        raise NotImplementedError('"UpdateInterface" is not yet implemented for "{}"!'.format(self._ClassName()))

    # protected methods
    def _MapInternal(self, variable_origin, variable_destination, mapper_flags):
        if mapper_flags.Is(Mapper.USE_TRANSPOSE):
            mapper_flags.Reset(Mapper.USE_TRANSPOSE)
            mapper_flags.Set(KM.VISITED, True)
            self.__GetInverseMapper().Map(variable_destination, variable_origin, mapper_flags)
            return
        elif mapper_flags.Is(KM.VISITED):
            self.__MapInternalTranspose(variable_origin, variable_destination, mapper_flags)
            return

        self.__CheckMapperExists()

        var_dim = GetVariableDimension(variable_origin)

        origin_data_size = self.model_part_origin.NumberOfNodes()*var_dim
        destination_data_size = self.model_part_destination.NumberOfNodes()*var_dim

        c_origin_array = KratosFieldToCArray(self.model_part_origin.Nodes, variable_origin)
        c_destination_array = (ctp.c_double * destination_data_size)(0.0)

        EmpireMortarMapper.mapper_lib.doConsistentMapping(
            ConvertToChar(self.mapper_name),
            ctp.c_int(var_dim),
            ctp.c_int(origin_data_size),
            c_origin_array,
            ctp.c_int(destination_data_size),
            c_destination_array
            )

        CArrayToKratosField(
            c_destination_array,
            destination_data_size,
            self.model_part_destination.Nodes,
            variable_destination,
            mapper_flags.Is(Mapper.ADD_VALUES), mapper_flags.Is(Mapper.SWAP_SIGN))

    def _InverseMapInternal(self, variable_origin, variable_destination, mapper_flags):
        if mapper_flags.Is(Mapper.USE_TRANSPOSE):
            self.__MapInternalTranspose(variable_origin, variable_destination, mapper_flags)
        else:
            self.__GetInverseMapper().Map(variable_destination, variable_origin, mapper_flags)

    @classmethod
    def _GetDefaultParameters(cls):
        this_defaults = KM.Parameters("""{
            "path_mapper_lib"           : "",
            "dual"                      : false,
            "enforce_consistency"       : false,
            "opposite_normals"          : false,
            "use_initial_configuration" : false
        }""")
        this_defaults.AddMissingParameters(super()._GetDefaultParameters())
        return this_defaults

    # private methods
    def __CreateEmpireFEMesh(self, model_part, mesh_name):
        if model_part.NumberOfNodes() < 1:
            raise Exception('No nodes exist in ModelPart "{}"!'.format(model_part.FullName()))

        num_elements = model_part.NumberOfElements()
        num_conditions = model_part.NumberOfConditions()

        if num_elements > 0 and num_conditions > 0:
            err_msg  = "Both Elements and Conditions are present which is not allowed!\n"
            err_msg += "Name of ModelPart: {}\n".format(model_part.FullName())
            err_msg += "Number of Elements: {}\n".format(num_elements)
            err_msg += "Number of Conditions: {}".format(num_elements)
            raise Exception(err_msg)

        if num_elements + num_conditions == 0:
            err_msg  = "No Elements and Conditions are present which is not allowed!\n"
            err_msg += "Name of ModelPart: {}\n".format(model_part.FullName())
            raise Exception(err_msg)

        if num_conditions > 0:
            entities_to_use = model_part.Conditions
        else:
            entities_to_use = model_part.Elements

        for ent in entities_to_use:
            if ent.GetGeometry().PointsNumber() not in [3,4]:
                raise Exception("The EmpireMortarMapper only works with Triangles and Quadrilaterals")

        c_mesh_name = ConvertToChar(mesh_name)

        if EmpireMortarMapper.mapper_lib.hasMesh(c_mesh_name):
            raise Exception('Mesh "{}" exists already in Empire!'.format(mesh_name))

        c_num_nodes          = ctp.c_int(model_part.NumberOfNodes())
        c_num_elems          = ctp.c_int(len(entities_to_use))
        c_node_ids           = (ctp.c_int * c_num_nodes.value) (0)
        c_node_coords        = (ctp.c_double * (3 * c_num_nodes.value))(0.0)
        c_num_nodes_per_elem = (ctp.c_int * c_num_elems.value) (0)

        if self.mapper_settings["use_initial_configuration"].GetBool():
            for i_node, node in enumerate(model_part.Nodes):
                c_node_ids[i_node] = node.Id
                c_node_coords[i_node*3]   = node.X0
                c_node_coords[i_node*3+1] = node.Y0
                c_node_coords[i_node*3+2] = node.Z0

        else:
            for i_node, node in enumerate(model_part.Nodes):
                c_node_ids[i_node] = node.Id
                c_node_coords[i_node*3]   = node.X
                c_node_coords[i_node*3+1] = node.Y
                c_node_coords[i_node*3+2] = node.Z

        elem_node_ctr = 0
        for elem_ctr, elem in enumerate(entities_to_use):
            c_num_nodes_per_elem[elem_ctr] = len(elem.GetNodes())
            elem_node_ctr += c_num_nodes_per_elem[elem_ctr]

        elem_index = 0
        c_elems = (ctp.c_int * elem_node_ctr) (0)
        for elem_ctr, elem in enumerate(entities_to_use):
            for elem_node_ctr, elem_node in enumerate(elem.GetNodes()):
                c_elems[elem_index + elem_node_ctr] = elem_node.Id
            elem_index += len(elem.GetNodes())

        triangulateAll = False
        EmpireMortarMapper.mapper_lib.initFEMesh(c_mesh_name, c_num_nodes, c_num_elems, triangulateAll)
        EmpireMortarMapper.mapper_lib.setNodesToFEMesh(c_mesh_name, c_node_ids, c_node_coords)
        EmpireMortarMapper.mapper_lib.setElementsToFEMesh(c_mesh_name, c_num_nodes_per_elem, c_elems)

        if self.echo_level > 1:
            KM.Logger.PrintInfo('EmpireMortarMapper', 'Printing Mesh "{}"'.format(model_part.FullName()))
            EmpireMortarMapper.mapper_lib.printMesh(c_mesh_name)

    def __GetInverseMapper(self):
        if not self.__inverse_mapper:
            self.__inverse_mapper = self.__class__(self.model_part_destination, self.model_part_origin, self.mapper_settings)
        return self.__inverse_mapper

    def __CreateMapper(self):
        dual               = int(self.mapper_settings["dual"].GetBool())
        enforceConsistency = int(self.mapper_settings["enforce_consistency"].GetBool())
        opposite_normals   = int(self.mapper_settings["opposite_normals"].GetBool())

        EmpireMortarMapper.mapper_lib.initFEMMortarMapper(
            ConvertToChar(self.mapper_name),
            ConvertToChar(self.mesh_name_origin),
            ConvertToChar(self.mesh_name_destination),
            ctp.c_int(opposite_normals),
            ctp.c_int(dual),
            ctp.c_int(enforceConsistency)
        )

        self.__CheckMapperExists()

        EmpireMortarMapper.mapper_lib.buildCouplingMatrices(ConvertToChar(self.mapper_name))

    def __CheckMapperExists(self):
        if not EmpireMortarMapper.mapper_lib.hasMapper(ConvertToChar(self.mapper_name)):
            raise Exception('Mapper "{}" does not exist!'.format(self.mapper_name))

    def __MapInternalTranspose(self, variable_origin, variable_destination, mapper_flags):
        self.__CheckMapperExists()

        var_dim = GetVariableDimension(variable_destination)

        origin_data_size = self.model_part_origin.NumberOfNodes()*var_dim
        destination_data_size = self.model_part_destination.NumberOfNodes()*var_dim

        c_destination_array = KratosFieldToCArray(self.model_part_destination.Nodes, variable_destination)
        c_origin_array = (ctp.c_double * origin_data_size)(0.0)

        EmpireMortarMapper.mapper_lib.doConservativeMapping(
            ConvertToChar(self.mapper_name),
            ctp.c_int(var_dim),
            ctp.c_int(destination_data_size),
            c_destination_array,
            ctp.c_int(origin_data_size),
            c_origin_array
            )

        CArrayToKratosField(
            c_origin_array,
            origin_data_size,
            self.model_part_origin.Nodes,
            variable_origin,
            mapper_flags.Is(Mapper.ADD_VALUES), mapper_flags.Is(Mapper.SWAP_SIGN))

# Helper functions
def GetVariableDimension(variable):
    var_type = KM.KratosGlobals.GetVariableType(variable.Name())
    if var_type == "Array":
        return 3
    elif var_type == "Double":
        return 1
    else:
        raise Exception('Wrong variable type: "{}". Only "Array", "Double" and "Component" are allowed'.format(var_type))

def ConvertToChar(string):
    return ctp.c_char_p(string.encode(encoding='UTF-8'))

def LoadEmpireMapperLib(path_mapper_lib):
    KM.Logger.PrintInfo("EmpireMapperLibLoader", "Determining path to mapper lib")
    # first try automatic detection using the environment that is set by Empire => startEMPIRE
    if ('EMPIRE_MAPPER_LIBSO_ON_MACHINE' in os.environ):
        KM.Logger.PrintInfo("EmpireMapperLibLoader", "EMPIRE_MAPPER_LIBSO_ON_MACHINE found in environment")
        mapper_lib_path = os.environ['EMPIRE_MAPPER_LIBSO_ON_MACHINE']

    else:
        KM.Logger.PrintInfo("EmpireMapperLibLoader", "EMPIRE_MAPPER_LIBSO_ON_MACHINE NOT found in environment, using manually specified path to load the mapper lib")
        mapper_lib_path = path_mapper_lib
        if mapper_lib_path == "":
            raise Exception('The automatic detection of the mapper lib failed, the path to the mapper lib has to be specified!')

    KM.Logger.PrintInfo("EmpireMapperLibLoader", "Attempting to load the mapper lib")
    # TODO check if still both are needed! (does the mapperlib link to MPI?)
    try:
        try: # OpenMPI
            loaded_mapper_lib = ctp.CDLL(mapper_lib_path, ctp.RTLD_GLOBAL)
            KM.Logger.PrintInfo('EmpireMapperLibLoader', 'Using standard OpenMPI')
        except: # Intel MPI or OpenMPI compiled with "–disable-dlopen" option
            loaded_mapper_lib = ctp.cdll.LoadLibrary(mapper_lib_path)
            KM.Logger.PrintInfo('EmpireMapperLibLoader', 'Using Intel MPI or OpenMPI compiled with "–disable-dlopen" option')
    except OSError:
        raise Exception("Mapper lib could not be loaded!")

    KM.Logger.PrintInfo("EmpireMapperLibLoader", 'Successfully loaded the mapper lib from "{}"'.format(mapper_lib_path))

    return loaded_mapper_lib


def SetSolutionStepValue(node, variable, value):
    return node.SetSolutionStepValue(variable, 0, value)

def KratosFieldToCArray(nodes, variable):
    dim = GetVariableDimension(variable)
    size = dim * len(nodes)
    c_array = (ctp.c_double * size)(0.0)

    if dim == 1:
        for i_node, node in enumerate(nodes):
            c_array[i_node] = node.GetSolutionStepValue(variable)
    else:
        for i_node, node in enumerate(nodes):
            node_value = node.GetSolutionStepValue(variable)
            for i_dim in range(dim):
                c_array[i_node*dim + i_dim] = node_value[i_dim]

    return c_array

def CArrayToKratosField(c_array, c_array_size, nodes, variable, add_values, swap_sign):
    dim = GetVariableDimension(variable)
    if c_array_size != dim * len(nodes):
        raise RuntimeError("Wrong size!")

    if swap_sign:
        for i in range(c_array_size):
            c_array[i] *= (-1)

    if add_values:
        current_values = KratosFieldToCArray(nodes, variable)
        for i in range(c_array_size):
            c_array[i] += current_values[i]

    if dim == 1:
        for i_node, node in enumerate(nodes):
            node.SetSolutionStepValue(variable, 0, c_array[i_node])
    else:
        for i_node, node in enumerate(nodes):
            values = [c_array[i_node*dim], c_array[i_node*dim+1], c_array[i_node*dim+2]]
            node.SetSolutionStepValue(variable, 0, values)

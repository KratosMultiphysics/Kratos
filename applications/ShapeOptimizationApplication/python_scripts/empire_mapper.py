import ctypes as ctp
import os

class EmpireMapper():
    def __init__(self, origin_model_part, destination_model_part, settings):
        try: # OpenMPI
            self.mapper_library = ctp.CDLL(os.environ['EMPIRE_MAPPER_LIBSO_ON_MACHINE'], ctp.RTLD_GLOBAL)
            print("::EMPIRE:: Using standard OpenMPI")
        except: # Intel MPI or OpenMPI compiled with "–disable-dlopen" option
            self.mapper_library = ctp.cdll.LoadLibrary(os.environ['EMPIRE_MAPPER_LIBSO_ON_MACHINE'])
            print("::EMPIRE:: Using Intel MPI or OpenMPI compiled with \"–disable-dlopen\" option")

        self.name = "Kratos_Empire_Mapper_Base".encode(encoding='UTF-8')

        self.origin_model_part = origin_model_part
        self.origin_model_part_name = origin_model_part.Name.encode(encoding='UTF-8')
        self.destination_model_part = destination_model_part
        if destination_model_part.Name == origin_model_part.Name:
            self.destination_model_part_name = (destination_model_part.Name + "_destination").encode(encoding='UTF-8')
        else:            
            self.destination_model_part_name = destination_model_part.Name.encode(encoding='UTF-8')

        self.__MakeEmpireFEMesh(self.origin_model_part_name, origin_model_part)
        self.__MakeEmpireFEMesh(self.destination_model_part_name, destination_model_part)

        self.consistent_mapping = settings["consistent_mapping"].GetBool()

    def __del__(self):
        print("Base Class Delete")
        
        if bool(self.mapper_library.hasMapper(ctp.c_char_p(self.name))):
            self.mapper_library.deleteMapper(ctp.c_char_p(self.name))

        if bool(self.mapper_library.hasMesh(ctp.c_char_p(self.origin_model_part_name))):
            self.mapper_library.deleteMesh(ctp.c_char_p(self.origin_model_part_name))

        if bool(self.mapper_library.hasMesh(ctp.c_char_p(self.destination_model_part_name))):
            self.mapper_library.deleteMesh(ctp.c_char_p(self.destination_model_part_name))

    def Initialize(self):

        # if bool(self.mapper_library.hasMesh(ctp.c_char_p(self.origin_model_part_name))):
        #     self.mapper_library.printMesh(ctp.c_char_p(self.origin_model_part_name))
        # if bool(self.mapper_library.hasMesh(ctp.c_char_p(self.destination_model_part_name))):
        #     self.mapper_library.printMesh(ctp.c_char_p(self.destination_model_part_name))

        if self.mapper_library.hasMapper(ctp.c_char_p(self.name)):
            self.mapper_library.buildCouplingMatrices(ctp.c_char_p(self.name))

    def Update(self):
        pass

    def Map(self, origin_variable, destination_variable):
        print("EmpireVertexMorphingMapper.Map")
        dim = 3

        origin_data_size = len(self.origin_model_part.Nodes)*dim
        destination_data_size = len(self.destination_model_part.Nodes)*dim

        c_origin_array = self.__KratosFieldToCArray(
                                      dim, self.origin_model_part, origin_variable)
        c_destination_array = (ctp.c_double * destination_data_size)(0.0)

        self.mapper_library.doConsistentMapping(
            ctp.c_char_p(self.name), 
            ctp.c_int(dim), 
            ctp.c_int(origin_data_size),
            c_origin_array,
            ctp.c_int(destination_data_size),
            c_destination_array
            )

        self.__CArrayToKratosField(
            dim, 
            c_destination_array, 
            destination_data_size, 
            self.destination_model_part, 
            destination_variable)
    
    def InverseMap(self, destination_variable, origin_variable):
        print("EmpireVertexMorphingMapper.InverseMap")
        if (self.consistent_mapping):
            print("Calling Map")
            self.Map(destination_variable, origin_variable)
            return
        dim = 3

        destination_data_size = len(self.destination_model_part.Nodes)*dim
        origin_data_size = len(self.origin_model_part.Nodes)*dim
        
        c_destination_array = self.__KratosFieldToCArray(
            dim, 
            self.destination_model_part, 
            destination_variable)
        c_origin_array = (ctp.c_double * origin_data_size)(0.0)
        
        self.mapper_library.doConservativeMapping(
            ctp.c_char_p(self.name), 
            ctp.c_int(dim),
            ctp.c_int(destination_data_size),
            c_destination_array,
            ctp.c_int(origin_data_size),
            c_origin_array
            )
        
        self.__CArrayToKratosField(
            dim, 
            c_origin_array, 
            origin_data_size, 
            self.origin_model_part, 
            origin_variable)

    def __MakeEmpireFEMesh(self, mesh_name, model_part):

        c_mesh_name       = ctp.c_char_p(mesh_name)
        c_num_nodes       = ctp.c_int(len(model_part.Nodes))
        c_num_elems       = ctp.c_int(len(model_part.Conditions))
        c_node_ids        = (ctp.c_int * c_num_nodes.value) (0)
        c_node_coords     = (ctp.c_double * (3 * c_num_nodes.value))(0.0)
        c_num_nodes_per_elem = (ctp.c_int * c_num_elems.value) (0)

        for node_ctr, node in enumerate(model_part.Nodes):
            c_node_ids[node_ctr] = node.Id
            c_node_coords[node_ctr*3]   = node.X
            c_node_coords[node_ctr*3+1] = node.Y
            c_node_coords[node_ctr*3+2] = node.Z

        elem_node_ctr = 0
        for elem_ctr, elem in enumerate(model_part.Conditions):
            c_num_nodes_per_elem[elem_ctr] = len(elem.GetNodes())
            elem_node_ctr += c_num_nodes_per_elem[elem_ctr]

        elem_index = 0
        c_elems = (ctp.c_int * elem_node_ctr) (0)
        for elem_ctr, elem in enumerate(model_part.Conditions):
            for elem_node_ctr, elem_node in enumerate(elem.GetNodes()):
                c_elems[elem_index + elem_node_ctr] = elem_node.Id
            elem_index += len(elem.GetNodes())
        
        self.mapper_library.initFEMesh(c_mesh_name, c_num_nodes, c_num_elems, False)
        self.mapper_library.setNodesToFEMesh(c_mesh_name, c_node_ids, c_node_coords)
        self.mapper_library.setElementsToFEMesh(c_mesh_name, c_num_nodes_per_elem, c_elems)

    def __KratosFieldToCArray(self, dim, model_part, variable):

        #TODO check the dimension of the field before
        size = dim * len(model_part.Nodes)

        c_array = (ctp.c_double * size)(0.0)
        for node_ctr, node in enumerate(model_part.Nodes):
            node_value = node.GetSolutionStepValue(variable)
            # print(node_value)
            if node_value.Size() != dim:
                raise RuntimeError("Wrong dimensions!")
            for iDim in range(dim):
                c_array[node_ctr*dim + iDim] = node_value[iDim]

        return c_array

    def __CArrayToKratosField(self, dim, c_array, size, model_part, variable):

        #TODO check the dimension of the field before
        if size != dim * len(model_part.Nodes):
            raise RuntimeError("Wrong size!")
        
        # TODO
        for node_ctr, node in enumerate(model_part.Nodes):
            values = [ c_array[node_ctr*dim], c_array[node_ctr*dim+1], c_array[node_ctr*dim+2]]
            node.SetSolutionStepValue(variable, values)

class EmpireVertexMorphingMapper(EmpireMapper):
    def __init__(self, origin_model_part, destination_model_part, settings):
        super(EmpireVertexMorphingMapper, self).__init__(origin_model_part, destination_model_part, settings)
        
        self.name = "Kratos_Empire_Mapper_VertexMorphing".encode(encoding='UTF-8')
        
        self.mapper_library.initVertexMorphingMapper(
            ctp.c_char_p(self.name),
            ctp.c_char_p(self.origin_model_part_name), 
            ctp.c_char_p(self.destination_model_part_name))

        self.filter_types = {
            "linear": 0,
            "gaussian": 1
        }

        filter_type_string = settings["filter_function_type"].GetString()
        filter_type_int = self.filter_types[filter_type_string]
        filter_radius = settings["filter_radius"].GetDouble()

        self.mapper_library.setVMParameters(
            ctp.c_char_p(self.name),
            ctp.c_int(filter_type_int),
            ctp.c_double(filter_radius))
    

class EmpireFEMortarMapper(EmpireMapper):
    def __init__(self, origin_model_part, destination_model_part, settings):
        super(EmpireFEMortarMapper, self).__init__(origin_model_part, destination_model_part, settings)
        
        self.name = "Kratos_Empire_Mapper_FEMortar".encode(encoding='UTF-8')
        
        opposite_normal = int(settings["empire_settings"]["opposite_normal"].GetBool())
        dual = int(False)
        enforceConsistency = int(settings["empire_settings"]["enforce_consistency"].GetBool())
        
        self.mapper_library.initFEMMortarMapper(
            ctp.c_char_p(self.name), 
            ctp.c_char_p(self.origin_model_part_name),
            ctp.c_char_p(self.destination_model_part_name),
            ctp.c_int(opposite_normal),
            ctp.c_int(dual),
            ctp.c_int(enforceConsistency)
            )

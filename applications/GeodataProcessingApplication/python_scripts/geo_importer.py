import KratosMultiphysics
from geo_processor import GeoProcessor

class GeoImporter( GeoProcessor ):

    ### --- overwriting the import function --- ###

    def SetGeoModelPart( self, modelPartIn ):

        KratosMultiphysics.Logger.PrintWarning("GeoImporter: SetGeoModelPart", "This function cannot be used. The model is created be reading a file.")


    """ FUNCTION UNDER CONSTRUCTION """
    def ObjToSplit(self, obj_file_name_input, name_model_part="ModelPart"):
        # with this function, we split the obj file
        # (is useful if we have a file with all buildings in the same group)

        self._InitializeModelPart(name_model_part)

        # import obj file and fill a ModelPart with nodes and elements
        with open (obj_file_name_input) as read_file:
            ID_vertex = 1
            ID_elem = 1

            for row in read_file.readlines():
                row = row.split()
                if (row[0] == "v"):
                    x = float(row[1])
                    y = float(row[2])
                    z = float(row[3])
                    self.ModelPart.CreateNewNode(ID_vertex, x, y, z)
                    ID_vertex += 1
                
                elif (row[0] == "f"):
                    n1 = int(row[1])
                    n2 = int(row[2])
                    n3 = int(row[3])
                    self.ModelPart.CreateNewElement("Element2D3N", ID_elem, [n1, n2, n3], self.ModelPart.GetProperties()[2])
                    ID_elem += 1
        
        for node in self.ModelPart.Nodes:
            node.Set(KratosMultiphysics.TO_ERASE,True)      # we set all nodes as TO_ERASE
            # node.Set(KratosMultiphysics.VISITED,False)      # we set all nodes as not VISITED
        
        # we mark all free nodes as TO_ERASE
        for elem in self.ModelPart.Elements:
            for node in elem.GetNodes():
                node.Set(KratosMultiphysics.TO_ERASE, False)
        
        # we erase unused nodes
        self.ModelPart.RemoveNodesFromAllLevels(KratosMultiphysics.TO_ERASE)

        self.HasModelPart = True

        #########################################################################################################
        # geom_id = 1         # a counter
        # dict_visited = []	# list with the elements in the dictionary already visited

        # num_elem = self.ModelPart.NumberOfElements()    # key_list = list(elem_map.keys())	# it is necessary because we delete elements in the dictionary; so we can't iterate on dictionary directly
        # """for faces in key_list:"""
        # for faces in self.ModelPart.Elements:
        #     if (faces.Id in dict_visited):
        #         continue
        #     dict_visited.append(faces.Id)

        #     next_faces = [faces.Id]	# the next faces that will be processed
        #     input(next_faces)

        #     node_visited = []	# list with nodes already visited in this geometry
        #     elem_visited = []	# list with elements already visited in this geometry
        #     while (next_faces):	# this loop ends when next_faces is empty
        #         current_face = next_faces[0]			# get the first element in next_faces up to the last element. In next operations we subtract elements in next_faces
        #         """for node in elem_map[current_face]:"""
        #         for node in self.ModelPart.GetElement(current_face):        # CHECK THIS IF IT IS CORRECT
        #             elem_visited.append(current_face)	# we update elem_visited
        #             dict_visited.append(current_face)
        #             if (node.Id in node_visited):
        #                 continue						# we go on if the node is already visited
        #             node_visited.append(node.Id)			# we update node_visited
        #             """for id, elem in elem_map.items():"""
        #             for elem in self.ModelPart.Elements:
        #                 if (node in elem):
        #                     next_faces.append(id)
        #         next_faces = list(set(next_faces) - set(elem_visited))	# we remove duplicates and the elements that are already visited

        #         del(elem_map[current_face])			# we delete current_face in elem_map

        #     geometries[geom_id] = list(set(elem_visited))			# key: id geometry; value: id elements of this geometry
        #     geom_id += 1

        dict_Building = {}      # key: Building ids; value: element ids
        dict_Node = {}          # key: Building ids; value: node ids
        Building_id = 0         # we initialize the Building id. It is a progressive number
        list_visited = []       # node visited
        #######################################################################################################################

        for elem in self.ModelPart.Elements:
            current_id = 0  # temporary key (useful to store nodes and elements in dictionaries)
            temp_list = []

            nnn = elem.GetNodes()
            print("\nelem.Id: {}; node: {}".format(elem.Id, [nnn[0].Id, nnn[1].Id, nnn[2].Id]))
            print("\ttemp_list before: ", temp_list)

            for node in elem.GetNodes():
                print("\t\tnode: ", node.Id)
                if not bool(dict_Node):
                    # dict_Node is empty
                    temp_list.append(node.Id)
                    print("\t\tdict_Node was empty")
                else:
                    new_node = True     # True only if the Node is not in the dictionary
                    for i_building, list_nodes in dict_Node.items():
                        if not (node.Id in list_nodes):
                            # node is not in the dict
                            # # temp_list.append(node.Id)
                            new_node = True
                            print("\t\tNode {} is not in Building {}".format(node.Id, i_building))
                        else:
                            # node is in dict
                            current_id = i_building
                            print("\t\tNode {} is in dict at Building: {}".format(node.Id, i_building))
                            new_node = False
                            break
                    
                    if new_node:
                        # we add node.Id in temp_list only if the node is not in the entire dictionary
                        temp_list.append(node.Id)

            


            print("\ttemp_list: ", temp_list)

            if (len(temp_list) == 3):
                # in this case all nodes are in a different Building
                Building_id += 1                                # id of the new Building
                dict_Building[Building_id] = [elem.Id]
                dict_Node[Building_id] = temp_list
            else:
                # if at least one node is in dict_Building, in current_id there is a Building id
                dict_Building[current_id].append(elem.Id)       # we append the Element id in the dictionary
                dict_Node[current_id].extend(temp_list)         # we extend the Node ids in the dictionary
            
            print("\tdict_Node: {}".format(dict_Node))

            
            # # we add "temp_list" into "list_visited"
            # list_visited.extend(temp_list)

            # if new_building:
            #     Building_id += 1
            #     dict_Building[Building_id] = [elem.Id]
            # else:
            #     # we append the Element id in the dictionary
            #     dict_Building[Building_id].append(elem.Id)

        #######################################################################################################################

        current_sub_model_building = self.ModelPart
        for building_id, elem_id in dict_Building.items():
            # we check if the sub model part already exists
            if not self.ModelPart.HasSubModelPart("Building_{}".format(building_id)):
                current_sub_model_building = self.ModelPart.CreateSubModelPart("Building_{}".format(building_id))
            
            list_nodes = []
            for id in elem_id:
                elem = self.ModelPart.GetElement(id)
                for node in elem.GetNodes():
                    list_nodes.append(node.Id)
            
            # we add the nodes in the sub model part
            current_sub_model_building.AddNodes(list(set(list_nodes)))

            # we add the element in the sub model part
            current_sub_model_building.AddElements(list(elem_id))


        #########################################################################################################


    def ObjImportBuildings(self, obj_file_name_input, name_model_part="ModelPart", change_coord=False, building_groups=True):
        """ function to import buildings from OBJ format

        Note:
            - the buildings must be divided in object like "o Building<building_id>"
            - OSM2World tool reference: http://osm2world.org/

        Args:
            obj_file_name_input: OBJ file path
            name_model_part: ModelPart name
            change_coord: True value to convert coord Y=-Z and Z=Y (useful for importing OBJ buildings obtained from OSM2World tool)
            building_groups: True value for a SubModelPart for each building; False value otherwise

        Returns:
            a ModelPart with the geometry of buildings
        """

        print("[DEBUG][geo_importer][ObjImportBuildings] START ObjImportBuildings\n")
        self._InitializeModelPart(name_model_part)

        if building_groups:
            bool_building = False   # boolean with value True if the geometry is a building
        else:
            bool_building = True    # in this case there are no groups. All Buildings are in a sub_model_part

        # read buildings model to extract verices and elements information
        with open (obj_file_name_input) as read_file:
            # vertex_coord = []		# here the coordinates of the node will be saved
            num_building = 1		# progressive number to count the number of buildings
            vertices_to_element = [0,0,0]

            # bool_building = False    # boolean with value True if the geometry is a building

            ID_vertex = 1
            ID_elem = 1

            current_sub_model_building = self.ModelPart # we set the ModelPart in current_sub_model_building if we haven't groups in obj file

            for row in read_file.readlines():
                row = row.split()
                if row:							# check if "row" is empty. To avoid an error about "out of range"
                    vertex_coord = []
                    if (row[0] == "v"):			# vertex: "v 0.0 0.0 0.0"
                        if ("#" in row):
                            # exclude the comment part
                            row = row[:4]
                        vertex_coord = [float(coord) for coord in row[1:]]
                        # for coord in row[1:]:
                        #     vertex_coord.append(float(coord))
                        
                        # we create new nodes in SubModelPart
                        if change_coord:
                            # when we importing OBJ file, there is a change in coordinate system: x = x; y = -z; z = y
                            node = self.ModelPart.CreateNewNode(ID_vertex, vertex_coord[0], -vertex_coord[2], vertex_coord[1])
                        else:
                            node = self.ModelPart.CreateNewNode(ID_vertex, vertex_coord[0], vertex_coord[1], vertex_coord[2])
                        
                        # we add the node in current sub model part
                        current_sub_model_building.AddNode(node, 0)
                        
                        vertex_coord = []
                        ID_vertex += 1
                    
                    elif (row[0] == "f") and bool_building:		# face: "f 1 2 3"
                        vertices_to_element[0] = int(row[1])
                        vertices_to_element[1] = int(row[2])
                        vertices_to_element[2] = int(row[3])
                        
                        # we add nodes into SubModelPart
                        # current_sub_model_building.AddNodes(vertices_to_element[0], vertices_to_element[1], vertices_to_element[2])
                        for node in vertices_to_element:
                            if (node not in current_sub_model_building.Nodes):
                                current_sub_model_building.AddNodes([node])
                    
                        # we create new elements in SubModelPart
                        current_sub_model_building.CreateNewElement("Element2D3N", ID_elem, [vertices_to_element[0], vertices_to_element[1], vertices_to_element[2]], self.ModelPart.GetProperties()[0])
                        ID_elem += 1
                    
                    elif (row[0] == "o"):			# when there is "o Building" we create a new sub_model_part because there is a new Building
                        if ("Building" in row[1]):
                            name_sub_model_building = "Building_{}".format(num_building)
                            current_sub_model_building = self.ModelPart.CreateSubModelPart(name_sub_model_building)
                            num_building += 1
                            bool_building = True
                        else:
                            current_sub_model_building = self.ModelPart
                            bool_building = False
        
        # we set all nodes as TO_ERASE
        for node in self.ModelPart.Nodes:
            node.Set(KratosMultiphysics.TO_ERASE,True)
        
        # for sub_model in self.ModelPart.SubModelParts:
        #     for node in sub_model.Nodes:
        #         node.Set(KratosMultiphysics.TO_ERASE, False)
        #         node.Set(KratosMultiphysics.VISITED, False)     # we set all nodes as not visited
        
        # we delete all free nodes
        for elem in self.ModelPart.Elements:
            for node in elem.GetNodes():
                node.Set(KratosMultiphysics.TO_ERASE, False)
                node.Set(KratosMultiphysics.VISITED, False)     # we set all nodes as not visited

        # we erase unused nodes
        self.ModelPart.RemoveNodesFromAllLevels(KratosMultiphysics.TO_ERASE)
        self.HasModelPart = True


    def ObjImport(self, obj_file_name_input, name_model_part="ModelPart"):
        """ import from OBJ file

        Args:
            obj_file_name_input: OBJ file path
            name_model_part: ModelPart name

        Returns:
            a ModelPart with the geometry
        """

        self._InitializeModelPart(name_model_part)

        with open (obj_file_name_input) as read_file:
            node_id = 1
            elem_id = 1

            for row in read_file.readlines():
                row = row.split()
                if row:     # check if "row" is empty
                    if (row[0] == "v"):
                        coord_x, coord_y, coord_z = [float(coord) for coord in row[1:]]
                        self.ModelPart.CreateNewNode(node_id, coord_x, coord_y, coord_z)
                        node_id += 1
                    
                    elif (row[0] == "f"):
                        node_1, node_2, node_3 = [int(nodes) for nodes in row[1:]]
                        self.ModelPart.CreateNewElement("Element2D3N", elem_id, [node_1, node_2, node_3], self.ModelPart.GetProperties()[1])
                        elem_id += 1
            
            self.HasModelPart = True


    ### --- function imports the nodes from a *.stl file --- ###

    def StlImport(self, stl_file_name_input, name_model_part="ModelPart"):

        self._InitializeModelPart(name_model_part)

        with open (stl_file_name_input) as read_file:

            node_dict = {}        # dictionary with all vertices. key = (x, y, z); value = node_id
            node_id = 1

            elem_dict = {}
            elem_id = 1

            for row in read_file.readlines():
                row = row.split()

                if (row[0] == "outer"):
                    elem_dict[elem_id] = []
                
                elif (row[0] == "endloop"):
                    elem_id += 1
                
                elif (row[0] == "vertex" and all(self._IsFloat(n) for n in row[1:])):
                    # row[0] = "vertex"; row[1] = x coordinate; row[2] = y coordinate; row[3] = z coordinate
                    X_coord, Y_coord, Z_coord = [float(coord) for coord in row[1:]]

                    if (X_coord, Y_coord, Z_coord) not in node_dict:
                        node_dict[X_coord, Y_coord, Z_coord] = node_id
                        elem_dict[elem_id].append(node_id)
                        node_id += 1

                    else:
                        # this vertex is already in dictionary
                        for coord, id in node_dict.items():
                            if (coord == (X_coord, Y_coord, Z_coord)):
                                elem_dict[elem_id].append(id)           # if the "coord" are already in the "node_dict", we append just his Id in "elem_dict" and we leave the loop
                                break

        # node_dict = self._DeleteNodeOnBase(node_dict)

        for coord, node_id in node_dict.items():
            self.ModelPart.CreateNewNode(node_id, coord[0], coord[1], coord[2])
        
        for id_elem, ids_node in elem_dict.items():
            self.ModelPart.CreateNewElement("Element2D3N", id_elem, [ids_node[0], ids_node[1], ids_node[2]], self.ModelPart.GetProperties()[1])
            # self.ModelPart.CreateNewElement("Element3D3N", id_elem, [ids_node[0], ids_node[1], ids_node[2]], self.ModelPart.GetProperties()[1])

        self.HasModelPart = True


    def AddStlImport (self, stl_file_name_input, current_building_model_part):
        """ This function is able to add an STL file into an existent ModelPart """

        with open (stl_file_name_input) as read_file:

            node_dict = {}        # dictionary with all vertices. key = (x, y, z); value = node_id
            node_id = current_building_model_part.NumberOfNodes() + 1

            elem_dict = {}
            elem_id = current_building_model_part.NumberOfElements() + 1

            for row in read_file.readlines():
                row = row.split()

                if (row[0] == "outer"):
                    elem_dict[elem_id] = []
                
                elif (row[0] == "endloop"):
                    elem_id += 1
                
                elif (row[0] == "vertex" and all(self._IsFloat(n) for n in row[1:])):
                    # row[0] = "vertex"; row[1] = x coordinate; row[2] = y coordinate; row[3] = z coordinate
                    X_coord, Y_coord, Z_coord = [float(coord) for coord in row[1:]]

                    if (X_coord, Y_coord, Z_coord) not in node_dict:
                        node_dict[X_coord, Y_coord, Z_coord] = node_id
                        elem_dict[elem_id].append(node_id)
                        node_id += 1

                    else:
                        # this vertex is already in dictionary
                        for coord, id in node_dict.items():
                            if (coord == (X_coord, Y_coord, Z_coord)):
                                elem_dict[elem_id].append(id)           # if the "coord" are already in the "node_dict", we append just his Id in "elem_dict" and we leave the loop
                                break

        # node_dict = self._DeleteNodeOnBase(node_dict)

        for coord, node_id in node_dict.items():
            current_building_model_part.CreateNewNode(node_id, coord[0], coord[1], coord[2])
        
        for id_elem, ids_node in elem_dict.items():
            current_building_model_part.CreateNewElement("Element2D3N", id_elem, [ids_node[0], ids_node[1], ids_node[2]], self.ModelPart.GetProperties()[1])

        self.HasModelPart = True


    ### --- function imports the nodes from a *.xyz file --- ###

    def XyzImport( self, xyz_file_name_input ):

        self._InitializeModelPart()

        with open (xyz_file_name_input) as read_file:
            node_id = 1
            for row in read_file.readlines():
                row = row.split()
                if ( all( self._IsFloat(n) for n in row) ):
                    # row[0] = x coordinate; row[1] = y coordinate; row[2] = z coordinate
                    X_coord, Y_coord, Z_coord = [float(coord) for coord in row[0:3]]
                    n = self.ModelPart.CreateNewNode(node_id, X_coord, Y_coord, Z_coord)
                    node_id += 1

        self.HasModelPart = True


########################################################################
# auxiliary functions
#########################################################################

    def _IsFloat( self, value):
        """ detecting if a conversion in float is possible """
        try:
            float( value )
            return True
        except ValueError:
            return False


    def _InitializeModelPart( self, name_model_part="ModelPart" ):
        """ this function adds all variables that will be needed during the pipeline (please add) """

        current_model = KratosMultiphysics.Model()
        self.ModelPart = current_model.CreateModelPart(name_model_part)
        self.ModelPart.SetBufferSize(1)

        self.ModelPart.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, 3)
        self.ModelPart.ProcessInfo.SetValue(KratosMultiphysics.TIME, 0.0)
        self.ModelPart.ProcessInfo.SetValue(KratosMultiphysics.DELTA_TIME, 1.0)

        self.ModelPart.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_H)
        self.ModelPart.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL)
        self.ModelPart.AddNodalSolutionStepVariable(KratosMultiphysics.IS_FLUID)
        self.ModelPart.AddNodalSolutionStepVariable(KratosMultiphysics.IS_STRUCTURE)
        self.ModelPart.AddNodalSolutionStepVariable(KratosMultiphysics.IS_FREE_SURFACE)
        self.ModelPart.AddNodalSolutionStepVariable(KratosMultiphysics.IS_BOUNDARY)
        self.ModelPart.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE)
        self.ModelPart.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE_GRADIENT)
        self.ModelPart.AddNodalSolutionStepVariable(KratosMultiphysics.FLAG_VARIABLE)

        self.ModelPart.Nodes.clear()
        self.ModelPart.Elements.clear()
        self.ModelPart.Conditions.clear()

        self.HasModelPart = True


    def _DeleteNodeOnBase(self, node_dict):
        """ this function delete all nodes with Z == 0 """

        for coords, _ in list(node_dict.items()):
            if (coords[2] == 0.0):    # if (Z == 0.0) the node is excluded
                del (node_dict[coords])

        return node_dict


    def _find_max_node_id(self):
        return max((node.Id for node in self.ModelPart.Nodes))


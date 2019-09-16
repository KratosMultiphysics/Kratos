"""
    # edit:     30 August 2019 -> add bool_building variable in ObjImport function
    # edit:     30 August 2019 -> in ObjImport, split [elif (row[0] == "o") and ("Building" in row[1]):]
    # edit:     03 September 2019 -> added change_coord variable to avoid/allow the change of the coordinates in ObjImport function
    # edit:     03 September 2019 -> create nodes and elements in self.ModelPart in ObjToPyMap function
"""

import KratosMultiphysics
from geo_processor import GeoProcessor

class GeoImporter( GeoProcessor ):

    ### --- overwriting the import function --- ###

    def SetGeoModelPart( self, modelPartIn ):

        KratosMultiphysics.Logger.PrintWarning("GeoImporter: SetGeoModelPart", "This function cannot be used. The model is created be reading a file.")


    def CreateBoxTetra(self, name_model_part, Xmin, Xmax, Ymin, Ymax, Zmin, Zmax):
        
        self._InitializeModelPart(name_model_part)
        # sub model part
        bottom = self.ModelPart.CreateSubModelPart("BottomModelPart")
        top = self.ModelPart.CreateSubModelPart("TopModelPart")
        lateral = self.ModelPart.CreateSubModelPart("LateralModelPart")
        # node
        node = self.ModelPart.CreateNewNode(1, Xmin, Ymin, Zmin)
        node.SetSolutionStepValue(KratosMultiphysics.DISTANCE, 0, 1e-7)		# set distance value
        node = self.ModelPart.CreateNewNode(2, Xmax, Ymin, Zmin)
        node.SetSolutionStepValue(KratosMultiphysics.DISTANCE, 0, 1e-7)		# set distance value
        node = self.ModelPart.CreateNewNode(3, Xmin, Ymax, Zmin)
        node.SetSolutionStepValue(KratosMultiphysics.DISTANCE, 0, 1e-7)		# set distance value
        node = self.ModelPart.CreateNewNode(4, Xmax, Ymax, Zmin)
        node.SetSolutionStepValue(KratosMultiphysics.DISTANCE, 0, 1e-7)		# set distance value
        node = self.ModelPart.CreateNewNode(5, Xmin, Ymin, Zmax)
        node.SetSolutionStepValue(KratosMultiphysics.DISTANCE, 0, 1.0)		# set distance value
        node = self.ModelPart.CreateNewNode(6, Xmax, Ymin, Zmax)
        node.SetSolutionStepValue(KratosMultiphysics.DISTANCE, 0, 1.0)		# set distance value
        node = self.ModelPart.CreateNewNode(7, Xmin, Ymax, Zmax)
        node.SetSolutionStepValue(KratosMultiphysics.DISTANCE, 0, 1.0)		# set distance value
        node = self.ModelPart.CreateNewNode(8, Xmax, Ymax, Zmax)
        node.SetSolutionStepValue(KratosMultiphysics.DISTANCE, 0, 1.0)		# set distance value
        # elements 3D4N
        prop = self.ModelPart.GetProperties()[1]
        self.ModelPart.CreateNewElement("Element3D4N", 1, [1, 3, 8, 4], prop)
        self.ModelPart.CreateNewElement("Element3D4N", 2, [3, 8, 7, 1], prop)
        self.ModelPart.CreateNewElement("Element3D4N", 3, [5, 7, 8, 1], prop)
        self.ModelPart.CreateNewElement("Element3D4N", 4, [1, 5, 6, 8], prop)
        self.ModelPart.CreateNewElement("Element3D4N", 5, [1, 2, 4, 8], prop)
        self.ModelPart.CreateNewElement("Element3D4N", 6, [2, 6, 8, 1], prop)
        # conditions 3D3N
        self.ModelPart.CreateNewCondition("SurfaceCondition3D3N",  1, [1, 2, 4], prop)		# bottom sub model part
        self.ModelPart.CreateNewCondition("SurfaceCondition3D3N",  2, [1, 4, 3], prop)		# bottom sub model part
        self.ModelPart.CreateNewCondition("SurfaceCondition3D3N",  3, [6, 5, 8], prop)		# top sub model part
        self.ModelPart.CreateNewCondition("SurfaceCondition3D3N",  4, [8, 5, 7], prop)		# top sub model part
        self.ModelPart.CreateNewCondition("SurfaceCondition3D3N",  5, [1, 5, 6], prop)		# lateral sub model part
        self.ModelPart.CreateNewCondition("SurfaceCondition3D3N",  6, [1, 6, 2], prop)		# lateral sub model part
        self.ModelPart.CreateNewCondition("SurfaceCondition3D3N",  7, [2, 6, 8], prop)		# lateral sub model part
        self.ModelPart.CreateNewCondition("SurfaceCondition3D3N",  8, [2, 8, 4], prop)		# lateral sub model part
        self.ModelPart.CreateNewCondition("SurfaceCondition3D3N",  9, [4, 8, 3], prop)		# lateral sub model part
        self.ModelPart.CreateNewCondition("SurfaceCondition3D3N", 10, [3, 8, 7], prop)		# lateral sub model part
        self.ModelPart.CreateNewCondition("SurfaceCondition3D3N", 11, [3, 7, 1], prop)		# lateral sub model part
        self.ModelPart.CreateNewCondition("SurfaceCondition3D3N", 12, [1, 7, 5], prop)		# lateral sub model part

        # fill bottom model part
        bottom.AddNodes([1, 2, 3, 4])
        bottom.AddConditions([1, 2])
        # fill top model part
        top.AddNodes([5, 6, 7, 8])
        top.AddConditions([3, 4])
        # fill lateral model part
        lateral.AddNodes([1, 2, 3, 4, 5, 6, 7, 8])
        lateral.AddConditions([5, 6, 7, 8, 9, 10, 11, 12])

        self.HasModelPart = True


    """ TEST!   FUNCTION UNDER CONSTRUCTION """
    def ObjToPyMap(self, obj_file_name_input):

        self._InitializeModelPart("building_model_part")

        with open (obj_file_name_input) as read_file:
            # vertex
            vertex_map = {}         # key: vertex_id; value: [x, y, z]
            vertex_id = 1
            
            # element
            elem_map = {}           # key: elem_id; value: [node1, node2, node3]
            elem_id = 1
            for row in read_file.readlines():
                row = row.split()
                if (row[0] == "v"):
                    x = float(row[1])
                    y = float(row[2])
                    z = float(row[3])
                    vertex_map[vertex_id] = [x, y, z]
                    self.ModelPart.CreateNewNode(vertex_id, x, y, z)  # added on 02/09/2019
                    vertex_id += 1
                
                elif (row[0] == "f"):
                    n1 = int(row[1])
                    n2 = int(row[2])
                    n3 = int(row[3])
                    elem_map[elem_id] = [n1, n2, n3]
                    self.ModelPart.CreateNewElement("Element2D3N", elem_id, [n1, n2, n3], self.ModelPart.GetProperties()[1])   # added on 02/09/2019
                    elem_id += 1
        
        self.HasModelPart = True        # added on 02/09/2019
        
        # probably the best solution is to fill a ModelPart. I will try also this way!
        return (vertex_map, elem_map)

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
        
        self.HasModelPart = True

        #########################################################################################################
        geom_id = 1         # a counter
        dict_visited = []	# list with the elements in the dictionary already visited

        num_elem = self.ModelPart.NumberOfElements()    # key_list = list(elem_map.keys())	# it is necessary because we delete elements in the dictionary; so we can't iterate on dictionary directly
        """for faces in key_list:"""
        for faces in self.ModelPart.Elements:
            if (faces.Id in dict_visited):
                continue
            dict_visited.append(faces.Id)

            next_faces = [faces.Id]	# the netx faces that will be processed

            node_visited = []	# list with nodes already visited in this geometry
            elem_visited = []	# list with elements already visited in this geometry
            while (next_faces):	# this loop ends when next_faces is empty
                current_face = next_faces[0]			# get the first element in next_faces up to the last element. In next operations we subtract elements in next_faces
                """for node in elem_map[current_face]:"""
                for node in self.ModelPart.GetElement(current_face):        # CHECK THIS IF IT IS CORRECT
                    elem_visited.append(current_face)	# we update elem_visited
                    dict_visited.append(current_face)
                    if (node.Id in node_visited):
                        continue						# we go on if the node is already visited
                    node_visited.append(node.Id)			# we update node_visited
                    """for id, elem in elem_map.items():"""
                    for elem in self.ModelPart.Elements:
                        if (node in elem):
                            next_faces.append(id)
                next_faces = list(set(next_faces) - set(elem_visited))	# we remove duplicates and the elements that are already visited

                del(elem_map[current_face])			# we delete current_face in elem_map

            geometries[geom_id] = list(set(elem_visited))			# key: id geometry; value: id elements of this geometry
            geom_id += 1
        #########################################################################################################



    """ FUNCTION UNDER CONSTRUCTION """
    def ObjImport(self, obj_file_name_input, name_model_part="ModelPart", change_coord=True):
        # a sub_model_part is created for each Building

        self._InitializeModelPart(name_model_part)
        print("\n***** START ObjImport: geo_importer/ObjImport *****\n")

        # read buildings model to extract verices and elements information
        with open (obj_file_name_input) as read_file:
            vertex_coord = []		# here the coordinates of the node will be saved
            num_building = 1		# progressive number to count the number of buildings
            vertices_to_element = [0,0,0]

            bool_building = False    # boolean with value True if the geometry is a building

            ID_vertex = 1
            ID_elem = 1

            current_sub_model_building = self.ModelPart # we set the ModelPart in current_sub_model_building if we haven't groups in obj file

            # print("***** NUMBER OF LINES: ", len(read_file.readlines()))

            for row in read_file.readlines():
                row = row.split()
                if row:							# check if "row" is empty. To avoid an error about "out of range"
                    if (row[0] == "v"):			# vertex: "v 0.0 0.0 0.0"
                        for coord in row[1:]:
                            vertex_coord.append(float(coord))
                        
                        # we create new nodes in SubModelPart
                        if change_coord:
                            # when we importing OBJ file, there is a change in coordinate system: x = x; y = -z; z = y
                            self.ModelPart.CreateNewNode(ID_vertex, vertex_coord[0], -vertex_coord[2], vertex_coord[1])
                        else:
                            self.ModelPart.CreateNewNode(ID_vertex, vertex_coord[0], vertex_coord[1], vertex_coord[2])
                        
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
                        current_sub_model_building.CreateNewElement("Element2D3N", ID_elem, [vertices_to_element[0], vertices_to_element[1], vertices_to_element[2]], self.ModelPart.GetProperties()[2])
                        ID_elem += 1
                    
                    elif (row[0] == "o"):			# when there is 'o Building' we create a new sub_model_part because there is a new Building
                        if ("Building" in row[1]):
                            name_sub_model_building = "Building_{}".format(num_building)
                            current_sub_model_building = self.ModelPart.CreateSubModelPart(name_sub_model_building)
                            num_building += 1
                            bool_building = True
                        else:
                            bool_building = False

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


    """ This function is able to add an STL file into an existent ModelPart """
    def AddStlImport (self, stl_file_name_input, current_building_model_part):

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
        # detecting if a conversion in float is possible
        try:
            float( value )
            return True
        except ValueError:
            return False


    def _InitializeModelPart( self, name_model_part="ModelPart" ):
        # this function adds all variables that will be needed during the pipeline (please add)

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


    def _DeleteNodeOnBase(self, node_dict):
        # this function delete all nodes with Z == 0

        for coords, _ in list(node_dict.items()):
            if (coords[2] == 0.0):    # if (Z == 0.0) the node is excluded
                del (node_dict[coords])

        return node_dict

import KratosMultiphysics
from geo_processor import GeoProcessor

class GeoImporter( GeoProcessor ):

    ### --- overwriting the import function --- ###

    def SetGeoModelPart( self, modelPartIn ):

        KratosMultiphysics.Logger.PrintWarning("GeoImporter: SetGeoModelPart", "This function cannot be used. The model is created be reading a file.")


    ''' FUNCTION UNDER CONSTRUCTION '''
    def ObjImport(self, obj_file_name_input, name_model_part="ModelPart"):
        # a sub_model_part is created for each Building

        self._InitializeModelPart(name_model_part)

        # read buildings model to extract verices and elements information
        with open (obj_file_name_input) as read_file:
            vertex_coord = []		# here the coordinates of the node will be saved
            num_building = 1		# progressive number to count the number of buildings
            vertices_to_element = [0,0,0]

            ID_vertex = 1
            ID_elem = 1

            current_sub_model_building = self.ModelPart # we set the ModelPart if we haven't groups in obj file

            for row in read_file.readlines():
                row = row.split()
                if row:							# check if "row" is empty. To avoid an error about "out of range"
                    if (row[0] == "v"):			# vertex: "v 0.0 0.0 0.0"
                        for coord in row[1:]:
                            vertex_coord.append(float(coord))
                        
                        # we create new nodes in SubModelPart
                        # when we importing OBJ file, there is a change in coordinate system: x = x; y = -z; z = y
                        self.ModelPart.CreateNewNode(ID_vertex, vertex_coord[0], -vertex_coord[2], vertex_coord[1])
                        vertex_coord = []
                        ID_vertex += 1
                    
                    elif (row[0] == "f"):		# face: "f 1 2 3"
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
                    
                    elif (row[0] == "o") and ("Building" in row[1]):			# when there is 'o Building' we create a new sub_model_part because there is a new Building
                        name_sub_model_building = "Building_{}".format(num_building)
                        current_sub_model_building = self.ModelPart.CreateSubModelPart(name_sub_model_building)
                        num_building += 1

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

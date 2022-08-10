# CoSimulation imports
from re import I
import KratosMultiphysics as KM
import KratosMultiphysics.CoSimulationApplication as KratosCoSim
import KratosMultiphysics.StructuralMechanicsApplication # needed for some variables

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.solver_wrappers.external.remote_controlled_mdpa_solver_wrapper import RemoteControlledWithModalPartSolverWrapper

# Other imports
from KratosMultiphysics.CoSimulationApplication.utilities import model_part_utilities
import os

def Create(settings, model, solver_name):
    return PEDFlowSolverWrapper(settings, model, solver_name)

class PEDFlowSolverWrapper(RemoteControlledWithModalPartSolverWrapper):
    """This class is a generic wrapper for connecting external solvers that are being remote controlled
    """
    def __init__(self, settings, model, solver_name):
        super().__init__(settings, model, solver_name)
        self.boundary_part_names = settings["boundary_part_names"].GetStringArray()
        self.boundary_part_types = settings["boundary_part_types"].GetStringArray()
        self.coupling_part_name = settings["coupling_part_name"].GetString()
        self.coupling_part_file = settings["coupling_part_file"].GetString()

        self.auxiliary_model = KM.Model()
        self.auxiliary_model_part = self.auxiliary_model.CreateModelPart(self.model_part_name)
        KM.ModelPartIO(self.coupling_part_file).ReadModelPart(self.auxiliary_model_part)

    def Initialize(self):
        super().Initialize()
        # self.__TranslateMDPAToFECAD()
        # print("New pedflow.inp.domn file created.")
        

    @classmethod
    def _GetDefaultParameters(cls):
        return KM.Parameters("""{
            "type"                    : "",
            "solver_wrapper_settings" : {},
            "io_settings"             : {},
            "data"                    : {},
            "echo_level"              : 0,
            "model_part_file"         : "",
            "model_part_name"         : "",
            "mpi_settings"            : {},
            "boundary_part_names"     : [],
            "boundary_part_types"     : [],
            "coupling_part_name"      : "",
            "coupling_part_file"      : "",
            "interface_submodel_part" : ""
        }""")

    def __TranslateMDPAToFECAD(self):
        '''
            Creates the FECAD file from a ModelPart
        '''
        # File management
        filePath = os.getcwd() + '/pedflow.inp.domn'
        coupling_model_part = self.auxiliary_model_part.GetSubModelPart(self.coupling_part_name)

        with open(filePath, "w") as f:

            # Header 
            f.write("     nelem     npoin     nboun\n")
            nelem = coupling_model_part.NumberOfElements()
            npoin = coupling_model_part.NumberOfNodes()
            nboun = 0
            for submodel_name in self.boundary_part_names:
                nboun += self.auxiliary_model_part.GetSubModelPart(submodel_name).NumberOfNodes()
            nboun -= 8 # For rectangular plates with 4 corners, find better solution

            f.write("     {}     {}     {}  {:.3f}\n".format(nelem, npoin, nboun, 0.000))

            # Element information
            # ID(ielem) 3-Nodes(ip1-ip3) inreme material
            f.write(" ielem ip1-ip3 inreme material\n")

            ## Node and element numbering maps
            element_map = dict()
            node_map = dict()
            i = 1
            for node in coupling_model_part.GetNodes():
                node_map[node.Id] = i
                i += 1
            inv_node_map = {v: k for k, v in node_map.items()}

            j = 1
            for element in coupling_model_part.GetElements():
                element_map[element.Id] = j
                entry = str(element_map[element.Id])
                ## Change node numbering to clockwise
                nodes_element = [None] * 3
                k = 2
                for node in element.GetNodes():
                    nodes_element[k] = node_map[node.Id]
                    k -= 1
                for node in nodes_element:
                    entry = entry + " {}".format(node)
                entry = entry + " 0 0 0 0 0 0 0 1\n"
                f.write(entry)
                j += 1
            
            # Node information
            # ID(ipoin) x y z
            f.write(" ipoin x y z\n")

            for node in coupling_model_part.GetNodes():
                entry = " {} {} {} {}\n".format(node_map[node.Id], node.X, node.Y, node.Z)
            
                f.write(entry)
            
            # ID(ipoin) unkno
            f.write(" ipoin unkno\n")

            for node in coupling_model_part.GetNodes():
                entry = " {} 0.00 0.00 0.00 0.00 0.00 0.00\n".format(node_map[node.Id])
                f.write(entry)

            # Boundary information
            # Heavily limited by the geometry, it needs attention and manual input to keep the signs coherent.
            f.write(" boundary points\n")

            dict_of_corners = dict()

            # Save and write the corners

            for submodel_name, submodel_type in zip(self.boundary_part_names, self.boundary_part_types):
                if submodel_type == "corner_input":
                    submodel_part = self.auxiliary_model_part.GetSubModelPart(submodel_name)
                    f.write(" {} 4 1 0 {} 0 0.000\n".format(node_map[[node for node in submodel_part.GetNodes()][0].Id], submodel_name.split("_")[-1]))
                    dict_of_corners["Line_" + submodel_name.split("_")[-1]] = node_map[[node for node in submodel_part.GetNodes()][0].Id]
                elif submodel_type == "corner_output":
                    submodel_part = self.auxiliary_model_part.GetSubModelPart(submodel_name)
                    f.write(" {} 4 2 0 {} 0 0.000\n".format(node_map[[node for node in submodel_part.GetNodes()][0].Id], submodel_name.split("_")[-1]))
                    dict_of_corners["Line_" + submodel_name.split("_")[-1]] = node_map[[node for node in submodel_part.GetNodes()][0].Id]

            # Save and write the lines

            for submodel_name, submodel_type in zip(self.boundary_part_names, self.boundary_part_types):
                if not (submodel_type == "corner_input" or submodel_type == "corner_output"):
                    submodel_part = self.auxiliary_model_part.GetSubModelPart(submodel_name)
                    line_number =  submodel_name.split("_")[-1]
                    end_corner = dict_of_corners[submodel_name]

                    ## Calculate line parameters for xi values

                    dict_of_nodes  = dict()
                    for node in submodel_part.GetNodes():
                        dict_of_nodes[node_map[node.Id]] = [node.X, node.Y, node.Z]
                    max_X_id = max(dict_of_nodes, key= lambda x: dict_of_nodes[x][0])
                    min_X_id = min(dict_of_nodes, key= lambda x: dict_of_nodes[x][0])
                    len_X = abs(dict_of_nodes[max_X_id][0]-dict_of_nodes[min_X_id][0])
                    if len_X > 1E-5:
                        xi_calculation_mode = 0
                        if min_X_id == end_corner:
                            xi_direction = 1 # positive = 1, negative = -1
                        elif max_X_id == end_corner:
                            xi_direction = -1
                        else:
                            raise Exception("End point of " + submodel_name + " do not coincide with indicated corner.") 

                    else:
                        max_Y_id = max(dict_of_nodes, key= lambda x: dict_of_nodes[x][1])
                        min_Y_id = min(dict_of_nodes, key= lambda x: dict_of_nodes[x][1])
                        len_Y = abs(dict_of_nodes[max_Y_id][1]-dict_of_nodes[min_Y_id][1])
                        if len_Y > 1E-5:
                            xi_calculation_mode = 1
                            if min_Y_id == end_corner:
                                xi_direction = 1 # positive = 1, negative = -1
                            elif max_Y_id == end_corner:
                                xi_direction = -1
                            else:
                                raise Exception("End point of " + submodel_name + " do not coincide with indicated corner.") 
                        else:
                            max_Z_id= max(dict_of_nodes, key= lambda x: dict_of_nodes[x][2])
                            min_Z_id = min(dict_of_nodes, key= lambda x: dict_of_nodes[x][2])
                            len_Z = abs(dict_of_nodes[max_Z_id][2]-dict_of_nodes[min_Z_id][2])

                            if len_Z < 1E-5:
                                raise Exception("Error generating FECAD file: " + submodel_name + " is empty.")
                            xi_calculation_mode = 2
                            if min_Z_id == end_corner:
                                xi_direction = 1 # positive = 1, negative = -1
                            elif max_Z_id == end_corner:
                                xi_direction = -1
                            else:
                                raise Exception("End points of " + submodel_name + " do not coincide with indicated corner.") 

                    ## Calculate xi values

                    list_of_nodes = [[node_map[node.Id], node.X, node.Y, node.Z] for node in submodel_part.GetNodes() if not node_map[node.Id] in dict_of_corners.values()]
                    list_of_xi = [None] * len(list_of_nodes)
                    if xi_calculation_mode == 0:
                        end_corner_X = submodel_part.GetNode(inv_node_map[end_corner]).X
                        for i in range(len(list_of_xi)):
                            list_of_xi[i] = (list_of_nodes[i][1]-end_corner_X)/len_X*xi_direction

                    elif xi_calculation_mode == 1:
                        end_corner_Y = submodel_part.GetNode(inv_node_map[end_corner]).Y
                        for i in range(len(list_of_xi)):
                            list_of_xi[i] = (list_of_nodes[i][2]-end_corner_Y)/len_Y*xi_direction

                    elif xi_calculation_mode == 2:
                        end_corner_Z = submodel_part.GetNode(inv_node_map[end_corner]).Z
                        for i in range(len(list_of_xi)):
                            list_of_xi[i] = (list_of_nodes[i][3]-end_corner_Z)/len_Z*xi_direction
                    
                    # Print entries

                    if submodel_type == "input":
                        for node, xi in zip(list_of_nodes,list_of_xi):
                            f.write(" {} 4 1 1 {} 0 {}\n".format(node[0], line_number, xi))

                    elif submodel_type == "output":
                        for node, xi in zip(list_of_nodes,list_of_xi):
                            f.write(" {} 4 2 1 {} 0 {}\n".format(node[0], line_number, xi)) 

                    elif submodel_type == "wall":
                        for node, xi in zip(list_of_nodes,list_of_xi):
                            f.write(" {} 0 0 1 {} 0 {}\n".format(node[0], line_number, xi))


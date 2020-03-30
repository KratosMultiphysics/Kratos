import KratosMultiphysics as KratosMultiphysics
import KratosMultiphysics.CableNetApplication as CableNetApplication

from KratosMultiphysics import Logger

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return EdgeCableElementProcess(Model, settings["Parameters"])

class custom_node:
    def __init__(self,start_distance,kratos_node):
        self.start_distance = start_distance
        self.kratos_node = kratos_node
    def return_distance_to_line_start(self):
        return self.start_distance

def return_node_distance_to_line_start(node):
    return node.return_distance_to_line_start()



class EdgeCableElementProcess(KratosMultiphysics.Process):

    def __init__(self, Model, settings ):
        KratosMultiphysics.Process.__init__(self)
        default_settings = KratosMultiphysics.Parameters("""
        {
            "edge_sub_model_part_name"  : "Structure.example_part",
            "element_type"              : "cable",
            "node_id_order"             : [1,2,3],
            "element_id"                : 1,
            "property_id"               : 1
        }
        """)
        default_settings.ValidateAndAssignDefaults(settings)

        self.edge_model_part = Model[settings["edge_sub_model_part_name"].GetString()]

        node_list = settings["node_id_order"].GetVector()
        if len(node_list)==0:
            node_list = self.CreateCorrectNodeOrder()
            settings["node_id_order"].SetVector(node_list)

        self.edge_cable_element_process = CableNetApplication.EdgeCableElementProcess(self.edge_model_part, settings)


    def ExecuteInitialize(self):
        self.edge_cable_element_process.ExecuteInitialize()
        Logger.PrintInfo("Initialized","EdgeCableElementProcess")

    def CreateCorrectNodeOrder(self):
        ## find start/end nodes and calculate total distance
        max_distance,end_points = 0, []
        for node_i in self.edge_model_part.Nodes:
            for node_j in self.edge_model_part.Nodes:
                distance_i = (node_i.X0 - node_j.X0)*(node_i.X0 - node_j.X0)
                distance_i += (node_i.Y0 - node_j.Y0)*(node_i.Y0 - node_j.Y0)
                distance_i += (node_i.Z0 - node_j.Z0)*(node_i.Z0 - node_j.Z0)
                distance_i = distance_i**0.5

                if distance_i>max_distance:
                    max_distance=distance_i
                    end_points = [node_i,node_j]

        ## create sorted node_list
        custom_node_list = []
        for node_i in self.edge_model_part.Nodes:
            distance_i = (node_i.X0 - end_points[0].X0)*(node_i.X0 - end_points[0].X0)
            distance_i += (node_i.Y0 - end_points[0].Y0)*(node_i.Y0 - end_points[0].Y0)
            distance_i += (node_i.Z0 - end_points[0].Z0)*(node_i.Z0 - end_points[0].Z0)
            distance_i = distance_i**0.5
            custom_node_i = custom_node(distance_i,node_i)
            custom_node_list.append(custom_node_i)

        sorted_node_list = sorted(custom_node_list, key=return_node_distance_to_line_start)

        return [node.kratos_node.Id for node in sorted_node_list]
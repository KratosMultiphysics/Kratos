import KratosMultiphysics as Kratos

class PouliotRecoveryTools:
    def __init__(self):
        pass

    def MakeRecoveryModelPart(self, model_part):
        edges_model_part = Kratos.ModelPart("Edges")
        set_of_all_edges = set()
        for elem in model_part.Elements:
            for i, first_node in enumerate(elem.Nodes[:-1]):
                for j, second_node in enumerate(elem.Nodes[i:]):
                    edge_ids = (first_node.Id, second_node.Id)
                    set_of_all_edges.add(edge_ids)
        for i, edge in enumerate(set_of_all_edges):
            edges_model_part.CreateNewElement("Element3D1N", i, edge, edges_model_part.GetProperties()[0])

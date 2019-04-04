class EmpireMapper():
    def __init__(self, origin_model_part, destination_model_part, settings):
        self.origin_model_part = origin_model_part
        self.destination_model_part = destination_model_part

    def Initialize(self):
        pass

    def Update(self):
        pass

    def Map(self, origin_variable, destination_variable):
        origin_nodes = self.origin_model_part.Nodes
        destination_nodes = self.destination_model_part.Nodes
        for node_origin, node_destination in zip(origin_nodes, destination_nodes):
            value = node_origin.GetSolutionStepValue(origin_variable)
            node_destination.SetSolutionStepValue(destination_variable, value)

    def InverseMap(self, destination_variable, origin_variable):
        origin_nodes = self.origin_model_part.Nodes
        destination_nodes = self.destination_model_part.Nodes
        for node_origin, node_destination in zip(origin_nodes, destination_nodes):
            value = node_destination.GetSolutionStepValue(destination_variable)
            node_origin.SetSolutionStepValue(origin_variable, value)

import KratosMultiphysics

class CopyValueModeler(object):

    def __init__(self, model, settings):
        self.model = model
        self.settings = settings

    def SetupGeometryModel(self):
        pass

    def PrepareGeometryModel(self):
        pass

    def SetupModelPart(self):
        orig_model_part = self.model.GetModelPart(self.settings["origin_model_part"].GetString())
        dest_model_part = self.model.GetModelPart(self.settings["destination_model_part"].GetString())

        n_nodes_orig = orig_model_part.NumberOfNodes()
        n_nodes_dest = dest_model_part.NumberOfNodes()
        if n_nodes_orig != n_nodes_dest:
            err_msg = "Origin model part number of nodes is {} while destination model part one is {}. Number of nodes must match.".format(n_nodes_orig, n_nodes_dest)
            raise Exception(err_msg)

        var_list = []
        for var_name in self.settings["variable_names"].GetStringArray():
            var_list.append(KratosMultiphysics.KratosGlobals.GetVariable(var_name))

        for orig_node, dest_node in zip(orig_model_part.Nodes, dest_model_part.Nodes):
            for var in var_list:
                dest_node.SetSolutionStepValue(var, orig_node.GetSolutionStepValue(var))

def Factory(model, settings):
    return CopyValueModeler(model, settings)

import numpy as np
import copy

import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools
cs_data_structure = cs_tools.cs_data_structure


def Create(parameters):
    return MapperPermutation(parameters)


# Class MapperPermutation: Permutation of coordinates.
class MapperPermutation(object):
    def __init__(self, parameters):
        """
        TODO: description
        """
        super().__init__()

        self.settings = parameters['settings']
        self.interpolator = False

        self.permutation = []
        for par in self.settings['permutation']:
            self.permutation.append(par.GetInt())

    def Initialize(self, model_part_from, model_part_to):
        # *** idea: make it more elegant by asking model_part_in and bool 'forward' as input

        permutation = self.permutation
        if model_part_to is None and model_part_from is not None:
            model_part_in = model_part_from
        elif model_part_to is not None and model_part_from is None:
            model_part_in = model_part_to
            permutation = np.argsort(self.permutation)
        else:
            raise ValueError('exactly one ModelPart must be None')

        model = cs_data_structure.Model()
        model_part_out = model.CreateModelPart('no_name')
        model_part_out._ModelPart__hist_variables = model_part_in._ModelPart__hist_variables
        for node in model_part_in.Nodes:
            coords = np.array([node.X, node.Y, node.Z])
            model_part_out.CreateNewNode(node.Id, *tuple(coords[permutation]))

        return model_part_out

    def Finalize(self):
        pass

    def __call__(self, args_from, args_to):
        model_part_from, var_from = args_from
        model_part_to, var_to = args_to

        for node_from, node_to in zip(model_part_from.Nodes, model_part_to.Nodes):
            value = node_from.GetSolutionStepValue(var_from)
            node_to.node.SetSolutionStepValue(var_to, 0, value)

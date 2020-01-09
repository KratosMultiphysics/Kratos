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
        This is not an interpolator, but a transformer.
        This is denoted by setting the self.interpolator
        attribute to False in the __init__.

        The difference is that a transformer is
        initialized from only one side (from or to)
        and returns the ModelPart corresponding
        to the forward or backward transformation.
        An interpolator is initialized from both
        sides and returns nothing.

        It can be initialized from both sides,
        based on the forward parameter.
        If forward == True, then the model_part_from
        is expected as input.
        Else, the model_part_to is expected.

        The historical variables are not changed,
        simply copied from input to output ModelPart.
        The coordinates are permutated according to the
        permutation parameter (list of ints) in the
        JSON file.
        """
        super().__init__()

        self.settings = parameters['settings']
        self.interpolator = False

        self.permutation = []
        for par in self.settings['permutation']:
            self.permutation.append(par.GetInt())

    def Initialize(self, model_part_in, forward):
        permutation = self.permutation
        if not forward:
            permutation = np.argsort(permutation)

        model = cs_data_structure.Model()
        model_part_out = model.CreateModelPart('no_name')
        model_part_out._ModelPart__hist_variables = model_part_in._ModelPart__hist_variables
        for node in model_part_in.Nodes:
            coords = np.array([node.X0, node.Y0, node.Z0])
            model_part_out.CreateNewNode(node.Id, *tuple(coords[permutation]))

        return model_part_out

    def Finalize(self):
        pass

    def __call__(self, args_from, args_to):
        model_part_from, var_from = args_from
        model_part_to, var_to = args_to

        for node_from, node_to in zip(model_part_from.Nodes, model_part_to.Nodes):
            value = node_from.GetSolutionStepValue(var_from)

            if var_from.Type() == 'Double':
                pass
            elif var_from.Type() == 'Array':
                value = list(np.array(value)[self.permutation])
            else:
                raise NotImplementedError(f'Mapping not yet implemented for Variable of Type {var_from.Type()}.')

            node_to.SetSolutionStepValue(var_to, 0, value)

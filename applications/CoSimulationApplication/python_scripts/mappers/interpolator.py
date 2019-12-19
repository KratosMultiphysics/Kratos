from scipy.spatial import cKDTree
import numpy as np

import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools
cs_data_structure = cs_tools.cs_data_structure


def Create(parameters):
    raise NotImplementedError('this class can only be used as super-class')


# Class MapperInterpolator: base class for interpolators.
class MapperInterpolator():
    def __init__(self, parameters):
        # store settings
        self.settings = parameters['settings']
        self.interpolator = True
        self.balanced_tree = self.settings['balanced_tree'].GetBool()
        self.n_nearest = 0  # must be set in sub-class!

        # get list with directions
        self.directions = []
        for direction in self.settings['directions'].list():
            tmp = direction.GetString().upper()
            if tmp not in ['X', 'Y', 'Z']:
                raise ValueError(f'"{tmp}" is not a valid direction.')
            self.directions.append(tmp + '0')

    def Initialize(self, model_part_from, model_part_to):
        # get coords_from
        self.n_from = model_part_from.NumberOfNodes()
        self.coords_from = np.zeros((self.n_from, len(self.directions)))
        for i, node in enumerate(model_part_from.Nodes):
            for j, direction in enumerate(self.directions):
                self.coords_from[i, j] = getattr(node, direction)

        # get coords_to
        self.n_to = model_part_to.NumberOfNodes()
        self.coords_to = np.zeros((self.n_to, len(self.directions)))
        for i, node in enumerate(model_part_to.Nodes):
            for j, direction in enumerate(self.directions):
                self.coords_to[i, j] = getattr(node, direction)

        # check if n_from is large enough
        if self.n_from < self.n_nearest:
            raise ValueError(f'not enough from-points: {self.n_from} < {self.n_nearest}')

        # build and query tree
        if self.balanced_tree:  # time-intensive
            tree = cKDTree(self.coords_from)
        else:  # less stable
            tree = cKDTree(self.coords_from, balanced_tree=False)
        self.distances, self.nearest = tree.query(self.coords_to, k=self.n_nearest)

    def Finalize(self):
        pass

    def __call__(self, args_from, args_to):
        model_part_from, var_from = args_from
        model_part_to, var_to = args_to

        # check if both Variables have same Type
        if var_from.Type() != var_to.Type():
            raise TypeError('Variables to be mapped have different Type.')

        # scalar interpolation
        if var_from.Type() == 'Double':
            hist_var_from = np.zeros(self.n_from)
            for i, node in enumerate(model_part_from.Nodes):
                hist_var_from[i] = node.GetSolutionStepValue(var_from)

            for i, node in enumerate(model_part_to.Nodes):
                hist_var_to = np.dot(self.coeffs[i], hist_var_from[self.nearest[i, :]])
                node.SetSolutionStepValue(var_to, 0, hist_var_to)

        # vector interpolation
        elif var_from.Type() == 'Array':
            hist_var_from = np.zeros((self.n_from, 3))
            for i, node in enumerate(model_part_from.Nodes):
                hist_var_from[i] = node.GetSolutionStepValue(var_from)

            for i, node in enumerate(model_part_to.Nodes):
                hist_var_to = [0., 0., 0.]
                for j in range(3):
                    hist_var_to[j] = np.dot(self.coeffs[i],
                                            hist_var_from[self.nearest[i, :], j])
                node.SetSolutionStepValue(var_to, 0, hist_var_to)

        # other types of Variables
        else:
            raise NotImplementedError(f'Mapping not yet implemented for Variable of Type {var_from.Type()}.')

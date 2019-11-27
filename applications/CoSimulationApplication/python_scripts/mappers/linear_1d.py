from scipy.spatial import cKDTree
import numpy as np

import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools
cs_data_structure = cs_tools.cs_data_structure


def Create(parameters):
    return MapperLinear1D(parameters)


# Class MapperLinear: Linear interpolation in 1D.
class MapperLinear1D(object):
    def __init__(self, parameters):
        """
        This mapper uses the two closest points; this does
        not mean that the point to be interpolated lies in
        between them, e.g. when there is grid stretching.
        This can lead to discontinuities if grid spacing
        changes quickly.

        At the boundaries of the domain, linear
        extrapolation is used.
        """
        super().__init__()

        self.settings = parameters['settings']

        self.balanced_tree = self.settings['balanced_tree'].GetBool()
        self.coord = self.settings['direction'].GetString().upper()
        if self.coord not in ['X', 'Y', 'Z']:
            raise ValueError(f'{self.coord} is not a valid direction.')

    def Initialize(self, model_part_from, model_part_to):
        self.n_from = model_part_from.NumberOfNodes()
        coords_from = np.zeros((self.n_from, 1))
        for i, node in enumerate(model_part_from.Nodes):
            coords_from[i] = getattr(node, self.coord)

        self.n_to = model_part_to.NumberOfNodes()
        coords_to = np.zeros((self.n_to, 1))
        for i, node in enumerate(model_part_to.Nodes):
            coords_to[i] = getattr(node, self.coord)

        # build and query tree
        if self.balanced_tree:  # time-intensive
            tree = cKDTree(coords_from)
        else:  # less stable
            tree = cKDTree(coords_from, balanced_tree=False)
        _, self.nearest = tree.query(coords_to, k=2, n_jobs=-1)

        # linear interpolation/extrapolation from nearest points
        self.coeffs = np.zeros((self.n_to, 2))
        for i in range(self.n_to):
            x = coords_to[i]
            x_1 = coords_from[self.nearest[i, 0]]
            x_2 = coords_from[self.nearest[i, 1]]
            if x_1 == x_2:
                print(x, x_1, x_2)

            c_1 = (x - x_2) / (x_1 - x_2)
            c_2 = (x - x_1) / (x_2 - x_1)
            self.coeffs[i, :] = [c_1, c_2]

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

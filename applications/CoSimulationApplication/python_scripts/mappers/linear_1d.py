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

    def Initialize(self, model_part_from, model_part_to):
        coord = self.settings['direction'].GetString().upper()
        if coord not in ['X', 'Y', 'Z']:
            raise ValueError(f'{coord} is not a valid direction.')

        coords_from = np.zeros((model_part_from.NumberOfNodes(), 1))
        for i, node in enumerate(model_part_from.Nodes):
            coords_from[i] = getattr(node, coord)

        n_to = model_part_to.NumberOfNodes()
        coords_to = np.zeros((n_to, 1))
        for i, node in enumerate(model_part_to.Nodes):
            coords_to[i] = getattr(node, coord)

        tree = cKDTree(coords_from)  # time-intensive part
        distance, self.nearest = tree.query(coords_to, k=2, n_jobs=-1)

        # linear interpolation/extrapolation from nearest points
        x = coords_to
        x0 = coords_from[self.nearest[:, 0]]
        x1 = coords_from[self.nearest[:, 1]]
        c0 = (x - x1) / (x0 - x1)
        c1 = (x - x0) / (x1 - x0)
        self.coeffs = np.hstack((c0.reshape(-1, 1), c1.reshape(-1, 1)))

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
            hist_var_from = np.zeros(model_part_from.NumberOfNodes())
            for i, node in enumerate(model_part_from.Nodes):
                hist_var_from[i] = node.GetSolutionStepValue(var_from)

            for i, node in enumerate(model_part_to.Nodes):
                hist_var_to = np.dot(self.coeffs[i], hist_var_from[self.nearest[i, :]])
                node.SetSolutionStepValue(var_to, 0, hist_var_to)

        # vector interpolation
        elif var_from.Type() == 'Array':
            hist_var_from = np.zeros((model_part_from.NumberOfNodes(), 3))
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

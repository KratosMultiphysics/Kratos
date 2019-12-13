from scipy.spatial import cKDTree
import numpy as np

import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools
from KratosMultiphysics.CoSimulationApplication.mappers.nearest import MapperNearest
cs_data_structure = cs_tools.cs_data_structure


def Create(parameters):
    return MapperLinear1D(parameters)


# Class MapperLinear: Linear interpolation in 1D.
class MapperLinear1D(MapperNearest):
    def __init__(self, parameters):
        """
        This mapper uses the two closest points; this does
        not mean that the point to be interpolated lies in
        between them, e.g. when there is grid stretching.
        This can lead to discontinuities if grid spacing
        changes quickly.

        At the boundaries of the domain, linear
        extrapolation is used.

        The __call__ and Finalize methods are
        inherited from MapperNearest.
        """

        self.settings = parameters['settings']
        self.interpolator = True
        self.balanced_tree = self.settings['balanced_tree'].GetBool()
        self.coord = self.settings['direction'].GetString().upper() + '0'
        if self.coord not in ['X0', 'Y0', 'Z0']:
            raise ValueError(f'{self.coord[:-1]} is not a valid direction.')

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
        _, self.nearest = tree.query(coords_to, k=2)

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

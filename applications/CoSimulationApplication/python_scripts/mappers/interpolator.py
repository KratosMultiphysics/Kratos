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
            if len(self.directions) > 3:
                raise ValueError(f'too many directions given')

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
            self.tree = cKDTree(self.coords_from)
        else:  # less stable
            self.tree = cKDTree(self.coords_from, balanced_tree=False)
        self.distances, self.nearest = self.tree.query(self.coords_to, k=self.n_nearest)
        self.nearest = self.nearest.reshape(-1, self.n_nearest)

        # check for duplicate points
        self.check_duplicate_points()


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

    def check_bounding_box(self):
        # *** test this function

        # *** implement function


        # *** call in Initialize
        pass

    def check_duplicate_points(self):
        # checks only from-points  *** because tree is available, to-points are checked in interpolation in other direction anyway
        # *** test this function

        # calculate reference distance (diagonal of bounding box)
        diagonal = np.zeros(len(self.directions))
        for i in range(diagonal.size):
            diagonal[i] = self.coords_from[:, i].max() - self.coords_from[:, i].min()
        d_ref = np.linalg.norm(diagonal)
        if d_ref == 0.:
            raise ValueError('all from-points coincide')

        # check for duplicate points
        dist, _ = self.tree.query(self.coords_from, k=2)
        duplicate = (dist[:, 1] / d_ref < 1e-15)
        if duplicate.any():
            raise ValueError(f'{np.sum(duplicate)} duplicate points found in from-points, ' +
                             f'first duplicate = {self.coords_from[np.argmax(duplicate)]}')

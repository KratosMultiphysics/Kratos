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
        self.balanced_tree = self.settings['balanced_tree'].GetBool()  # *** optional parameter?
        self.n_nearest = 0  # must be set in sub-class!

        # get list with directions
        self.directions = []
        for direction in self.settings['directions'].GetArray():
            if direction not in ['X', 'Y', 'Z']:
                raise ValueError(f'"{direction}" is not a valid direction.')
            self.directions.append(direction + '0')
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

        # check bounding boxes  *** best position in code for this?
        self.check_bounding_box(model_part_from, model_part_to)

        # build and query tree
        if self.balanced_tree:  # time-intensive
            self.tree = cKDTree(self.coords_from)
        else:  # less stable
            self.tree = cKDTree(self.coords_from, balanced_tree=False)
        self.distances, self.nearest = self.tree.query(self.coords_to, k=self.n_nearest)
        self.nearest = self.nearest.reshape(-1, self.n_nearest)

        # check for duplicate points
        self.check_duplicate_points(model_part_from)


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

    def check_bounding_box(self, model_part_from, model_part_to):
        # set tolerances  #*** overwrite tolerances from parameters?
        tol_center_warning = 0.02
        tol_center_error = 0.1
        tol_minmax_warning = 0.1
        tol_minmax_error = 0.3

        # get bounding boxes
        coords_from = self.coords_from
        coords_to = self.coords_to
        from_min = coords_from.min(axis=0)
        from_max = coords_from.max(axis=0)
        to_min = coords_to.min(axis=0)
        to_max = coords_to.max(axis=0)
        from_center = (from_min + from_max) / 2
        to_center = (to_min + to_max) / 2

        # get reference distance (= average length of bounding box diagonal)
        diag_from = from_max - from_min
        diag_to = to_max - to_min
        d_ref = np.linalg.norm((diag_from + diag_to) / 2)

        # calculate errors on bounding boxes
        error_center = np.linalg.norm(from_center - to_center) / d_ref
        error_min = np.linalg.norm(from_min - to_min) / d_ref
        error_max = np.linalg.norm(from_max - to_max) / d_ref

        # raise warning or error if necessary
        msg_1 = f'ModelParts "{model_part_from.Name}", "{model_part_to.Name}": '
        msg_2 = ' values differ by '

        msg = f'{msg_1}center{msg_2}{100 * error_center:.1f}%'
        if error_center > tol_center_error:
            raise ValueError(msg)
        if error_center > tol_center_warning:
            raise Warning(msg)

        msg = f'{msg_1}min{msg_2}{100 * error_min:.1f}%'
        if error_min > tol_minmax_error:
            raise ValueError(msg)
        if error_min > tol_minmax_warning:
            raise Warning(msg)

        msg = f'{msg_1}max{msg_2}{100 * error_max:.1f}%'
        if error_max > tol_minmax_error:
            raise ValueError(msg)
        if error_max > tol_minmax_warning:
            raise Warning(msg)

    def check_duplicate_points(self, model_part_from):
        # checks only from-points  *** because tree is available, to-points are checked in interpolation in other direction anyway
        # *** test this function
        tol_warning = 1e-8
        tol_error = 1e-12

        # calculate reference distance (diagonal of bounding box)
        d_ref = np.linalg.norm(self.coords_from.max(axis=0) - self.coords_from.min(axis=0))
        if d_ref == 0.:
            raise ValueError('all from-points coincide')

        # check for duplicate points
        dist, _ = self.tree.query(self.coords_from, k=2)

        msg_1 = f'ModelPart {model_part_from.Name}: '
        msg_2 = f' duplicate points found, first duplicate: '

        duplicate = (dist[:, 1] / d_ref < tol_error)
        if duplicate.any():
            raise ValueError(f'{msg_1}{np.sum(duplicate)}{msg_2}' +
                             f'{self.coords_from[np.argmax(duplicate)]}.')

        duplicate = (dist[:, 1] / d_ref < tol_warning)
        if duplicate.any():
            raise Warning(f'{msg_1}{np.sum(duplicate)}{msg_2}' +
                          f'{self.coords_from[np.argmax(duplicate)]}.')

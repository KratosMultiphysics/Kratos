from scipy.spatial import cKDTree
import numpy as np
from multiprocessing import Pool

import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools
cs_data_structure = cs_tools.cs_data_structure

import time
from contextlib import contextmanager
@contextmanager
def timer(name=None, t=0, n=0, ms=False):
    startTime = time.time()
    yield
    elapsedTime = time.time() - startTime
    if ms:
        s = '\n' * n + '\t' * t + f'{elapsedTime * 1000:.2f}ms'
        s.replace(',', ' ')
    else:
        s = '\n' * n + '\t' * t + f'{elapsedTime:.1f}s'
    if name is not None:
        s += f' - {name}'
    s += '\n' * n
    print(s)

def Create(parameters):
    return MapperLinear3D(parameters)


# Class MapperLinear: Linear interpolation in 3D.
class MapperLinear3D(object):
    def __init__(self, parameters):
        """
        This mapper uses the 3 nearest points.
        Barycentric interpolation is done, based on
        the projection of the point onto the plane
        through the 3 points.
        If this projected point lies outside the
        triangle, the data is extrapolated.

        If the 3 closest points are co-linear,
        a (specified) number of next nearest points
        are checked to find a triangle.

        If no triangle is found, the point is
        projected on the line through the 2 nearest
        points, and interpolation/extrapolation
        is done based on only 2 points.

        Parallellization can be used for calculating
        the coefficients. Speed improvement on a
        36-core machine is only ~4.4, so there is
        probably room for improvement here.

        On 36 cores, 1e6 points took around
        70s to Initialize, and 10s to do the
        mapping with __call__.
        """
        super().__init__()

        self.settings = parameters['settings']
        self.interpolator = True
        self.balanced_tree = self.settings['balanced_tree'].GetBool()
        self.parallel = self.settings['parallel'].GetBool()

    def Initialize(self, model_part_from, model_part_to):
        self.n_from = model_part_from.NumberOfNodes()
        self.coords_from = np.zeros((self.n_from, 3))
        for i, node in enumerate(model_part_from.Nodes):
            self.coords_from[i, :] = [node.X0, node.Y0, node.Z0]

        self.n_to = model_part_to.NumberOfNodes()
        self.coords_to = np.zeros((self.n_to, 3))
        for i, node in enumerate(model_part_to.Nodes):
            self.coords_to[i, :] = [node.X0, node.Y0, node.Z0]

        # build and query tree
        if self.balanced_tree:  # time-intensive
            tree = cKDTree(self.coords_from)
        else:  # less stable
            tree = cKDTree(self.coords_from, balanced_tree=False)
        _, self.nearest_bis = tree.query(self.coords_to, k=20)
        self.nearest = self.nearest_bis[:, :3].copy()

        # calculate coefficients
        with timer('coeffs', ms=True):
            if self.parallel:
                with Pool() as self.pool:
                    out = self.pool.map(self.fun, range(self.n_to))
                self.coeffs = np.vstack(tuple([_[0] for _ in out]))
                self.nearest = np.vstack(tuple([_[1] for _ in out]))
            else:
                self.coeffs = np.zeros((self.n_to, 3))
                for i in range(self.n_to):
                    self.coeffs[i, :] = self.fun(i)[0]

    def Finalize(self):
        pass

    def __call__(self, args_from, args_to):
        # exactly the same as linear_1d

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

    def __getstate__(self):
        # necessary for multiprocessing
        self_dict = self.__dict__.copy()
        del self_dict['pool']
        return self_dict

    def fun(self, i):
        """
        indices:
            0 = point coords_to[i]
            1, 2, 3 = 3 nearest points from coords_from
            no index = projection of point 0 on
                       triangle 1-2-3 or line 1-2

        names:
            p = point
            v_ab = vector from point a to point b
            c = interpolation coefficients
        """
        def normalize(v):
            # return normalized a vector
            n = v / np.linalg.norm(v)
            return n

        def degenerate(p_1, p_2, p_3):
            # check if triangle is degenerate, return boolean
            v_12 = p_2 - p_1
            v_13 = p_3 - p_1
            v_23 = p_3 - p_2
            area_rect = np.linalg.norm(np.cross(v_12, v_13))
            if area_rect == 0.:
                return True

            length = max([np.linalg.norm(v_12), np.linalg.norm(v_13), np.linalg.norm(v_23)])
            aspect_ratio = length ** 2 / area_rect

            if aspect_ratio > 30:
                return True
            return False

        p_0 = self.coords_to[i, :]
        p_1 = self.coords_from[self.nearest[i, 0], :]
        p_2 = self.coords_from[self.nearest[i, 1], :]
        p_3 = self.coords_from[self.nearest[i, 2], :]

        triangle_interpolation = True
        if degenerate(p_1, p_2, p_3):
            triangle_interpolation = False
            for j in range(3, self.nearest_bis.shape[1]):
                p_3 = self.coords_from[self.nearest_bis[i, j], :]
                if not degenerate(p_1, p_2, p_3):
                    triangle_interpolation = True
                    self.nearest[i, 2] = self.nearest_bis[i, j]
                    break

        v_01 = p_1 - p_0
        v_12 = p_2 - p_1
        v_13 = p_3 - p_1
        v_23 = p_3 - p_2

        if triangle_interpolation:
            # project on triangle
            tmp = np.cross(v_12, v_13)
            v_n = normalize(tmp)
            p = p_0 + v_n * np.dot(v_01, v_n)

            # calculate weights
            v_p1 = p_1 - p
            v_p2 = p_2 - p
            v_p3 = p_3 - p
            ref = np.linalg.norm(tmp)

            v_a = np.cross(v_p2, v_p3)
            v_b = np.cross(v_12, v_13)
            c_1 = np.sign(np.dot(v_a, v_b)) * np.linalg.norm(v_a) / ref

            v_a = np.cross(v_p3, v_p1)
            v_b = np.cross(v_23, -v_12)
            c_2 = np.sign(np.dot(v_a, v_b)) * np.linalg.norm(v_a) / ref

            c_3 = 1. - c_1 - c_2
        else:
            # project on line
            v_t = normalize(v_12)
            p = p_0 + (v_01 - v_t * np.dot(v_01, v_t))

            # calculate weights
            v_p2 = p_2 - p
            ref = np.linalg.norm(v_12)
            c_1 = np.sign(np.dot(v_p2, v_12)) * np.linalg.norm(v_p2) / ref
            c_2 = 1. - c_1
            c_3 = 0.

        return (np.array([[c_1, c_2, c_3]]), self.nearest[i:i+1, :])
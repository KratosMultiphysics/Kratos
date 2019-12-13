from scipy.spatial import cKDTree

import numpy as np

import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools
from KratosMultiphysics.CoSimulationApplication.mappers.nearest import MapperNearest
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
    return MapperLinear2D(parameters)


# Class MapperLinear: Linear interpolation in 2D.
class MapperLinear2D(MapperNearest):
    def __init__(self, parameters):
        """
        This mapper uses the two closest points; this does
        not mean that the projection of the point to be
        interpolated lies in between them, e.g. when
        there is grid stretching.
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
        self.coord1 = self.settings['direction_1'].GetString().upper() + '0'
        self.coord2 = self.settings['direction_2'].GetString().upper() + '0'
        for coord in [self.coord1, self.coord2]:
            if coord not in ['X0', 'Y0', 'Z0']:
                raise ValueError(f'{coord[:-1]} is not a valid direction.')

    def Initialize(self, model_part_from, model_part_to):

        self.n_from = model_part_from.NumberOfNodes()
        coords_from = np.zeros((self.n_from, 2))
        for i, node in enumerate(model_part_from.Nodes):
            coords_from[i, 0] = getattr(node, self.coord1)
            coords_from[i, 1] = getattr(node, self.coord2)

        self.n_to = model_part_to.NumberOfNodes()
        coords_to = np.zeros((self.n_to, 2))
        for i, node in enumerate(model_part_to.Nodes):
            coords_to[i, 0] = getattr(node, self.coord1)
            coords_to[i, 1] = getattr(node, self.coord2)

        # build and query tree
        if self.balanced_tree:  # time-intensive
            tree = cKDTree(coords_from)
        else:  # less stable
            tree = cKDTree(coords_from, balanced_tree=False)
        _, self.nearest = tree.query(coords_to, k=2)

        # linear interpolation/extrapolation from nearest points
        def normalize(v):
            n = v / np.linalg.norm(v)
            return n

        # loop over all points (could be parallellized...)
        self.coeffs = np.zeros((self.n_to, 2))
        for i in range(self.n_to):
            p_0 = coords_to[i, :]  # point to be projected
            p_1 = coords_from[self.nearest[i, 0], :]  # nearest point 1
            p_2 = coords_from[self.nearest[i, 1], :]  # nearest point 2
            v_01 = p_1 - p_0
            v_12 = p_2 - p_1  # vector from 1 to 2
            v_t = normalize(v_12)  # tangential vector
            v_para = v_t * np.dot(v_t, v_01)
            v_perp = v_01 - v_para
            p = p_0 + v_perp  # projected point

            x, x_1, x_2 = p[0], p_1[0], p_2[0]
            c_1 = (x - x_2) / (x_1 - x_2)
            c_2 = (x - x_1) / (x_2 - x_1)
            self.coeffs[i, :] = [c_1, c_2]

            error = np.abs(1 - (c_1 + c_2))
            if error > 1e-8:
                raise ValueError(f'sum of coefficients is not 1 (error = {error:.1e})')

from scipy.spatial import cKDTree, distance
from scipy.linalg import solve
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
    return MapperRadialBasis2D(parameters)


# Class MapperLinear: Radial basis function interpolation in 2D.
class MapperRadialBasis2D(MapperNearest):
    def __init__(self, parameters):
        """
        TODO
        interpolation is done on lines

        The __call__ and Finalize methods are
        inherited from MapperNearest.


        HOW TO GET THIS THINGS TO WORK?
        --> do much more tests first, in crap.py file

        - 1D interpolation:
            - low, fixed number of from-points, start with uniform spacing
            - draw simple function through points, e.g. linear
            - get radial basis coeffs for given function values,
            draw interpolated function
            - now calculate the other way, get coeffs for fixed to-points,
            then add function value and see if it's the same

            ..



        at the end: write jupyter notebook


        """

        self.settings = parameters['settings']
        self.interpolator = True
        self.balanced_tree = self.settings['balanced_tree'].GetBool()
        self.coord1 = self.settings['direction_1'].GetString().upper() + '0'
        self.coord2 = self.settings['direction_2'].GetString().upper() + '0'
        for coord in [self.coord1, self.coord2]:
            if coord not in ['X0', 'Y0', 'Z0']:
                raise ValueError(f'{coord[:-1]} is not a valid direction.')
        self.ref_distance = None

    def Initialize(self, model_part_from, model_part_to):
        # *** check if there are enough nodes on from-side

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
        distances, self.nearest = tree.query(coords_to, k=9)

        self.coeffs = np.zeros_like(self.nearest)
        # loop over all to-points
        for index in range(self.n_to):
            d_ref = distances[index, -1]
            nearest = self.nearest[index, :]

            # create vector Phi_to, based on distances to from-points
            d_to = distances[index, :].reshape(-1, 1)
            Phi_to = self.phi(d_to, d_ref)

            # create matrix Phi, based on distances between from-points
            d = distance.squareform(distance.pdist(coords_from[nearest, :]))
            Phi = self.phi(d, d_ref)

            # solve system Phi^T c = Phi_t for c (Phi is symmetric)
            c = solve(Phi, Phi_to, sym_pos=True)  # , sym_pos=True

            # store c in coeffs
            self.coeffs[index] = c.flatten()

            print(f'sum of coeffs = {np.sum(c)}')  # *** check if coeffs add up to 1

    def phi(self, d, d_ref):
        r = d / d_ref
        return (1 - r) ** 4 * (1 + 4 * r) * np.heaviside(1 - r, 0)

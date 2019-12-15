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
    return MapperRadialBasis(parameters)


# Class MapperRadialBasis: Radial basis function interpolation in 2D.
class MapperRadialBasis(MapperNearest):
    def __init__(self, parameters):
        """
        TODO
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
            - ...
        - ...
        - at the end: make generalized function for n-dimensions,
        because most of the code would overlap anyways

        """
        # store settings
        self.settings = parameters['settings']
        self.interpolator = True
        self.balanced_tree = self.settings['balanced_tree'].GetBool()

        # get list with directions
        self.directions = []
        for direction in self.settings['directions'].list():
            tmp = direction.GetString().upper()
            if tmp not in ['X', 'Y', 'Z']:
                raise ValueError(f'"{tmp}" is not a valid direction.')
            self.directions.append(tmp + '0')

        # determine number of nearest neighbours
        if len(self.directions) == 3:
            self.n_nearest = 81
        else:
            self.n_nearest = 9

    def Initialize(self, model_part_from, model_part_to):
        self.n_from = model_part_from.NumberOfNodes()
        if self.n_from < self.n_nearest:
            raise ValueError('not enough from-points for radial basis interpolation')
        coords_from = np.zeros((self.n_from, len(self.directions)))
        for i, node in enumerate(model_part_from.Nodes):
            for j, direction in enumerate(self.directions):
                coords_from[i, j] = getattr(node, direction)

        self.n_to = model_part_to.NumberOfNodes()
        coords_to = np.zeros((self.n_to, len(self.directions)))
        for i, node in enumerate(model_part_to.Nodes):
            for j, direction in enumerate(self.directions):
                coords_to[i, j] = getattr(node, direction)

        # build and query tree
        if self.balanced_tree:  # time-intensive
            tree = cKDTree(coords_from)
        else:  # less stable
            tree = cKDTree(coords_from, balanced_tree=False)
        distances, self.nearest = tree.query(coords_to, k=self.n_nearest)

        self.coeffs = np.zeros(self.nearest.shape)  #*** problem was here: zeros_like made ndarray of ints...
        # loop over all to-points
        for i_to in range(self.n_to):
            d_ref = distances[i_to, -1] * 2
            nearest = self.nearest[i_to, :]

            # create vector Phi_to, based on distances to from-points
            d_to = distances[i_to, :].reshape(-1, 1)
            Phi_to = self.phi(d_to/ d_ref)

            # create matrix Phi, based on distances between from-points
            d = distance.squareform(distance.pdist(coords_from[nearest, :]))
            Phi = self.phi(d / d_ref)

            # solve system Phi^T c = Phi_t for c (Phi is symmetric)
            c = solve(Phi, Phi_to, sym_pos=True)  # , sym_pos=True

            # store c in coeffs
            self.coeffs[i_to, :] = c.flatten().copy()

            print(f'sum of coeffs = {np.sum(c)}')  # *** check if coeffs add up to 1

    def phi(self, r):
        return (1 - r) ** 4 * (1 + 4 * r) * np.heaviside(1 - r, 0)

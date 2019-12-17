from scipy.spatial import cKDTree, distance
from scipy.linalg import solve
import numpy as np
from multiprocessing import Pool, cpu_count


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

        """
        # store settings
        self.settings = parameters['settings']
        self.interpolator = True
        self.balanced_tree = self.settings['balanced_tree'].GetBool()
        self.parallel = self.settings['parallel'].GetBool()

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
        self.coords_from = np.zeros((self.n_from, len(self.directions)))
        for i, node in enumerate(model_part_from.Nodes):
            for j, direction in enumerate(self.directions):
                self.coords_from[i, j] = getattr(node, direction)

        self.n_to = model_part_to.NumberOfNodes()
        self.coords_to = np.zeros((self.n_to, len(self.directions)))
        for i, node in enumerate(model_part_to.Nodes):
            for j, direction in enumerate(self.directions):
                self.coords_to[i, j] = getattr(node, direction)

        if self.n_from < self.n_nearest:
            raise ValueError('not enough from-points for radial basis interpolation')

        # build and query tree
        with timer('tree', ms=True, t=1):
            if self.balanced_tree:  # time-intensive
                tree = cKDTree(self.coords_from)
            else:  # less stable
                tree = cKDTree(self.coords_from, balanced_tree=False)
            self.distances, self.nearest = tree.query(self.coords_to, k=self.n_nearest)

        with timer('coeffs', ms=True):
            iterable = []
            for i_to in range(self.n_to):
                nearest = self.nearest[i_to, :]
                iterable.append((self.distances[i_to, :],
                                 self.coords_from[nearest, :]))

            if self.parallel:
                processes = cpu_count()
                with Pool(processes=processes) as pool:
                    out = pool.starmap(get_coeffs, iterable)  # optimal chunksize automatically calculated
                self.coeffs = np.vstack(tuple(out))
            else:
                self.coeffs = np.zeros(self.nearest.shape)
                for tup in iterable:
                    self.coeffs[i_to, :] = get_coeffs(*tup).flatten()

        # *** warning: coeffs don't always add up to 1!

def get_coeffs(distances, coords_from):
    def phi(r):
        return (1 - r) ** 4 * (1 + 4 * r) * np.heaviside(1 - r, 0)

    d_ref = distances[-1] * 2

    # create column Phi_to, based on distances to from-points
    d_to = distances.reshape(-1, 1)
    Phi_to = phi(d_to / d_ref)

    # create matrix Phi, based on distances between from-points
    d = distance.squareform(distance.pdist(coords_from))
    Phi = phi(d / d_ref)

    # solve system Phi^T c = Phi_t for c (Phi is symmetric)
    coeffs = solve(Phi, Phi_to, sym_pos=True)

    return coeffs.reshape(1, -1)

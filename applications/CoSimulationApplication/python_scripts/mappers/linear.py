from scipy.spatial import cKDTree
import numpy as np
from multiprocessing import Pool, cpu_count

import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools
from KratosMultiphysics.CoSimulationApplication.mappers.interpolator import MapperInterpolator
cs_data_structure = cs_tools.cs_data_structure


def Create(parameters):
    return MapperLinear(parameters)

# *** should be adapted! doesn't work right now

# Class MapperLinear: linear interpolation.
class MapperLinear(MapperInterpolator):
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

        # *** TO DO: update this, put part in mappers.md
        """
        super().__init__(parameters)

        # store settings
        self.parallel = self.settings['parallel'].GetBool()

        # determine number of nearest neighbours *** ADAPT:1D, 2D, 3D
        if len(self.directions) == 3:
            self.n_nearest = 3
        else:
            self.n_nearest = 2

    def Initialize(self, model_part_from, model_part_to):
        super().Initialize(model_part_from, model_part_to)

        n_directions = len(self.directions)
        if n_directions == 1:
            get_coeffs = get_coeffs_1d
        elif n_directions == 2:
            get_coeffs = get_coeffs_2d
        else:
            get_coeffs = get_coeffs_3d

        # calculate coefficients  *** new
        with cs_tools.quicktimer('coeffs', ms=True):
            iterable = []
            for i_to in range(self.n_to):
                nearest = self.nearest[i_to, :]
                iterable.append((self.coords_from[nearest, :], self.coords_to[i_to, :]))

            if self.parallel:
                processes = cpu_count()
                with Pool(processes=processes) as pool:
                    out = pool.starmap(get_coeffs, iterable)  # optimal chunksize automatically calculated
                self.coeffs = np.vstack(tuple(out))
            else:
                self.coeffs = np.zeros(self.nearest.shape)
                for i_to, tup in enumerate(iterable):
                    self.coeffs[i_to, :] = get_coeffs(*tup).flatten()


def get_coeffs_1d(coords_from, coord_to):
    """
    conventions:
        P_a = point a

    subscripts:
        0 = to-point
        1 = from-point 1
        2 = from-point 2
    """
    P_0 = coord_to[0]
    P_1 = coords_from[0, 0]
    P_2 = coords_from[1, 0]

    # check if to-point lies between from-points
    coeffs = np.zeros(2)
    if (P_1 - P_0) * (P_2 - P_0) <= 0:
        coeffs[0] = np.abs((P_2 - P_0) / (P_2 - P_1))
    else:
        # nearest neighbour
        coeffs[0] = 1.
    coeffs[1] = 1. - coeffs[0]

    return coeffs.reshape(1, -1)

def get_coeffs_2d(coords_from, coord_to):
    """
    conventions:
        P_a = point a
        v_ab = vector from a to b

    subscripts:
        0 = to-point
        1 = from-point 1
        2 = from-point 2
        p = to-point projected line through from-points
        n = normalized vector
    """
    P_0 = coord_to
    P_1 = coords_from[0, :]
    P_2 = coords_from[1, :]

    v_01 = P_1 - P_0
    v_12n = normalize(P_2 - P_1)

    P_p = P_0 + (v_01 - v_12n * np.dot(v_01, v_12n))

    # check if to-point lies between from-points
    coeffs = np.zeros(2)
    if np.dot(P_2 - P_1, P_p - P_1) >= 0:
        # linear interpolation
        coeffs[0] = np.linalg.norm(P_2 - P_p) / np.linalg.norm(P_2 - P_1)
    else:
        # nearest neighbour
        coeffs[0] = 1.
    coeffs[1] = 1. - coeffs[0]

    return coeffs.reshape(1, -1)

def get_coeffs_3d():
    coeffs = np.zeros(3)  # *** TO DO
    return coeffs.reshape(1, -1)

def normalize(v):
    n = v / np.linalg.norm(v)
    return n

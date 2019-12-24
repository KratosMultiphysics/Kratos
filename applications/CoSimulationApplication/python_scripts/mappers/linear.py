import numpy as np
from multiprocessing import Pool, cpu_count

import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools
from KratosMultiphysics.CoSimulationApplication.mappers.interpolator import MapperInterpolator
cs_data_structure = cs_tools.cs_data_structure
import matplotlib.pyplot as plt

def Create(parameters):
    return MapperLinear(parameters)


# Class MapperLinear: linear interpolation.
class MapperLinear(MapperInterpolator):
    def __init__(self, parameters):
        """
        geometrical calculations:
            P_a = point a
            v_ab = vector from a to b

            subscripts:
                0 = to-point
                1 = 1st from-point
                2 = 2nd from-point
                2 = 3rd from-point
                p = to-point projected on line/triangle
                n = normal unit vector
                t = tangential unit vector
        """
        super().__init__(parameters)

        # store settings
        self.parallel = self.settings['parallel'].GetBool()

        # determine number of nearest neighbours
        if len(self.directions) == 3:
            self.n_nearest = 3
        else:
            self.n_nearest = 2

    def Initialize(self, model_part_from, model_part_to):
        super().Initialize(model_part_from, model_part_to)

        if len(self.directions) == 3:
            get_coeffs = get_coeffs_3d
        else:
            get_coeffs = get_coeffs_1d_2d

        # calculate coefficients
        with cs_tools.quicktimer('coeffs', ms=True):
            iterable = []
            for i_to in range(self.n_to):
                nearest = self.nearest[i_to, :]
                iterable.append((self.coords_from[nearest, :], self.coords_to[i_to, :]))

            if self.parallel:
                processes = cpu_count()
                with Pool(processes=processes) as pool:
                    # optimal chunksize automatically calculated
                    out = pool.starmap(get_coeffs, iterable)
                self.coeffs = np.vstack(tuple(out))
            else:
                self.coeffs = np.zeros(self.nearest.shape)
                for i_to, tup in enumerate(iterable):
                    self.coeffs[i_to, :] = get_coeffs(*tup).flatten()


def get_coeffs_1d_2d(coords_from, coord_to):
    coeffs = np.zeros(2)
    P_0 = coord_to
    P_1 = coords_from[0, :]
    P_2 = coords_from[1, :]

    coeffs[0] = line_interpolation_coeff(P_0, P_1, P_2)
    coeffs[1] = 1. - coeffs[0]
    return coeffs.reshape(1, -1)

def get_coeffs_3d(coords_from, coord_to):
    coeffs = np.zeros(3)
    P_0 = coord_to
    P_1 = coords_from[0, :]
    P_2 = coords_from[1, :]
    P_3 = coords_from[2, :]

    # check if triangle is degenerate (e.g. co-linear points)
    if degenerate_triangle(P_1, P_2, P_3):
        coeffs[0] = line_interpolation_coeff(P_0, P_1, P_2)
        coeffs[1] = 1. - coeffs[0]
    else:
        P_p = project_on_triangle(P_0, P_1, P_2, P_3)
        if point_on_triangle(P_p, P_1, P_2, P_3):
            # barycentric interpolation
            a_ref = triangle_area(P_1, P_2, P_3)
            coeffs[0] = triangle_area(P_p, P_2, P_3) / a_ref
            coeffs[1] = triangle_area(P_p, P_1, P_3) / a_ref
            coeffs[2] = 1. - coeffs[0] - coeffs[1]
        else:
            coeffs[0] = line_interpolation_coeff(P_0, P_1, P_2)
            coeffs[1] = 1. - coeffs[0]

    return coeffs.reshape(1, -1)

def line_interpolation_coeff(P_0, P_1, P_2):
    # project P_0 on line
    # *** only necessary if 2D??
    if P_0.size == 1:
        P_p = P_0
    else:
        v_01 = P_1 - P_0
        v_t = (P_2 - P_1) / np.linalg.norm(P_2 - P_1)
        P_p = P_0 + (v_01 - v_t * np.dot(v_01, v_t))

    # check if point lies on line
    if np.dot(P_1 - P_p, P_2 - P_p) <= 0:
        # linear interpolation
        c = np.linalg.norm(P_2 - P_p) / np.linalg.norm(P_2 - P_1)
    else:
        # nearest neighbour interpolation
        c = 1.
    return c

def degenerate_triangle(P_1, P_2, P_3):
    v_12 = P_2 - P_1
    v_13 = P_3 - P_1
    v_23 = P_3 - P_2

    # check if area of triangle is 0
    area_rect = np.linalg.norm(np.cross(v_12, v_13))
    if area_rect == 0.:
        return True

    # check if aspect ratio of bounding box is too high
    length = max([np.linalg.norm(v_12), np.linalg.norm(v_13), np.linalg.norm(v_23)])
    aspect_ratio = length ** 2 / area_rect
    if aspect_ratio > 30:
        return True
    return False

def project_on_triangle(P_0, P_1, P_2, P_3):
    v_n = np.cross(P_2 - P_1, P_3 - P_1)
    v_n /= np.linalg.norm(v_n)
    P_p = P_0 + v_n * np.dot(P_1 - P_0, v_n)
    return P_p

def point_on_triangle(P_p, P_1, P_2, P_3):
    a = np.cross(P_p - P_1, P_2 - P_1)
    b = np.cross(P_p - P_2, P_3 - P_2)
    c = np.cross(P_p - P_3, P_1 - P_3)
    # check if three vectors point in same direction
    if np.dot(a, b) >= 0 and np.dot(a, c) >= 0:
        return True
    return False

def triangle_area(P_1, P_2, P_3):
    return np.linalg.norm(np.cross(P_2 - P_1, P_3 - P_1) / 2)

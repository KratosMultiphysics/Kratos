from scipy.spatial import distance
from scipy.linalg import solve
import numpy as np
from multiprocessing import Pool, cpu_count


import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools
from KratosMultiphysics.CoSimulationApplication.mappers.interpolator import MapperInterpolator
cs_data_structure = cs_tools.cs_data_structure


def Create(parameters):
    return MapperRadialBasis(parameters)


# Class MapperRadialBasis: radial basis function interpolation.
class MapperRadialBasis(MapperInterpolator):
    def __init__(self, parameters):
        super().__init__(parameters)

        # store settings
        self.parallel = self.settings['parallel'].GetBool()

        # determine number of nearest neighbours
        if len(self.directions) == 3:
            self.n_nearest = 81
        else:
            self.n_nearest = 9

    def Initialize(self, model_part_from, model_part_to):
        super().Initialize(model_part_from, model_part_to)

        # calculate coefficients
        with cs_tools.quicktimer('coeffs', ms=True):
            iterable = []
            for i_to in range(self.n_to):
                nearest = self.nearest[i_to, :]
                iterable.append((self.distances[i_to, :],
                                 self.coords_from[nearest, :]))

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

def get_coeffs(distances, coords_from):
    def phi(r):
        return (1 - r) ** 4 * (1 + 4 * r) * np.heaviside(1 - r, 0)

    d_ref = distances[-1] * 2  # *** add parameter for this?

    # create column Phi_to, based on distances to from-points
    d_to = distances.reshape(-1, 1)
    Phi_to = phi(d_to / d_ref)

    # create matrix Phi, based on distances between from-points
    d = distance.squareform(distance.pdist(coords_from))
    Phi = phi(d / d_ref)

    # solve system Phi^T c = Phi_t for c (Phi is symmetric)
    coeffs = solve(Phi, Phi_to, sym_pos=True)

    return coeffs.reshape(1, -1)

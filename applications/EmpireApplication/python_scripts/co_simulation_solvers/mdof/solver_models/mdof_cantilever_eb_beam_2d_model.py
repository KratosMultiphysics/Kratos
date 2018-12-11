from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

# Importing the base class
from mdof_base_model import MDoFBaseModel
from co_simulation_tools import RecursivelyValidateAndAssignDefaults

# Other imports
import numpy as np
from scipy import linalg
from scipy.optimize import minimize
from functools import partial
import json
import os

# For numerically evaluating a symbolic input
from sympy import symbols
from sympy.core import sympify

def CreateSolver(model_settings):
    return MDoFCantileverEBBeam2DModel(model_settings)

class MDoFCantileverEBBeam2DModel(MDoFBaseModel):
    """
    A multi-degree-of-freedom MDoF model assuming
    bending-type deformations using the Euler-Bernoulli
    beam theory.

    ATTENTION:
    For this model a homogenous distribution of mass,
    stiffness and damping is a premise. For other cases
    this model is not adequate and changes need to be done.
    """
    def __init__(self, model_settings):

        default_settings = {
                "type": "cantilever_eb_beam_2d",
                "system_parameters":{
                    "density"           : 5.0,
                    "area"              : 10.0,
                    "target_frequency"  : 1.0,
                    "target_mode"       : 1,
                    "damping_ratio"     : 0.05,
                    "level_height"      : 3.5,
                    "number_of_levels"  : 3
                },
                "initial_conditions":{
                    "displacement"  : "none",
                    "velocity"      : "none",
                    "acceleration"  : "none",
                    "external_force": "none"
                }
            }

        RecursivelyValidateAndAssignDefaults(default_settings, model_settings)

        rho = model_settings["system_parameters"]["density"]
        area = model_settings["system_parameters"]["area"]
        target_freq = model_settings["system_parameters"]["target_frequency"]
        # adjust index
        target_mode = model_settings["system_parameters"]["target_mode"] - 1
        zeta = model_settings["system_parameters"]["damping_ratio"]
        level_height = model_settings["system_parameters"]["level_height"]
        num_of_levels = model_settings["system_parameters"]["number_of_levels"]

        self.m = self._CalculateMass(rho, area, level_height, num_of_levels)
        self.k = self._CalculateStiffness(self.m, level_height, num_of_levels, target_freq, target_mode)
        self.b = self._CalculateDamping(self.m, self.k, zeta)

        height_coordinates = self._GetNodalCoordinates(level_height, num_of_levels)

        self.nodal_coordinates = {'x0': np.zeros(len(height_coordinates)),
                             'y0': height_coordinates,
                             'x': None,
                             'y': None}

        initial_values = self._SetupInitialValues(model_settings['initial_conditions'],
                                                  self.nodal_coordinates['y0'])
        self.u0 = initial_values["displacement"]
        self.v0 = initial_values["velocity"]
        self.a0 = initial_values["acceleration"]
        self.f0 = initial_values["external_force"]

    def _GetNodalCoordinates(self, level_height, num_of_levels):
        nodal_coordinates = level_height * np.arange(1,num_of_levels+1)
        return nodal_coordinates

    def _CalculateMass(self, rho, area, level_height, num_of_levels):
        """
        Getting the consistant mass matrix
        """
        # mass values for one level
        length = level_height
        m_const = rho * area * length / 2
        m_elem = np.array([[1.0, 0.0],
                           [0.0, 1.0]])

        # global mass matrix initialization with zeros
        m_glob = np.zeros((num_of_levels + 1, num_of_levels + 1))
        # fill global mass matrix entries
        for i in range(num_of_levels):
            m_temp = np.zeros((num_of_levels +1, num_of_levels + 1))
            m_temp[i:i + 2, i:i + 2] = m_elem
            m_glob += m_const * m_temp

        # remove the fixed degrees of freedom -> applying Dirichlet BC implicitly
        for i in [0, 1]:
            m_glob = np.delete(m_glob, 0, axis=i)

        # return mass matrix
        return m_glob

    def _CalculateStiffness(self, m, level_height, num_of_levels, target_freq, target_mode):
        """
        Calculate uniform stiffness k_scalar. A uniform stiffness is assumed for all
        the elements and the value is calculated using an optimization (or "tuning")
        for a target frequency of a target mode.
        """
        print("Calculating stiffness k in MDoFBeamModel derived class \n")

        # setup k_scalar_guess as the input for the standard k for a shear-type
        # MDoF
        k_scalar_guess = 1000.

        # using partial to fix some parameters for the
        # self._calculate_frequency_for_current_scalar_k()
        optimizable_function = partial(self._CalculateFrequencyErrorForCurrentKScalar,
                                       m,
                                       level_height,
                                       num_of_levels,
                                       target_freq,
                                       target_mode)

        #print("Optimization for the target k matrix in MDoFBeamModel \n")
        minimization_result = minimize(optimizable_function,
                                       k_scalar_guess, method='Powell',
                                       options={'disp': False})
        #                               options={'disp': True})

        # returning only one value!
        k_scalar_opt = minimization_result.x

        return self._AssembleK(level_height, num_of_levels, k_scalar_opt)

    def _AssembleK(self, level_height, num_of_levels, k_scalar):
        """
        For the MDoFBeam model stiffness distribution according to beam theory is assumed
        the stiffness matrix is asembled with the k_scalar calculated.
        """
        # k_scalar = EI
        length = level_height
        k_const = k_scalar / pow(length, 3)
        # stifness values for one level
        k_elem = np.array([[        12,     6 * length,         -12,     6 * length],
                           [6 * length, 4 * length **2, -6 * length, 2 * length **2],
                           [       -12,    -6 * length,          12,    -6 * length],
                           [6 * length, 2 * length **2, -6 * length, 4 * length **2]])

        # global stiffness matrix initialization with zeros
        k_glob = np.zeros((2 * num_of_levels + 2, 2 * num_of_levels + 2))
        # fill global stiffness matix entries
        for i in range(num_of_levels):
            k_temp = np.zeros(
                (2 * num_of_levels + 2, 2 * num_of_levels + 2))
            k_temp[2 * i:2 * i + 4, 2 * i:2 * i + 4] = k_elem
            k_glob += k_const * k_temp

        # remove the fixed degrees of freedom
        for dof in [1, 0]:
            for i in [0, 1]:
                k_glob = np.delete(k_glob, dof, axis=i)

        # return stiffness matrix
        return k_glob

    def _CalculateDamping(self, m, k, zeta):
        """
        Calculate damping b based upon the Rayleigh assumption
        using the first 2 eigemodes - here generically i and i
        """
        print("Calculating damping b in MDoFShearModel derived class \n")
        mode_i = 0
        mode_j = 1
        zeta_i = zeta
        zeta_j = zeta

        # TODO: try to avoid this code duplication
        # raw results
        eig_values_raw, eigen_modes_raw = linalg.eigh(k, m)
        # rad/s
        eig_values = np.sqrt(np.real(eig_values_raw))
        # 1/s = Hz
        eig_freqs = eig_values / 2. / np.pi
        # sort eigenfrequencies
        eig_freqs_sorted_indices = np.argsort(eig_freqs)
        #

        a = np.linalg.solve(0.5 *
                            np.array(
                                [[1 / eig_values[eig_freqs_sorted_indices[mode_i]],
                                  eig_values[
                                  eig_freqs_sorted_indices[mode_i]]],
                                    [1 / eig_values[eig_freqs_sorted_indices[mode_j]],
                                     eig_values[
                                     eig_freqs_sorted_indices[
                                     mode_j]]]]),
                            [zeta_i, zeta_j])
        return a[0] * m + a[1] * k

    def _CalculateFrequencyErrorForCurrentKScalar(self, m, level_height, num_of_levels, target_freq, target_mode, k_scalar):
        k = self._AssembleK(level_height, num_of_levels, k_scalar)

        # TODO: try to avoid this code duplication
        # raw results
        eig_values_raw, eigen_modes_raw = linalg.eigh(k, m)
        # rad/s
        eig_values = np.sqrt(np.real(eig_values_raw))
        # 1/s = Hz
        eig_freqs = eig_values / 2. / np.pi
        # sort eigenfrequencies
        eig_freqs_sorted_indices = np.argsort(eig_freqs)
        #

        current_target_freq = eig_freqs[eig_freqs_sorted_indices[target_mode-1]]

        return (target_freq - current_target_freq) **2 / target_freq**2

    def _GetIOName(self):
        return "mdof_cantilever_eb_beam_2d_model"

    def _Name(self):
        return self.__class__.__name__

    # PMT: to be implemented
    def _DofList(self):
        '''
        A DoF list saying which DoF entry
        what kind of deformation it represents
        In this case probably:
        ["DeltaX","ThethaY","DeltaX","ThethaY",...]
        '''
        pass

    # PMT: to extend to be more rebust, generic, permit multiple variables
    # for now suitable for a '1d' (so line-type) model

    # wil NOT be ok for EB Beam - update, taken from Shear model
    def _SetupInitialValues(self, initial_values, height_coordinates):
        '''
        Loops over the prescribed values for
        displacement, velocity, acceleration, external load
        and generates an array depending on the input:
        'none' or a symbolic expression as a function of height (param x)
        '''
        for key, value in initial_values.items():
            if value.lower() == "none":
                value = np.zeros(len(height_coordinates))
            else:
                # evaluate numerically the symbolic function
                input_str = value
                x = symbols('x')
                sympified_expr = sympify(input_str)
                value = np.array([sympified_expr.evalf(subs={x:num_val}) for num_val in height_coordinates], dtype='float64')

            initial_values[key] = value

        return initial_values
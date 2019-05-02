from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

# Importing the base class
from mdof_base_model import MDoFBaseModel
from co_simulation_tools import RecursivelyValidateAndAssignDefaults

# Other imports
import numpy as np
import json
import os

def CreateModel(model_settings):
    return MDoFBridge2DoFModel(model_settings)

class MDoFBridge2DoFModel(MDoFBaseModel):
    """
    MDoF model - 2DoF for a bridge section

    ATTENTION:
    For this model it is assumed that
    inertia, stiffness and damping are decoupled
    """
    def __init__(self, model_settings):

        default_settings = {
                "type" : "bridge_2dof",
                "system_parameters":{
                    "length_of_section" : 0.75,
                    "help"              : "1st value - translational dof, 2nd value - rotational dof",
                    "mass_per_length"   : [5.0, 0.025],
                    "target_frequency"  : [2.75, 2],
                    "damping_log_decr"  : [0.05, 0.10]
                },
                "initial_conditions":{
                    "displacement"      : [0.05, 0.35],
                    "velocity"          : [0.0, 0.0],
                    "acceleration"      : [0.0, 0.0],
                    "external_force"    : [0.0, 0.0]
                }
            }

        RecursivelyValidateAndAssignDefaults(default_settings, model_settings)

        l = model_settings["system_parameters"]["length_of_section"]

        # heave (translation)
        # mass over unit length as input

        m = model_settings["system_parameters"]["mass_per_length"][0]
        m_h = m * l
        f_h = model_settings["system_parameters"]["target_frequency"][0]
        k_h = m_h * ((2 * np.pi * f_h) ** 2)
        logd_h = model_settings["system_parameters"]["damping_log_decr"][0]

        # pitch (rotation)
        # inertia over unit length as input
        I = model_settings["system_parameters"]["mass_per_length"][1]
        m_r = I * l
        f_r = model_settings["system_parameters"]["target_frequency"][1]
        k_r = m_r * ((2 * np.pi * f_r) ** 2)
        logd_r = model_settings["system_parameters"]["damping_log_decr"][1]

        self.m = self._CalculateMass(m_h, m_r)
        self.k = self._CalculateStiffness(k_h, k_r)
        self.b = self._CalculateDamping(m_h, m_r,
                                   k_h, k_r,
                                   logd_h, logd_r)

        # needed as placeholder
        self.nodal_coordinates = {"x0": None,
                            "y0": None,
                            "x": None,
                            "y": None}

        initial_values = self._SetupInitialValues(model_settings['initial_conditions'])
        self.u0 = initial_values["displacement"]
        self.v0 = initial_values["velocity"]
        self.a0 = initial_values["acceleration"]
        self.f0 = initial_values["external_force"]

    def _CalculateMass(self, m_h, m_r):
        """
        Calculate mass m
        """
        return np.array([[m_h, 0],
                         [0, m_r]])

    def _CalculateStiffness(self, k_h, k_r):
        """
        Calculate mass k
        """
        return np.array([[k_h, 0],
                         [0, k_r]])

    def _CalculateDamping(self, m_h, m_r, k_h, k_r, logd_h, logd_r):
        """
        Calculate damping b
        """
        xi_h = logd_h / (2 * np.pi)
        xi_r = logd_r / (2 * np.pi)

        return np.array([[2 * xi_h * np.sqrt(m_h * k_h), 0],
                         [0, 2 * xi_r * np.sqrt(m_r * k_r)]])

    def _GetIOName(self):
        #return "mdof_bridge_2dof_model"
        return "mdof"

    def _Name(self):
        return self.__class__.__name__

    # PMT: to be implemented
    def _DofList(self):
        '''
        A DoF list saying which DoF entry
        what kind of deformation it represents
        In this case probably:
        ["DeltaY","ThethaZ"]
        '''
        pass

    def _SetupInitialValues(self, initial_values):
        '''
        From a list generate numpy array for compatibility
        '''
        for key, value in initial_values.items():

            value = np.array([num_val for num_val in value], dtype='float64')

            initial_values[key] = value

        return initial_values
from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

# Importing the base class
from mdof_base_model import MDoFBaseModel
from co_simulation_tools import RecursivelyValidateAndAssignDefaults

# Other imports
import numpy as np
import json
import os

def CreateModel(model_settings):
    return MDoFSDoFModel(model_settings)

class MDoFSDoFModel(MDoFBaseModel):
    """
    A single-degree-of-freedom SDoF model

    Uses the MDoF solver
    """
    def __init__(self, model_settings):

        default_settings = {
                "type" : "sdof",
                "system_parameters":{
                    "mass"              : 1.0,
                    "target_frequency"  : 2.0,
                    "damping_ratio"     : 0.0
                },
                "initial_conditions":{
                    "displacement"  : 0.5,
                    "velocity"      : 0.0,
                    "acceleration"  : 0.0,
                    "external_force": 0.0
                }
            }

        RecursivelyValidateAndAssignDefaults(default_settings, model_settings)

        m = model_settings["system_parameters"]["mass"]
        target_freq = model_settings["system_parameters"]["target_frequency"]
        zeta = model_settings["system_parameters"]["damping_ratio"]

        k = self._CalculateStiffness(m, target_freq)
        b = self._CalculateDamping(m, k, zeta)

        # needed as placeholder
        self.nodal_coordinates = {"x0": None,
                            "y0": None,
                            "x": None,
                            "y": None}

        # mass, stiffness and damping matrices
        self.m = np.array([[m]])
        self.k = np.array([[k]])
        self.b = np.array([[b]])

        # initial conditions - displacement, velocity, acceleration, external force
        self.u0 = np.array([model_settings["initial_conditions"]["displacement"]])
        self.v0 = np.array([model_settings["initial_conditions"]["velocity"]])
        self.a0 = np.array([model_settings["initial_conditions"]["acceleration"]])
        self.f0 = np.array([model_settings["initial_conditions"]["external_force"]])

    def _CalculateStiffness(self, m, target_freq):
        """
        Calculate stiffness k
        """
        return m * (target_freq * 2 * np.pi)**2

    def _CalculateDamping(self, m, k, zeta):
        """
        Calculate damping b
        """
        return zeta * 2.0 * np.sqrt(m * k)

    def _GetIOName(self):
        return "mdof_sdof_model"

    def _Name(self):
        return self.__class__.__name__

    # PMT: to be implemented
    def _DofList(self):
        '''
        A DoF list saying which DoF entry
        what kind of deformation it represents
        In this case probably:
        ["DeltaX"]
        '''
        pass
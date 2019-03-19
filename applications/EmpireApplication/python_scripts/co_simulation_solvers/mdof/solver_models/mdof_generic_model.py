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
    A generic MDoF model

    Mass, stiffness, damping matrices and
    vectors for initial displacement, velocity, acceleration,
    force have to be passed explicitly

    DEFAULT VALUES: the example A.K. Chopra - Dynamics of structures
    Example 11.1 with Rayleigh damping, all DoFs represent a lateral
    displacement
    """
    def __init__(self, model_settings):

        default_settings = {
                "type" : "generic",
                "system_parameters":{
                    "mass"      : [[1.036, 0.0, 0.0], [0.0, 1.036, 0.0], [0.0, 0.0, 0.518]],
                    "stiffness" : [[1220.0, -610.0, 0.0], [-610.0, 1220.0, -610.0], [0.0, -610.0, 610.0]],
                    "damping"   : [[3.55, -1.3, 0.0], [-1.3, 3.55, -1.3], [0.0, -1.3, 1.77]],
                    "height_coordinates" : [3.5, 7.0, 10.5]
                },
                "initial_conditions":{
                    "displacement"  : [0.1, 0.15, 0.45],
                    "velocity"      : [0.01, 0.015, 0.045],
                    "acceleration"  : [0.02, 0.03, 0.09],
                    "external_force": [1.2, 2.3, 4.1]
                }
            }

        RecursivelyValidateAndAssignDefaults(default_settings, model_settings)

        # needed as placeholder
        self.nodal_coordinates = {"x0": np.array(model_settings["system_parameters"]["height_coordinates"]),
                            "y0": None,
                            "x": None,
                            "y": None}

        # mass, stiffness and damping matrices
        self.m = np.array(model_settings["system_parameters"]["mass"])
        self.k = np.array(model_settings["system_parameters"]["stiffness"])
        self.b = np.array(model_settings["system_parameters"]["damping"])

        # initial conditions - displacement, velocity, acceleration, external force
        self.u0 = np.array(model_settings["initial_conditions"]["displacement"])
        self.v0 = np.array(model_settings["initial_conditions"]["velocity"])
        self.a0 = np.array(model_settings["initial_conditions"]["acceleration"])
        self.f0 = np.array(model_settings["initial_conditions"]["external_force"])

    def _GetIOName(self):
        return "mdof_generic_model"

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
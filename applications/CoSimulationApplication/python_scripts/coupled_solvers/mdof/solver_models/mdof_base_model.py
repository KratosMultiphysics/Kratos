
from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

# Other imports
import numpy as np
import json
import os

# PMT create a nice base class, cleanup, check if needed

class MDoFBaseModel(object):
    '''
    Provide a description related to what the mdof model
    is useful for, what are the limitations, underlying assumption
    for mass, stiffness, damping, initial conditions etc.
    '''
    def __init__(self, model_settings):

        # provide and validate default setttings
        default_settings = None

        # create mass, stiffness and damping matrices
        self.m = self._CalculateMass()
        self.k = self._CalculateStiffness()
        self.b = self._CalculateDamping()

        # get nodal coordinates if needed
        height_coordinates = self._GetNodalCoordinates()
        self.nodal_coordinates = {'x0': None,
                             'y0': None,
                             'x': None,
                             'y': None}

        # initial initial value in the correct format
        initial_values = self._SetupInitialValues()
        self.u0 = None #initial_values["displacement"]
        self.v0 = None #initial_values["velocity"]
        self.a0 = None #initial_values["acceleration"]
        self.f0 = None #initial_values["external_force"]

    def _GetNodalCoordinates(self):
        """
        """
        pass

    def _CalculateMass(self):
        """
        """
        pass

    def _CalculateStiffness(self):
        """
        """
        pass

    def _CalculateDamping(self):
        """
        """
        pass

    def _Name(self):
        return self.__class__.__name__

    def _DofList(self):
        """
        """
        pass

    def _SetupInitialValues(self):
        """
        """
        pass
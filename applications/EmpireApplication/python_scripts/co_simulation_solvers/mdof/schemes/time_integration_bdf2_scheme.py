from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

# Importing the base class
from time_integration_base_scheme import TimeIntegrationBaseScheme

# Importing tools
from co_simulation_tools import ValidateAndAssignDefaults

# Other imports
import numpy as np
import json
import os

def CreateScheme(scheme_settings):
    return TimeIntegrationBDF2Scheme(scheme_settings)

# PMT to be checked, seems to be written acceleration based, might need correction

class TimeIntegrationBDF2Scheme(TimeIntegrationBaseScheme):
    """
    A single-degree-of-freedom SDoF model

    Using for testing of the MDoF solver
    """
    def __init__(self, scheme_settings):

        default_settings = {
                "type"      : "bdf2",
                "settings"  : {}
            }

        ValidateAndAssignDefaults(default_settings, scheme_settings)

        # add buffer size - this is not user-specified
        # each derived scheme specifies it
        scheme_settings["settings"].update({"buffer_size":4})

        # base scheme settings
        super(TimeIntegrationBDF2Scheme, self).__init__(scheme_settings["settings"])

    def _AssembleLHS(self, model):
        """
        """
        pass

    def _AssembleRHS(self, model):
        """
        """
        pass

    def UpdateDerivedValues(self):
        """
        """
        pass
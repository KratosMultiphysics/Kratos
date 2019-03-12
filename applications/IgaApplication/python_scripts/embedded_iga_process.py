from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import numpy as np
# Importing the Kratos Library
import KratosMultiphysics
from KratosMultiphysics.IgaApplication import *

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return EmbeddedIgaProcess(Model, settings["Parameters"])

class EmbeddedIgaProcess(KratosMultiphysics.Process):

    def __init__(self, model, params):
        KratosMultiphysics.Process.__init__(self)

        default_parameters = KratosMultiphysics.Parameters('''{
            "tessellation_tolerance"   : 0.0,
            "triangulation_tolerance"  : 0.0,
            "max_triangulation_iterations"  : 0.0
        }''')

        ## Overwrite the default settings with user-provided parameters
        self.params = params
        self.params.ValidateAndAssignDefaults(default_parameters)

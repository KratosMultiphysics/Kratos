from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# Importing the Kratos Library
import KratosMultiphysics
from KratosMultiphysics import eigen_solver_factory

import KratosMultiphysics.IgaApplication as IGA


def Factory(settings, model):
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    return NitscheStabilizationProcess(model, settings["Parameters"])

class NitscheStabilizationProcess(KratosMultiphysics.Process):
    """This class is used in order to compute automatically the Nitsche stabilization factor.

    Only the member variables listed below should be accessed directly.

    Public member variables:
    model -- the container of the different model parts.
    params -- Kratos parameters containing the settings.
    """
    def __init__(self, model, params):
        """ The default constructor of the class

        Keyword arguments:
        self -- It signifies an instance of a class.
        model -- the container of the different model parts.
        params -- Kratos parameters containing solver settings.
        """
        KratosMultiphysics.Process.__init__(self)
        self.process = IGA.AssignIntegrationPointsToBackgroundElementsProcess(model, params)
        self.model = model

    def ExecuteBeforeOutputStep(self):
        # Get the model parts which divide the problem
        self.process.AssignIntegrationPoints()

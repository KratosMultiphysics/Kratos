from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Other imports
import numpy as np
from numpy import linalg as la

class CouplingInterfaceResidual(object):
    """This class computes residuals and norms related to an InterfaceData
    Maybe this will be extended in the future towards different residuals
    """
    def __init__(self, interface_data):
        self.interface_data = interface_data
        self.UpdateReferenceData()

    def UpdateReferenceData(self):
        self.ref_data = self.interface_data.GetData()

    def GetResidual(self):
        return self.interface_data.GetData() - self.ref_data

    def GetNorm(self):
        residual = self.GetResidual()
        return la.norm(residual) / np.sqrt(residual.size)

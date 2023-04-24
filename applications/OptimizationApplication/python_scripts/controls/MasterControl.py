from abc import ABC, abstractmethod
import KratosMultiphysics as Kratos
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import SupportedSensitivityFieldVariableTypes

class MasterControl(object):

    def __init__(self) -> None:
        self.SubControlList = []

    def Comp
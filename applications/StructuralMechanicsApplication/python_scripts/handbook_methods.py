import KratosMultiphysics
import KratosMultiphysics.StructuralMechanicsApplication as SMA
import numpy as np

class HandbookMethods:

    def __init__(self):
        pass

    @staticmethod
    def von_mises(von_mises_stress: float, yield_stress: float):
        RF = yield_stress/von_mises_stress
        return RF
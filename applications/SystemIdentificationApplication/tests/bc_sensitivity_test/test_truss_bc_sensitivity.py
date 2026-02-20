import numpy as np
import KratosMultiphysics as Kratos
import KratosMultiphysics.SystemIdentificationApplication as KratosSI
import KratosMultiphysics.KratosUnittest as UnitTest
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import StructuralMechanicsAnalysis

class TestTrussBoundaryConditionSensitivity(UnitTest.TestCase):
    def test_TrussBoundaryCondition(self):
        with open("ProjectParameters.json", "r") as file_input:
            parameters = Kratos.Parameters(file_input.read())



import KratosMultiphysics

import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.ShallowWaterApplication.shallow_water_analysis import ShallowWaterAnalysis

class ShallowWaterTestFactory(KratosUnittest.TestCase):

    def test_execution(self):
        with KratosUnittest.WorkFolderScope(self.execution_directory, __file__):

            with open(self.execution_file + "_parameters.json",'r') as parameter_file:
                ProjectParameters = KratosMultiphysics.Parameters(parameter_file.read())

            model = KratosMultiphysics.Model()
            test = ShallowWaterAnalysis(model, ProjectParameters)
            test.Run()

class TestLagrangianShallowWaterElement(ShallowWaterTestFactory):
    execution_directory = "element_tests"
    execution_file = "lagrangian_swe"

class TestShallowWaterElement(ShallowWaterTestFactory):
    execution_directory = "element_tests"
    execution_file = "swe"

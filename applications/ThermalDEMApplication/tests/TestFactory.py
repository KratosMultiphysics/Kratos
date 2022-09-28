import os
import KratosMultiphysics
import KratosMultiphysics.DEMApplication
import KratosMultiphysics.ThermalDEMApplication
import KratosMultiphysics.KratosUnittest as KratosUnittest
from   KratosMultiphysics.ThermalDEMApplication.thermal_dem_analysis import ThermalDEMAnalysis

# This utility will control the execution scope
class controlledExecutionScope:
    def __init__(self, scope):
        self.currentPath = os.getcwd()
        self.scope = scope
    
    def __enter__(self):
        os.chdir(self.scope)
    
    def __exit__(self, type, value, traceback):
        os.chdir(self.currentPath)

# General test factory
class TestFactory(KratosUnittest.TestCase):
    
    def setUp(self):
        with controlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):
            # Read ProjectParameters
            with open(self.file_parameters,'r') as parameter_file:
                ProjectParameters = KratosMultiphysics.Parameters(parameter_file.read())
            
            # Create Model
            model = KratosMultiphysics.Model()
			
            self.test = ThermalDEMAnalysis(model,ProjectParameters)

    def test_execution(self):
        with controlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):
            self.test.Run()

    def tearDown(self):
        pass

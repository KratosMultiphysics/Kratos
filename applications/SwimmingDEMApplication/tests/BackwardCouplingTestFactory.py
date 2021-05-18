import os

# Importing the Kratos Library
import KratosMultiphysics

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as KratosUnittest
import tests_python_scripts.backward_coupling_scripts.backward_coupling_test_analysis as backward_coupling_test_analysis
BackwardCouplingTestAnalysis = backward_coupling_test_analysis.BackwardCouplingTestAnalysis

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
            # Setting parameters

            with open(self.file_parameters,'r') as parameter_file:
                self.parameters = KratosMultiphysics.Parameters(parameter_file.read())

            # Create Model
            self.model = KratosMultiphysics.Model()

            self.test = BackwardCouplingTestAnalysis(self.model, self.parameters)

    def test_execution(self):
        with controlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):
            self.test.Run()
        self._check_results()

    def tearDown(self):
        pass

    def _check_results(self):
        results = self.parameters['results']
        particles_component_volume = results['particles_component_volume'].GetDouble()
        total_particles_volume = results['total_particles_volume'].GetDouble()
        self.assertNotAlmostEqual(particles_component_volume, 0.0)
        self.assertAlmostEqual(particles_component_volume, total_particles_volume)

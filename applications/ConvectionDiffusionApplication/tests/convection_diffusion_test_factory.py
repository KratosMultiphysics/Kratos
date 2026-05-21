# Importing the Kratos Library
import KratosMultiphysics

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.ConvectionDiffusionApplication.convection_diffusion_analysis import ConvectionDiffusionAnalysis
import KratosMultiphysics.kratos_utilities as kratos_utilities

# Other imports
import os

# This utility will control the execution scope in case we need to access files or we depend
# on specific relative locations of the files.

# TODO: Should we move this to KratosUnittest?
class controlledExecutionScope:
    def __init__(self, scope):
        self.currentPath = os.getcwd()
        self.scope = scope

    def __enter__(self):
        os.chdir(self.scope)

    def __exit__(self, type, value, traceback):
        os.chdir(self.currentPath)


class ConvectionDiffusionTestFactory(KratosUnittest.TestCase):
    def setUp(self):
        # Within this location context:
        with controlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):

            # Reading the ProjectParameters
            with open(self.file_name + "_parameters.json",'r') as parameter_file:
                ProjectParameters = KratosMultiphysics.Parameters(parameter_file.read())

            self.modify_parameters(ProjectParameters)

            # To avoid many prints
            if (ProjectParameters["problem_data"]["echo_level"].GetInt() == 0):
                KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)

            # Creating the test
            model = KratosMultiphysics.Model()
            self.test = ConvectionDiffusionAnalysis(model, ProjectParameters)
            self.test.Initialize()

    def modify_parameters(self, project_parameters):
        """This function can be used in derived classes to modify existing parameters
        before the execution of the test (e.g. switch to MPI)
        """
        pass

    def test_execution(self):
        # Within this location context:
        with controlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):
            self.test.RunSolutionLoop()

    def tearDown(self):
        # Within this location context:
        with controlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):
            self.test.Finalize()

@KratosUnittest.skipIfApplicationsNotAvailable("LinearSolversApplication")
class BasicConvectionDiffusionStationaryTest(ConvectionDiffusionTestFactory):
    file_name = "basic_conv_diffusion_test/basic_conv_diffusion_test_stationary"

@KratosUnittest.skipIfApplicationsNotAvailable("LinearSolversApplication")
class BasicConvectionDiffusionTransientTest(ConvectionDiffusionTestFactory):
    file_name = "basic_conv_diffusion_test/basic_conv_diffusion_test_transient"

@KratosUnittest.skipIfApplicationsNotAvailable("LinearSolversApplication")
class BasicConvectionDiffusionTransientSemiImplicitTest(ConvectionDiffusionTestFactory):
    file_name = "basic_conv_diffusion_test/basic_conv_diffusion_test_transient_semi_implicit"

@KratosUnittest.skipIfApplicationsNotAvailable("LinearSolversApplication")
class BasicDiffusionStationaryTest(ConvectionDiffusionTestFactory):
    file_name = "basic_conv_diffusion_test/basic_diffusion_test_stationary"

@KratosUnittest.skipIfApplicationsNotAvailable("StructuralMechanicsApplication")
class SimpleThermoMechanicalTest(ConvectionDiffusionTestFactory):
    file_name = "thermo_mechanical_tests/thermo_mechanical/coupled_problem_test"

@KratosUnittest.skipIfApplicationsNotAvailable("StructuralMechanicsApplication")
class SimpleThermoMechanicalTableAccessorTest(ConvectionDiffusionTestFactory):
    file_name = "thermo_mechanical_tests/thermo_mechanical_accessor/coupled_accessor_test"

@KratosUnittest.skipIfApplicationsNotAvailable("StructuralMechanicsApplication")
class SimpleThermoMechanicalDamageTest(ConvectionDiffusionTestFactory):
    file_name = "thermo_mechanical_tests/thermo_mechanical_damage/thermo_mech_damage_test"
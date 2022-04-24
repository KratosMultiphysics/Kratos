import KratosMultiphysics as km
import KratosMultiphysics.TopologyOptimizationApplication
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utilities


class TopologyOptimizationTestFactory(KratosUnittest.TestCase):
    def setUp(self):
        km.Logger.GetDefaultOutput().SetSeverity(km.Logger.Severity.WARNING)

    def test_execution(self):
        with KratosUnittest.WorkFolderScope(
            self.execution_directory, __file__
        ):
            __import__(self.execution_directory + '.' + self.execution_file)

    def tearDown(self):
        with KratosUnittest.WorkFolderScope(
            self.execution_directory, __file__
        ):
            kratos_utilities.DeleteDirectoryIfExisting('__pycache__')


class Small_Cantilever_test(TopologyOptimizationTestFactory):
    execution_directory = 'Small_Cantilever'
    execution_file = 'run_TopOpt'


class Small_Cantilever_RAMP_test(TopologyOptimizationTestFactory):
    execution_directory = 'Small_Cantilever_RAMP'
    execution_file = 'run_TopOpt'


def AssembleTestSuites():
    suites = KratosUnittest.KratosSuites

    smallSuite = suites['small']
    smallSuite.addTest(Small_Cantilever_test('test_execution'))
    smallSuite.addTest(Small_Cantilever_RAMP_test('test_execution'))

    nightSuite = suites['nightly']
    nightSuite.addTests(smallSuite)

    validationSuite = suites['validation']
    validationSuite.addTests(nightSuite)

    allSuite = suites['all']
    allSuite.addTests(validationSuite)

    return suites


if __name__ == '__main__':
    KratosUnittest.runTests(AssembleTestSuites())

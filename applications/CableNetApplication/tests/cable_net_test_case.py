
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utils


class CableNetTestCase(KratosUnittest.TestCase):
    '''This class is the basis for the testing the framework
    It can be used to test complete cases with the "CoSimulation-Analysis"
    '''
    def _runTest(self):
        self.test.Run()
        kratos_utils.DeleteTimeFiles(self.problem_dir_name)
    # called only once for this class, opposed of tearDown()
    @classmethod
    def tearDownClass(cls):
        # Cleaning
        kratos_utils.DeleteDirectoryIfExisting("__pycache__")


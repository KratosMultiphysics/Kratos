
from KratosMultiphysics import *

from KratosMultiphysics import process_factory
from KratosMultiphysics.kratos_utilities import DeleteFileIfExisting
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.FSIApplication import fsi_analysis

class WorkFolderScope:
    def __init__(self, work_folder):
        self.currentPath = os.getcwd()
        self.scope = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)),work_folder))

    def __enter__(self):
        os.chdir(self.scope)

    def __exit__(self, type, value, traceback):
        os.chdir(self.currentPath)

class MokBenchmarkTest(KratosUnittest.TestCase):

    def setUp(self):
        self.print_output = False
        self.work_folder = "MokBenchmarkTest"
        self.settings = "ProjectParameters.json"
        self.fluid_input_file = "mok_benchmark_Fluid"
        self.structure_input_file = "mok_benchmark_Fluid"

    def tearDown(self):
        DeleteFileIfExisting(self.fluid_input_file+'.time')
        DeleteFileIfExisting(self.structure_input_file+'.time')

    def testMokBenchmark(self):
        with WorkFolderScope(self.work_folder):
            model = Model()
            parameter_file = open(self.settings, 'r')
            project_parameters = Parameters(parameter_file.read())
            fsi_analysis.FSIAnalysis(model, project_parameters).Run()

if __name__ == '__main__':
    test = MokBenchmarkTest()
    test.setUp()
    test.print_output = True
    test.testMokBenchmark()
    test.tearDown()

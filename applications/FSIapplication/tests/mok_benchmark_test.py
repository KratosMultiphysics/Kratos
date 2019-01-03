from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

from KratosMultiphysics import *
from KratosMultiphysics.FSIApplication import *
try:
    from KratosMultiphysics.MeshMovingApplication import *
    from KratosMultiphysics.FluidDynamicsApplication import *
    from KratosMultiphysics.ExternalSolversApplication import *
    from KratosMultiphysics.StructuralMechanicsApplication import *
except ImportError:
    pass

from os import remove

import process_factory
import KratosMultiphysics.KratosUnittest as KratosUnittest

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
        self.deleteOutFile(self.fluid_input_file+'.time')
        self.deleteOutFile(self.structure_input_file+'.time')

    def deleteOutFile(self,filename):
        with WorkFolderScope(self.work_folder):
            try:
                remove(filename)
            except FileNotFoundError as e:
                pass

    def testMokBenchmark(self):
        with WorkFolderScope(self.work_folder):
            model = Model()
            parameter_file = open(self.settings, 'r')
            project_parameters = Parameters(parameter_file.read())
            import fsi_analysis
            fsi_analysis.FSIAnalysis(model, project_parameters).Run()

if __name__ == '__main__':
    test = MokBenchmarkTest()
    test.setUp()
    test.print_output = True
    test.testMokBenchmark()
    test.tearDown()

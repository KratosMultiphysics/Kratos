import os
from KratosMultiphysics import *
from KratosMultiphysics.FluidDynamicsApplication import *
try:
    from KratosMultiphysics.HDF5Application import *
    have_hdf5 = True
except ImportError:
    have_hdf5 = False

import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utilities

from fluid_dynamics_analysis import FluidDynamicsAnalysis
from fluid_analysis_without_solution import FluidAnalysisWithoutSolution

class ControlledExecutionScope:
    def __init__(self, scope):
        self.currentPath = os.getcwd()
        self.scope = scope

    def __enter__(self):
        os.chdir(self.scope)

    def __exit__(self, type, value, traceback):
        os.chdir(self.currentPath)

@KratosUnittest.skipUnless(have_hdf5,"Missing required application: HDF5Application.")
class HDF5IOTest(KratosUnittest.TestCase):

    def setUp(self):
        self.work_folder = os.path.join(os.path.dirname(os.path.realpath(__file__)),"Cavity")

    def _removeH5Files(self, model_part_name):
        for name in os.listdir():
            if name.find(model_part_name) == 0:
                kratos_utilities.DeleteFileIfExisting(name)

    def testInputOutput(self):
        with ControlledExecutionScope(self.work_folder):
            # run simulation and write to hdf5 file
            model_out = Model()
            with open("output_test_parameters.json", 'r') as parameter_file:
                project_parameters = Parameters(parameter_file.read())
            test = FluidDynamicsAnalysis(model_out,project_parameters)
            test.Run()

            # start new simulation and read from hdf5 file
            model_in = Model()
            with open("input_test_parameters.json", 'r') as parameter_file:
                project_parameters = Parameters(parameter_file.read())
            test = FluidAnalysisWithoutSolution(model_in,project_parameters)
            test.Run()

    def tearDown(self):
        with ControlledExecutionScope(self.work_folder):
            # remove hdf5 file
            self._removeH5Files("MainModelPart")
            # remove other generated files
            kratos_utilities.DeleteFileIfExisting("./square5.time")
            kratos_utilities.DeleteFileIfExisting("./reference_results.json")

if __name__ == '__main__':
    KratosUnittest.main()


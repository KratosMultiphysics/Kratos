import os
from KratosMultiphysics import *
import KratosMultiphysics.KratosUnittest as KratosUnittest
import test_MainKratos

class ControlledExecutionScope:
    def __init__(self, scope):
        self.currentPath = os.getcwd()
        self.scope = scope

    def __enter__(self):
        os.chdir(self.scope)

    def __exit__(self, type, value, traceback):
        os.chdir(self.currentPath)

class TestCase(KratosUnittest.TestCase):

    def setUp(self):
        pass
    
    def _remove_file(self, file_path):
        if os.path.isfile(file_path):
            os.remove(file_path)

    def _remove_h5_files(self, model_part_name):
        for name in os.listdir():
            if name.find(model_part_name) == 0:
                self._remove_file(name)

    def test_Execution(self):
        with ControlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):
            # run simulation and write to hdf5 file
            parameter_file = open("test_input_output/output_test_parameters.json", 'r')
            project_parameters = Parameters(parameter_file.read())
            parameter_file.close()
            test = test_MainKratos.MainKratos(project_parameters)
            test.Solve()
            # start new simulation and read from hdf5 file
            parameter_file = open("test_input_output/input_test_parameters.json", 'r')
            project_parameters = Parameters(parameter_file.read())
            parameter_file.close()
            test = test_MainKratos.MainKratos(project_parameters)
            test.Solve()
            # remove hdf5 file
            self._remove_h5_files("MainModelPart")
            # remove other generated files
            self._remove_file("./test_input_output/io_test.time")
            self._remove_file("./test_input_output/reference_results.json")

    def tearDown(self):
        pass

if __name__ == '__main__':
    KratosUnittest.main()

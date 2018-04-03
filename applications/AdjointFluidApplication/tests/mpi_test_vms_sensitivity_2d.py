import os
from KratosMultiphysics import *
import KratosMultiphysics.mpi as KratosMPI
import KratosMultiphysics.KratosUnittest as KratosUnittest
import test_MainKratosMPI

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

    def _create_test(self, parameter_file_name):
        with open(parameter_file_name + '_parameters.json', 'r') as parameter_file:
            project_parameters = Parameters(parameter_file.read())
            parameter_file.close()
        test = test_MainKratosMPI.MainKratos(project_parameters)
        return test

    def solve(self, parameter_file_name):
        test = self._create_test(parameter_file_name)
        test.Solve()

    def test_Cylinder(self):
        with ControlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):
            # solve fluid
            self.solve('test_vms_sensitivity_2d/mpi_cylinder_test')
            # solve adjoint
            test = self._create_test('test_vms_sensitivity_2d/mpi_cylinder_test_adjoint')
            test.Solve()
            KratosMPI.mpi.world.barrier()
            rank = test.main_model_part.GetCommunicator().MyPID()
            # remove files
            if rank == 0:
                self._remove_file("./test_vms_sensitivity_2d/mpi_cylinder_test_probe1.dat")
                self._remove_file("./test_vms_sensitivity_2d/mpi_cylinder_test_probe2.dat")
                self._remove_file("./test_vms_sensitivity_2d/mpi_cylinder_test_adjoint_probe1.dat")
                self._remove_file("./test_vms_sensitivity_2d/mpi_cylinder_test_adjoint_probe2.dat")
                self._remove_file("./test_vms_sensitivity_2d/mpi_cylinder_test_adjoint_probe3.dat")
                self._remove_file("./test_vms_sensitivity_2d/cylinder_test.time")
                self._remove_h5_files("MainModelPart")
            self._remove_file("./test_vms_sensitivity_2d/cylinder_test_" + str(rank) + ".time")
            self._remove_file("./test_vms_sensitivity_2d/cylinder_test_" + str(rank) + ".mdpa")

    def tearDown(self):
        pass

if __name__ == '__main__':
    KratosUnittest.main()

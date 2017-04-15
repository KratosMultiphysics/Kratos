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
    
    def removeFile(self, file_path):
        if os.path.isfile(file_path):
            os.remove(file_path)

    def createTest(self, parameter_file_name):
        with open(parameter_file_name + '_parameters.json', 'r') as parameter_file:
            project_parameters = Parameters(parameter_file.read())
            parameter_file.close()
        test = test_MainKratosMPI.MainKratos(project_parameters)
        return test

    def solve(self, parameter_file_name):
        test = self.createTest(parameter_file_name)
        test.Solve()

    def test_Cylinder(self):
        with ControlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):
            # solve fluid
            self.solve('test_vms_sensitivity_2d/mpi_cylinder_test')
            # solve adjoint
            test = self.createTest('test_vms_sensitivity_2d/mpi_cylinder_test_adjoint')
            test.Solve()
            KratosMPI.mpi.world.barrier()
            rank = test.main_model_part.GetCommunicator().MyPID()
            # remove files
            if rank == 0:
                self.removeFile("./test_vms_sensitivity_2d/mpi_cylinder_test_probe1.dat")
                self.removeFile("./test_vms_sensitivity_2d/mpi_cylinder_test_probe2.dat")
                self.removeFile("./test_vms_sensitivity_2d/mpi_cylinder_test_adjoint_probe1.dat")
                self.removeFile("./test_vms_sensitivity_2d/mpi_cylinder_test_adjoint_probe2.dat")
                self.removeFile("./test_vms_sensitivity_2d/mpi_cylinder_test_adjoint_probe3.dat")
                self.removeFile("./test_vms_sensitivity_2d/cylinder_test.time")
            self.removeFile("./test_vms_sensitivity_2d/cylinder_test_" + str(rank) + ".time")
            self.removeFile("./test_vms_sensitivity_2d/cylinder_test_" + str(rank) + ".mdpa")
            self.removeFile("./test_vms_sensitivity_2d/cylinder_test_" + str(rank) + ".h5")

    def tearDown(self):
        pass

if __name__ == '__main__':
    KratosUnittest.main()

import os
from KratosMultiphysics import *
import KratosMultiphysics.mpi as KratosMPI
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

    def removeFile(self, filepath):
        if os.path.isfile(filepath):
            os.remove(filepath)

    def createTest(self, parameter_file_name):
        with open(parameter_file_name + '_parameters.json', 'r') as parameter_file:
            project_parameters = Parameters(parameter_file.read())
            parameter_file.close()
        test = test_MainKratos.MainKratos(project_parameters)
        return test

    def test_Rectangle_2D3N(self):
        with ControlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):
            test = self.createTest('test_structural_mesh_motion_2d/mpi_rectangle_2D3N_test')
            test.Solve()
            KratosMPI.mpi.world.barrier()
            rank = test.main_model_part.GetCommunicator().MyPID()
            # remove files
            if rank == 0:
                self.removeFile("./test_structural_mesh_motion_2d/rectangle_2D3N_test_probe1.dat")
                self.removeFile("./test_structural_mesh_motion_2d/rectangle_2D3N_test_probe2.dat")
                self.removeFile("./test_structural_mesh_motion_2d/rectangle_2D3N_test.time")
            self.removeFile("./test_structural_mesh_motion_2d/rectangle_2D3N_test_" + str(rank) + ".time")
            self.removeFile("./test_structural_mesh_motion_2d/rectangle_2D3N_test_" + str(rank) + ".mdpa")

    def tearDown(self):
        pass

if __name__ == '__main__':
    KratosUnittest.main()

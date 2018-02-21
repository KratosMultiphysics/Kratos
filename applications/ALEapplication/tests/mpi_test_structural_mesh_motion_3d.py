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

    def test_Rectangle_3D8N(self):
        with ControlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):
            test = self.createTest('test_structural_mesh_motion_3d/mpi_rectangle_3D8N_test')
            test.Solve()
            KratosMPI.mpi.world.barrier()
            rank = KratosMPI.mpi.rank
            # remove files
            if rank == 0:
                self.removeFile("./test_mdpa_files/rectangle_3D8N_test.time")

    def test_Rectangle_3D4N(self):
        with ControlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):
            test = self.createTest('test_structural_mesh_motion_3d/mpi_rectangle_3D4N_test')
            test.Solve()
            KratosMPI.mpi.world.barrier()
            rank = KratosMPI.mpi.rank
            # remove files
            if rank == 0:
                self.removeFile("./test_mdpa_files/rectangle_3D4N_test.time")


    def tearDown(self):
        pass

if __name__ == '__main__':
    KratosUnittest.main()

import os
import KratosMultiphysics as Kratos
import KratosMultiphysics.FluidDynamicsApplication
import KratosMultiphysics.mpi as KratosMPI
import KratosMultiphysics.KratosUnittest as KratosUnittest
from fluid_dynamics_analysis import FluidDynamicsAnalysis


missing_applications_message = ["Missing required application(s):",]
have_required_applications = True

try:
    import KratosMultiphysics.HDF5Application as kh5
except ImportError:
    have_required_applications = False
    missing_applications_message.append("HDF5Application")

from fluid_dynamics_analysis import FluidDynamicsAnalysis
from adjoint_fluid_analysis import AdjointFluidAnalysis

class ControlledExecutionScope:
    def __init__(self, scope):
        self.currentPath = os.getcwd()
        self.scope = scope

    def __enter__(self):
        os.chdir(self.scope)

    def __exit__(self, type, value, traceback):
        os.chdir(self.currentPath)

@KratosUnittest.skipUnless(have_required_applications," ".join(missing_applications_message))

class AdjointMPIVMSSensitivity(KratosUnittest.TestCase):

    def setUp(self):
        pass

    def _remove_file(self, file_path):
        if os.path.isfile(file_path):
            os.remove(file_path)

    def _remove_h5_files(self, model_part_name):
        for name in os.listdir():
            if name.find(model_part_name) == 0:
                self._remove_file(name)

    def testCylinder(self):
        with ControlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):
            # solve fluid
            model = Kratos.Model()
            with open('AdjointVMSSensitivity2DTest/mpi_cylinder_test_parameters.json', 'r') as parameter_file:
                project_parameters = Kratos.Parameters(parameter_file.read())
                parameter_file.close()
            primal_simulation = FluidDynamicsAnalysis(model,project_parameters)
            primal_simulation.Run()
            KratosMPI.mpi.world.barrier()

            # solve adjoint
            with open('AdjointVMSSensitivity2DTest/mpi_cylinder_test_adjoint_parameters.json', 'r') as parameter_file:
                project_parameters = Kratos.Parameters(parameter_file.read())
                parameter_file.close()

            adjoint_model = Kratos.Model()
            adjoint_simulation = AdjointFluidAnalysis(adjoint_model,project_parameters)
            adjoint_simulation.Run()
            KratosMPI.mpi.world.barrier()
            rank = adjoint_simulation._solver.main_model_part.GetCommunicator().MyPID()
            # remove files
            if rank == 0:
                self._remove_file("./AdjointVMSSensitivity2DTest/mpi_cylinder_test_probe1.dat")
                self._remove_file("./AdjointVMSSensitivity2DTest/mpi_cylinder_test_probe2.dat")
                self._remove_file("./AdjointVMSSensitivity2DTest/mpi_cylinder_test_adjoint_probe1.dat")
                self._remove_file("./AdjointVMSSensitivity2DTest/mpi_cylinder_test_adjoint_probe2.dat")
                self._remove_file("./AdjointVMSSensitivity2DTest/mpi_cylinder_test_adjoint_probe3.dat")
                self._remove_file("./AdjointVMSSensitivity2DTest/cylinder_test.time")
                self._remove_h5_files("primal")
            self._remove_file("./AdjointVMSSensitivity2DTest/cylinder_test_" + str(rank) + ".time")
            self._remove_file("./AdjointVMSSensitivity2DTest/cylinder_test_" + str(rank) + ".mdpa")

    def tearDown(self):
        pass

if __name__ == '__main__':
    KratosUnittest.main()

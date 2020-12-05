import os
import KratosMultiphysics as Kratos
import KratosMultiphysics.FluidDynamicsApplication
import KratosMultiphysics.mpi as KratosMPI
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utils
from  KratosMultiphysics.kratos_utilities import CheckIfApplicationsAvailable

missing_applications_message = ["Missing required application(s):",]
have_required_applications = CheckIfApplicationsAvailable("HDF5Application")
if have_required_applications:
    import KratosMultiphysics.HDF5Application as kh5
else:
    missing_applications_message.append("HDF5Application")

from KratosMultiphysics.FluidDynamicsApplication.fluid_dynamics_analysis import FluidDynamicsAnalysis
from KratosMultiphysics.FluidDynamicsApplication.adjoint_fluid_analysis import AdjointFluidAnalysis

class ControlledExecutionScope:
    def __init__(self, scope):
        self.currentPath = os.getcwd()
        self.scope = scope

    def __enter__(self):
        os.chdir(self.scope)

    def __exit__(self, type, value, traceback):
        os.chdir(self.currentPath)

@KratosUnittest.skipUnless(have_required_applications," ".join(missing_applications_message))
class AdjointMPIQSVMSSensitivity(KratosUnittest.TestCase):

    def setUp(self):
        pass

    def _removeH5Files(self, model_part_name):
        for name in os.listdir():
            if name.find(model_part_name) == 0:
                kratos_utils.DeleteFileIfExisting(name)

    def testCylinder(self):
        with ControlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):
            # solve fluid
            model = Kratos.Model()
            with open('AdjointQSVMSSensitivity2DTest/cylinder_test_parameters.json', 'r') as parameter_file:
                project_parameters = Kratos.Parameters(parameter_file.read())
            project_parameters["problem_data"]["parallel_type"].SetString("MPI")

            primal_simulation = FluidDynamicsAnalysis(model,project_parameters)
            primal_simulation.Run()
            Kratos.DataCommunicator.GetDefault().Barrier()

            # solve adjoint
            with open('AdjointQSVMSSensitivity2DTest/cylinder_test_adjoint_parameters.json', 'r') as parameter_file:
                project_parameters = Kratos.Parameters(parameter_file.read())
            project_parameters["problem_data"]["parallel_type"].SetString("MPI")

            adjoint_model = Kratos.Model()
            adjoint_simulation = AdjointFluidAnalysis(adjoint_model,project_parameters)
            adjoint_simulation.Run()
            Kratos.DataCommunicator.GetDefault().Barrier()

            self._removeH5Files("MainModelPart")
            kratos_utils.DeleteFileIfExisting("./tests.post.lst")
            kratos_utils.DeleteDirectoryIfExisting("AdjointVMSSensitivity2DTest/cylinder_test_partitioned")

    def tearDown(self):
        pass

if __name__ == '__main__':
    KratosUnittest.main()

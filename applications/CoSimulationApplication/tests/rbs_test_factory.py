import KratosMultiphysics as KM
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utils
from KratosMultiphysics.CoSimulationApplication.co_simulation_analysis import CoSimulationAnalysis
from KratosMultiphysics.CoSimulationApplication.solver_wrappers.rigid_body.rigid_body_solver import RigidBodySolver


class RBSTestFactory(KratosUnittest.TestCase):
    def setUp(self):
        # Within this location context:
        with KratosUnittest.WorkFolderScope(".", __file__):

            # Reading the ProjectParameters
            with open(self.file_name + "_parameters.json",'r') as parameter_file:
                ProjectParameters = KM.Parameters(parameter_file.read())

            # To avoid many prints
            if ProjectParameters["problem_data"]["echo_level"].GetInt() == 0:
                KM.Logger.GetDefaultOutput().SetSeverity(KM.Logger.Severity.WARNING)
            else:
                KM.Logger.GetDefaultOutput().SetSeverity(KM.Logger.Severity.INFO)

            # Creating the test
            if self.cosim_run:
                self.test = CoSimulationAnalysis(ProjectParameters)
            else:
                model = KM.Model()
                self.test = RigidBodySolver(model, ProjectParameters)

            # self.addCleanup(kratos_utils.DeleteTimeFiles, self.file_name)

    def test_execution(self):
        # Within this location context:
        with KratosUnittest.WorkFolderScope(".", __file__):
            self.test.Run()


class TestRBSStandalone(RBSTestFactory):
    file_name = "rbs_test/RBS_standalone/RBS_standalone"
    cosim_run = False

class TestRBSRBS(RBSTestFactory):
    file_name = "rbs_test/RBS_RBS/RBS_RBS"
    cosim_run = True

@KratosUnittest.skipIfApplicationsNotAvailable("FluidDynamicsApplication", "MappingApplication", "LinearSolversApplication")
class TestBarc2DRigidBody(RBSTestFactory):
    file_name = "rbs_test/Barc2DRigidBody/Barc2DRigidBody"
    cosim_run = True


if __name__ == '__main__':
    KratosUnittest.main()

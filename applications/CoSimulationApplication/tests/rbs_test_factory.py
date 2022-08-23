import KratosMultiphysics as KM
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utils
from KratosMultiphysics.CoSimulationApplication.co_simulation_analysis import CoSimulationAnalysis
from KratosMultiphysics.CoSimulationApplication.solver_wrappers.rigid_body.rigid_body_solver import RigidBodySolver


data_comm = KM.Testing.GetDefaultDataCommunicator()

import os

try:
    import numpy
    numpy_available = True
except ImportError:
    numpy_available = False

have_cfd_dependencies = kratos_utils.CheckIfApplicationsAvailable("FluidDynamicsApplication", "MappingApplication", "LinearSolversApplication")

# def GetFilePath(fileName):
#     return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)


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
            model = KM.Model()
            self.test = RigidBodySolver(model, ProjectParameters)

            # self.addCleanup(kratos_utils.DeleteTimeFiles, self.file_name)

    def test_execution(self):
        # Within this location context:
        with KratosUnittest.WorkFolderScope(".", __file__):
            self.test.Run()


class CoSimTestFactory(KratosUnittest.TestCase):
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
            model = KM.Model()
            self.test = CoSimulationAnalysis(ProjectParameters)

            # self.addCleanup(kratos_utils.DeleteTimeFiles, self.file_name)

    def test_execution(self):
        # Within this location context:
        with KratosUnittest.WorkFolderScope(".", __file__):
            self.test.Run()

    # def tearDown(self):
    #     # Within this location context:
    #     with KratosUnittest.WorkFolderScope(".", __file__):
    #         self.test.Finalize()



class RBSStandaloneTest(RBSTestFactory):
    if not numpy_available:
        CoSimRBSTestFactory.skipTest("Numpy not available")
    file_name = "rbs_test/RBS_standalone/RBS_standalone"

class RBSRBSTest(CoSimTestFactory):
    if not numpy_available:
        CoSimRBSTestFactory.skipTest("Numpy not available")
    file_name = "rbs_test/RBS_RBS/RBS_RBS"

@KratosUnittest.skipIfApplicationsNotAvailable("FluidDynamicsApplication")
class Barc2DRigidBodyTest(CoSimTestFactory):
    # if not numpy_available:
    #     self.skipTest("Numpy not available")
    # if not False:
    #     self.skipTest("CFD dependencies are not available!")
    file_name = "rbs_test/Barc2DRigidBody/Barc2DRigidBody"


if __name__ == '__main__':
    KratosUnittest.main()

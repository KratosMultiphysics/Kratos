from KratosMultiphysics import *
from KratosMultiphysics.MultilevelMonteCarloApplication import *

try:
    import KratosMultiphysics.ConvectionDiffusionApplication
    have_convection_diffusion_application = True
except ImportError:
    have_convection_diffusion_application = False

import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utilities

import os

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(__file__), fileName)

@KratosUnittest.skipUnless(have_convection_diffusion_application,"Missing required application: ConvectionDiffusionApplication")
class KratosMultilevelMonteCarloGeneralTestsAuxiliary(KratosUnittest.TestCase):

    def MonteCarloTest(self):
        with KratosUnittest.WorkFolderScope(os.path.join(self.folder_name),__file__,add_to_path=True):
            import KratosMultiphysics.MultilevelMonteCarloApplication.mc_utilities as mc_utilities
            from simulation_definition import SimulationScenario

            # set the ProjectParameters.json path
            project_parameters_path = "problem_settings/parameters_poisson_square_2d_coarse.json"
            # set parameters of the MC simulation"""
            parameters_x_monte_carlo_path = "problem_settings/parameters_x_monte_carlo.json"
            # contruct MonteCarlo or MultilevelMonteCarlo class
            mc_manager = mc_utilities.MonteCarlo(parameters_x_monte_carlo_path,project_parameters_path,SimulationScenario)
            # execute algorithm
            mc_manager.Run()
            """delete files"""
            kratos_utilities.DeleteFileIfExisting(os.path.join(self.folder_name,"poisson_square_2d.post.bin"))
            kratos_utilities.DeleteFileIfExisting(os.path.join(self.folder_name,"poisson_square_2d.post.lst"))

    def MultilevelMonteCarloTest(self):
        with KratosUnittest.WorkFolderScope(os.path.join(self.folder_name),__file__,add_to_path=True):
            import KratosMultiphysics.MultilevelMonteCarloApplication.mlmc_utilities as mlmc_utilities
            from simulation_definition import SimulationScenario

            # set the ProjectParameters.json path
            project_parameters_path = "problem_settings/parameters_poisson_square_2d_coarse.json"
            # set parameters of the MLMC simulation
            parameters_x_monte_carlo_path = "problem_settings/parameters_x_monte_carlo.json"
            # set parameters of the metric of the adaptive refinement utility and set parameters of the remesh of the adaptive refinement utility
            parameters_refinement_path = "problem_settings/parameters_refinement.json"

            # contruct MultilevelMonteCarlo class
            mlmc_manager = mlmc_utilities.MultilevelMonteCarlo(parameters_x_monte_carlo_path,project_parameters_path,parameters_refinement_path,SimulationScenario)
            mlmc_manager.Run()
            """delete files"""
            kratos_utilities.DeleteFileIfExisting(os.path.join(self.folder_name,"poisson_square_2d.post.bin"))
            kratos_utilities.DeleteFileIfExisting(os.path.join(self.folder_name,"poisson_square_2d.post.lst"))


    def _runTest(self):
        # Code here that runs the test.
        pass

class KratosMultilevelMonteCarloGeneralTests(KratosMultilevelMonteCarloGeneralTestsAuxiliary):

    def setUp(self):
        self.folder_name_case_MC = "poisson_square_2d"
        self.folder_name_case_MLMC = "poisson_square_2d"
        pass

    def testMonteCarlo(self):
        self.folder_name = self.folder_name_case_MC
        self.MonteCarloTest()

    def testMultilevelMonteCarlo(self):
        self.folder_name = self.folder_name_case_MLMC
        self.MultilevelMonteCarloTest()

    def tearDown(self):
        # Code here will be placed AFTER every test in this TestCase.
        pass

if __name__ == '__main__':
    KratosUnittest.main()


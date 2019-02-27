from KratosMultiphysics import *
from KratosMultiphysics.MultilevelMonteCarloApplication import *
try:
    import KratosMultiphysics.ConvectionDiffusionApplication
    have_convection_diffusion_application = True
except ImportError:
    have_convection_diffusion_application = False

# import test_mc_utilities as mc_utilities
# import test_cmlmc_utilities as cmlmc_utilities
# import adaptive_refinement_utilities as refinement

import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utilities

import os

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(__file__), fileName)

@KratosUnittest.skipUnless(have_convection_diffusion_application,"Missing required application: ConvectionDiffusionApplication")
class KratosMultilevelMonteCarloGeneralTests(KratosUnittest.TestCase):

    def setUp(self):
        # Code here will be placed BEFORE every test in this TestCase.
        pass

    def tearDown(self):
        # Code here will be placed AFTER every test in this TestCase.
        pass

    def testSmallExample(self):
        self.assertEqual(True, True)

    def testNightlyFirstExample(self):
        self.assertEqual(True, True)

    def testNightlySecondExample(self):
        self.assertEqual(True, True)

    def testMonteCarloAnalysis(self):
        print(dir(KratosUnittest))
        self.test_case_folder_name = "poisson_square_2d"
        with KratosUnittest.WorkFolderScope(os.path.join("..","test_examples", self.test_case_folder_name),__file__):
            import test_mc_utilities as mc_utilities
            from simulation_definition import SimulationScenario

            # set the ProjectParameters.json path
            project_parameters_path = "problem_settings/parameters_poisson_square_2d_coarse.json"
            # customize setting parameters of the MC simulation"""
            settings_MC_simulation = KratosMultiphysics.Parameters("""
            {
                "tolerance" : 0.1,
                "cphi" : 5e-1,
                "batch_size" : 20,
                "convergence_criteria" : "MC_higher_moments_sequential_stopping_rule"
            }
            """)
            # contruct MonteCarlo or MultilevelMonteCarlo class
            mc_manager = mc_utilities.MonteCarlo(settings_MC_simulation,project_parameters_path,SimulationScenario)
            # execute algorithm
            mc_manager.Run()
            '''delete .time, .bin files'''
            kratos_utilities.DeleteFileIfExisting("PoissonSquareTest/square_coarse_2d.time")
            kratos_utilities.DeleteFileIfExisting("tests.post.lst")
            kratos_utilities.DeleteFileIfExisting("MLMCLaplacian.post.bin")

    # def testMultilevelMonteCarloAnalysis(self):
    #     # set the ProjectParameters.json path
    #     project_parameters_path = "problem_settings/parameters_poisson_square_2d_coarse.json"
    #     # customize setting parameters of the MLMC simulation
    #     settings_MLMC_simulation = KratosMultiphysics.Parameters("""
    #     {
    #         "tol0"                            : 0.25,
    #         "tolF"                            : 0.1,
    #         "cphi"                            : 1.0,
    #         "number_samples_screening"        : 25,
    #         "Lscreening"                      : 2,
    #         "Lmax"                            : 4,
    #         "initial_mesh_size"               : 0.5
    #     }
    #     """)
    #     # customize setting parameters of the metric of the adaptive refinement utility
    #     custom_metric_refinement_parameters = KratosMultiphysics.Parameters("""
    #         {
    #             "hessian_strategy_parameters"           :{
    #                     "metric_variable"               : ["TEMPERATURE"],
    #                     "estimate_interpolation_error"  : false,
    #                     "interpolation_error"           : 0.004
    #             },
    #             "anisotropy_remeshing"                  : true,
    #             "anisotropy_parameters":{
    #                 "reference_variable_name"           : "TEMPERATURE",
    #                 "hmin_over_hmax_anisotropic_ratio"  : 0.15,
    #                 "boundary_layer_max_distance"       : 1.0,
    #                 "interpolation"                     : "Linear"
    #             },
    #             "local_gradient_variable"               : "TEMPERATURE"
    #         }
    #     """)
    #     # customize setting parameters of the remesh of the adaptive refinement utility
    #     custom_remesh_refinement_settings = KratosMultiphysics.Parameters("""
    #         {
    #             "echo_level"                            : 0
    #         }
    #     """)
    #     # contruct MultilevelMonteCarlo class
    #     mlmc_manager = cmlmc_utilities.MultilevelMonteCarlo(settings_MLMC_simulation,project_parameters_path,custom_metric_refinement_parameters,custom_remesh_refinement_settings,SimulationScenario)
    #     mlmc_manager.Run()
    #     '''delete .time, .bin files'''
    #     kratos_utilities.DeleteFileIfExisting("PoissonSquareTest/square_coarse_2d.time")
    #     kratos_utilities.DeleteFileIfExisting("tests.post.lst")
    #     kratos_utilities.DeleteFileIfExisting("MLMCLaplacian.post.bin")

    def _runTest(self):
        # Code here that runs the test.
        pass

if __name__ == '__main__':
    KratosUnittest.main()


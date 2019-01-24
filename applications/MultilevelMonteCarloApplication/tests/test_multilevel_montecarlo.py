from KratosMultiphysics import *
from KratosMultiphysics.MultilevelMonteCarloApplication import *
try:
    import KratosMultiphysics.ConvectionDiffusionApplication
    have_convection_diffusion_application = True
except ImportError:
    have_convection_diffusion_application = False

import test_montecarlo_analysis as MC
import test_multilevel_montecarlo_analysis as MLMC
import test_cmlmc_utilities as mlmc
import test_mc_utilities as mc
import adaptive_refinement_utilities as refinement

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
        '''set the ProjectParameters.json path'''
        parameter_file_name = "PoissonSquareTest/parameters_poisson_coarse.json"
        '''create a serialization of the model and of the project parameters'''
        pickled_model,pickled_parameters = MC.SerializeModelParameters_Task(parameter_file_name)
        '''evaluate the exact expected value of Q (sample = 1.0)'''
        ExactExpectedValueQoI = MC.ExecuteExactMonteCarloAnalysis_Task(pickled_model,pickled_parameters)
        '''customize setting parameters of the ML simulation'''
        settings_MC_simulation = KratosMultiphysics.Parameters("""
        {
            "tolerance" : 0.1,
            "cphi" : 1e-1,
            "batch_size" : 20,
            "convergence_criteria" : "MC_higher_moments_sequential_stopping_rule"
        }
        """)
        '''contruct MonteCarlo class'''
        mc_class = mc.MonteCarlo(settings_MC_simulation)
        '''start MC algorithm'''
        while mc_class.convergence is not True:
            mc_class.InitializeMCPhase()
            for instance in range (mc_class.difference_number_samples[0]):
                mc_class.AddResults(MC.ExecuteMonteCarloAnalysis_Task(pickled_model,pickled_parameters))
            mc_class.FinalizeMCPhase()
        print("\nMC mean = ",mc_class.QoI.mean,"exact mean = ",ExactExpectedValueQoI)
        relative_error = (mc_class.QoI.mean[0]-ExactExpectedValueQoI)/ExactExpectedValueQoI
        print("relative error = ",relative_error)
        '''delete .time, .bin files'''
        kratos_utilities.DeleteFileIfExisting("PoissonSquareTest/square_coarse_2d.time")
        kratos_utilities.DeleteFileIfExisting("tests.post.lst")
        kratos_utilities.DeleteFileIfExisting("MLMCLaplacian.post.bin")


    def testMultilevelMonteCarloAnalysis(self):
        '''set the ProjectParameters.json path'''
        parameter_file_name = "PoissonSquareTest/parameters_poisson_coarse.json"
        '''create a serialization of the model and of the project parameters'''
        pickled_model,pickled_parameters = MLMC.SerializeModelParameters_Task(parameter_file_name)
        '''customize setting parameters of the ML simulation'''
        settings_ML_simulation = KratosMultiphysics.Parameters("""
        {
            "tol0"                            : 0.25,
            "tolF"                            : 0.1,
            "cphi"                            : 1.0,
            "number_samples_screening"        : 15,
            "Lscreening"                      : 2,
            "Lmax"                            : 4,
            "initial_mesh_size"               : 0.5
        }
        """)
        '''customize setting parameters of the metric of the adaptive refinement utility'''
        settings_metric_refinement = KratosMultiphysics.Parameters("""
            {
                "hessian_strategy_parameters"              :{
                        "metric_variable"                  : ["TEMPERATURE"],
                        "estimate_interpolation_error"     : false,
                        "interpolation_error"              : 0.004
                },
                "anisotropy_remeshing"              : true,
                "anisotropy_parameters":{
                    "reference_variable_name"          : "TEMPERATURE",
                    "hmin_over_hmax_anisotropic_ratio" : 0.15,
                    "boundary_layer_max_distance"      : 1.0,
                    "interpolation"                    : "Linear"
                },
                "local_gradient_variable"           : "TEMPERATURE"
            }
        """)
        '''customize setting parameters of the remesh of the adaptive refinement utility'''
        settings_remesh_refinement = KratosMultiphysics.Parameters("""
            {
                "echo_level"                       : 0
            }
        """)
        pickled_settings_metric_refinement,pickled_settings_remesh_refinement = MLMC.SerializeRefinementParameters(settings_metric_refinement,settings_remesh_refinement)

        '''contruct MultilevelMonteCarlo class'''
        mlmc_class = mlmc.MultilevelMonteCarlo(settings_ML_simulation)
        ''''start screening phase'''
        for lev in range(mlmc_class.current_number_levels+1):
            for instance in range (mlmc_class.number_samples[lev]):
                mlmc_class.AddResults(MLMC.ExecuteMultilevelMonteCarloAnalisys(lev,pickled_model,pickled_parameters,mlmc_class.sizes_mesh,pickled_settings_metric_refinement,pickled_settings_remesh_refinement))
        '''finalize screening phase'''
        mlmc_class.FinalizeScreeningPhase()
        '''start MLMC phase'''
        while mlmc_class.convergence is not True:
            '''initialize MLMC phase'''
            mlmc_class.InitializeMLMCPhase()
            '''MLMC execution phase'''
            for lev in range (mlmc_class.current_number_levels+1):
                for instance in range (mlmc_class.difference_number_samples[lev]):
                    mlmc_class.AddResults(MLMC.ExecuteMultilevelMonteCarloAnalisys(lev,pickled_model,pickled_parameters,mlmc_class.sizes_mesh,pickled_settings_metric_refinement,pickled_settings_remesh_refinement))
            '''finalize MLMC phase'''
            mlmc_class.FinalizeMLMCPhase()
        print("\niterations = ",mlmc_class.current_iteration,\
        "total error TErr computed = ",mlmc_class.TErr,"mean MLMC QoI = ",mlmc_class.mean_mlmc_QoI)
        '''delete .time, .bin files'''
        kratos_utilities.DeleteFileIfExisting("PoissonSquareTest/square_coarse_2d.time")
        kratos_utilities.DeleteFileIfExisting("tests.post.lst")
        kratos_utilities.DeleteFileIfExisting("MLMCLaplacian.post.bin")



    def _runTest(self):
        # Code here that runs the test.
        pass

if __name__ == '__main__':
    KratosUnittest.main()


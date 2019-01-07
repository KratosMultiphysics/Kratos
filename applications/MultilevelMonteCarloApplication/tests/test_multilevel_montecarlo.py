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
import adaptive_refinement_utilities as refinement

import KratosMultiphysics.KratosUnittest as KratosUnittest

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
        parameter_file_name = "Poisson2dTest/poisson_2d_project_parameters.json"
        '''create a serialization of the model and of the project parameters'''
        pickled_model,pickled_parameters = MC.SerializeModelParameters_Task(parameter_file_name)
        '''evaluate the exact expected value of Q (sample = 1.0)'''
        ExactExpectedValueQoI = MC.ExecuteExactMonteCarlo_Task(pickled_model,pickled_parameters) 
        '''define number samples'''
        number_samples = 10
        '''initialize Quantity of Interest'''
        QoI = mlmc.StatisticalVariable(0) # number of levels = 0 (we only have one level), needed using this class
        '''to exploit StatisticalVariable UpdateOnePassMeanVariance function we need to initialize a level 0 in values, mean, sample variance and second moment
        and store in this level the informations'''
        QoI.values = [[] for i in range (1)]
        QoI.mean = [[] for i in range (1)]
        QoI.second_moment = [[] for i in range (1)]
        QoI.sample_variance = [[] for i in range (1)]
        for instance in range (0,number_samples):
            QoI.values[0].append(MC.ExecuteMonteCarlo_Task(pickled_model,pickled_parameters))
        '''Compute mean, second moment and sample variance'''
        for i_sample in range (0,number_samples):
            QoI.UpdateOnepassMeanVariance(0,i_sample)
        '''Evaluation of the relative error between the computed mean value and the expected value of the QoI'''
        relative_error = MC.CompareMean_Task(QoI.mean[0],ExactExpectedValueQoI)
        print("relative error: ",relative_error)
        '''delete .time file'''
        kratos_utilities.DeleteFileIfExisting("SquareCoarse/square_coarse_2d.time")
        kratos_utilities.DeleteFileIfExisting("tests.post.lst")
        kratos_utilities.DeleteFileIfExisting("MLMCLaplacian.post.bin")

    def testMultilevelMonteCarloAnalysis(self):
        '''set the ProjectParameters.json path'''
        parameter_file_name = "Poisson2dTest/poisson_2d_project_parameters.json"
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
        '''contruct MultilevelMonteCarlo class'''
        mlmc_class = mlmc.MultilevelMonteCarlo(settings_ML_simulation)
        ''''start screening phase'''
        for lev in range(mlmc_class.current_number_levels+1):
            for instance in range (mlmc_class.number_samples[lev]):
                mlmc_class.AddResults(MLMC.ExecuteMultilevelMonteCarloAnalisys(lev,pickled_model,pickled_parameters,mlmc_class.sizes_mesh))
        '''finalize screening phase'''
        mlmc_class.FinalizeScreeningPhase()
        '''start MLMC phase'''
        while mlmc_class.convergence is not True:
            '''initialize MLMC phase'''
            mlmc_class.InitializeMLMCPhase()
            '''MLMC execution phase'''
            for lev in range (mlmc_class.current_number_levels+1):
                for instance in range (mlmc_class.difference_number_samples[lev]):
                    mlmc_class.AddResults(MLMC.ExecuteMultilevelMonteCarloAnalisys(lev,pickled_model,pickled_parameters,mlmc_class.sizes_mesh))
            '''finalize MLMC phase'''
            mlmc_class.FinalizeMLMCPhase()
        print("\niterations = ",mlmc_class.current_iteration,\
        "total error TErr computed = ",mlmc_class.TErr,"mean MLMC QoI = ",mlmc_class.mean_mlmc_QoI)
        kratos_utilities.DeleteFileIfExisting("SquareCoarse/square_coarse_2d.time")
        kratos_utilities.DeleteFileIfExisting("tests.post.lst")
        kratos_utilities.DeleteFileIfExisting("MLMCLaplacian.post.bin")



    def _runTest(self):
        # Code here that runs the test.
        pass

if __name__ == '__main__':
    KratosUnittest.main()


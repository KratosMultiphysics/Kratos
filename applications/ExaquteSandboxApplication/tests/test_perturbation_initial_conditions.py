import KratosMultiphysics
import KratosMultiphysics.ConvectionDiffusionApplication
from  KratosMultiphysics.ConvectionDiffusionApplication.convection_diffusion_analysis import ConvectionDiffusionAnalysis
import KratosMultiphysics.KratosUnittest as UnitTest

from math import sqrt
try:
    import scipy.fftpack
    is_fft_loaded = True
except:
    is_fft_loaded = False

def ComputeDivergence(analysis):
    # compute velocity gradients
    KratosMultiphysics.ComputeNodalGradientProcess(analysis._GetSolver().main_model_part,KratosMultiphysics.VELOCITY_X,KratosMultiphysics.VELOCITY_X_GRADIENT).Execute()
    KratosMultiphysics.ComputeNodalGradientProcess(analysis._GetSolver().main_model_part,KratosMultiphysics.VELOCITY_Y,KratosMultiphysics.VELOCITY_Y_GRADIENT).Execute()
    KratosMultiphysics.ComputeNodalGradientProcess(analysis._GetSolver().main_model_part,KratosMultiphysics.VELOCITY_Z,KratosMultiphysics.VELOCITY_Z_GRADIENT).Execute()
    # compute l2 divergence norm
    main_model_part_name = analysis.project_parameters["problem_data"]["model_part_name"].GetString()
    squared_divergence = 0
    for node in analysis.model.GetModelPart(main_model_part_name).Nodes:
        divergence = node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_X_GRADIENT)[0] + node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_Y_GRADIENT)[1] + node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_Z_GRADIENT)[2]
        squared_divergence += divergence**2
    return sqrt(squared_divergence)

@UnitTest.skipIfApplicationsNotAvailable("ConvectionDiffusionApplication","MappingApplication")
class PerturbedDivergenceFreeInitialConditionProcessTest(UnitTest.TestCase):

    def testPerturbedDivergenceFreeInitialConditionProcess(self):
        if not is_fft_loaded:
            self.skipTest("Missing required Python library: scipy.fftpack")

        parameter_file_name = "test_perturbation_initial_conditions/parameters_poisson_rectangle_2d.json"
        with open(parameter_file_name,'r') as parameter_file:
            parameters = KratosMultiphysics.Parameters(parameter_file.read())
        model = KratosMultiphysics.Model()
        simulation = ConvectionDiffusionAnalysis(model, parameters)
        simulation._GetSolver().main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY_X_GRADIENT)
        simulation._GetSolver().main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY_Y_GRADIENT)
        simulation._GetSolver().main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY_Z_GRADIENT)
        simulation.Initialize()
        l2norm_divergence = ComputeDivergence(simulation)
        self.assertAlmostEqual(l2norm_divergence,0.304470818597,delta=1e-12)

if __name__ == '__main__':
    UnitTest.main()
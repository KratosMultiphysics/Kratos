import KratosMultiphysics

import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.ShallowWaterApplication.shallow_water_analysis import ShallowWaterAnalysis

try:
    import scipy
    scipy_available = True
except ImportError:
    scipy_available = False

class ShallowWaterTestFactory(KratosUnittest.TestCase):
    need_numpy = False
    need_scipy = False
    def test_execution(self):
        if self.need_scipy and not scipy_available:
            self.skipTest("scipy not available")
        with KratosUnittest.WorkFolderScope(self.execution_directory, __file__):
            with open(self.execution_file + "_parameters.json",'r') as parameter_file:
                ProjectParameters = KratosMultiphysics.Parameters(parameter_file.read())
            model = KratosMultiphysics.Model()
            test = ShallowWaterAnalysis(model, ProjectParameters)
            test.Run()

class TestConservativeResidualViscosity2D3NElement(ShallowWaterTestFactory):
    execution_directory = "elements_tests"
    execution_file = "conservative_residual_viscosity_2d_3n"

class TestConservativeGradientJump2D3NElement(ShallowWaterTestFactory):
    execution_directory = "elements_tests"
    execution_file = "conservative_gradient_jump_2d_3n"

class TestConservativeFluxCorrected2D3NElement(ShallowWaterTestFactory):
    execution_directory = "elements_tests"
    execution_file = "conservative_flux_corrected_2d_3n"

class TestPrimitive2D3NElement(ShallowWaterTestFactory):
    execution_directory = "elements_tests"
    execution_file = "primitive_2d_3n"

class TestBoussinesq2D3NElement(ShallowWaterTestFactory):
    execution_directory = "elements_tests"
    execution_file = "boussinesq_2d_3n"

class TestSetTopographyProcess(ShallowWaterTestFactory):
    execution_directory = "processes_tests"
    execution_file = "set_topography_process"

class TestVisualizationMeshProcess(ShallowWaterTestFactory):
    execution_directory = "processes_tests"
    execution_file = "visualization_mesh_process"

class TestMacDonaldShockBenchmark(ShallowWaterTestFactory):
    execution_directory = "processes_tests"
    execution_file = "mac_donald_shock_benchmark"
    need_scipy = True
    need_numpy = True

class TestMacDonaldTransitionBenchmark(ShallowWaterTestFactory):
    execution_directory = "processes_tests"
    execution_file = "mac_donald_transition_benchmark"
    need_scipy = True
    need_numpy = True

class TestDamBreakBenchmark(ShallowWaterTestFactory):
    execution_directory = "processes_tests"
    execution_file = "dam_break_benchmark"
    need_scipy = True

class TestDryDamBreakBenchmark(ShallowWaterTestFactory):
    execution_directory = "processes_tests"
    execution_file = "dry_dam_break_benchmark"
    need_scipy = True

class TestPlanarSurfaceInParabolaBenchmark(ShallowWaterTestFactory):
    execution_directory = "processes_tests"
    execution_file = "planar_surface_in_parabola_benchmark"

class TestSolitaryWaveBenchmark(ShallowWaterTestFactory):
    execution_directory = "processes_tests"
    execution_file = "solitary_wave_benchmark"
    need_scipy = True
    need_numpy = True

class TestMeshMovingStrategy(ShallowWaterTestFactory):
    execution_directory = "nightly_tests"
    execution_file = "mesh_moving_strategy"

class TestDamBreakValidation(ShallowWaterTestFactory):
    execution_directory = "validation_tests"
    execution_file = "dam_break_validation"
    need_scipy = True

class TestMacDonaldShockValidation(ShallowWaterTestFactory):
    execution_directory = "validation_tests"
    execution_file = "mac_donald_shock_validation"
    need_scipy = True
    need_numpy = True

class TestSolitaryWaveValidation(ShallowWaterTestFactory):
    execution_directory = "validation_tests"
    execution_file = "solitary_wave_validation"
    need_scipy = True
    need_numpy = True

import KratosMultiphysics

import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.ShallowWaterApplication.shallow_water_analysis import ShallowWaterAnalysis

try:
    import numpy
    numpy_available = True
except ImportError:
    numpy_available = False

try:
    import scipy
    scipy_available = True
except ImportError:
    scipy_available = False

try:
    import mpmath
    mpmath_available = True
except ImportError:
    mpmath_available = False

class ShallowWaterTestFactory(KratosUnittest.TestCase):
    need_numpy = False
    need_scipy = False
    need_mpmath = False
    def test_execution(self):
        if self.need_numpy and not numpy_available:
            self.skipTest("numpy not available")
        if self.need_scipy and not scipy_available:
            self.skipTest("scipy not available")
        if self.need_mpmath and not mpmath_available:
            self.skipTest("mpmath not available")
        with KratosUnittest.WorkFolderScope(self.execution_directory, __file__):
            with open(self.execution_file + "_parameters.json",'r') as parameter_file:
                ProjectParameters = KratosMultiphysics.Parameters(parameter_file.read())
            model = KratosMultiphysics.Model()
            test = ShallowWaterAnalysis(model, ProjectParameters)
            test.Run()

class TestSemiLagrangianShallowWaterElement(ShallowWaterTestFactory):
    execution_directory = "elements_tests"
    execution_file = "semi_lagrangian_swe"

class TestShallowWaterElement(ShallowWaterTestFactory):
    execution_directory = "elements_tests"
    execution_file = "swe"

class TestShallowWater2D3NElement(ShallowWaterTestFactory):
    execution_directory = "elements_tests"
    execution_file = "shallow_water_2d_3n"

class TestMonotonicShallowWater2D3NElement(ShallowWaterTestFactory):
    execution_directory = "elements_tests"
    execution_file = "monotonic_shallow_water_2d_3n"

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

class TestMacDonaldTransitionBenchmark(ShallowWaterTestFactory):
    execution_directory = "processes_tests"
    execution_file = "mac_donald_transition_benchmark"
    need_scipy = True

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
    need_scipy = True

class TestSolitaryWaveBenchmark(ShallowWaterTestFactory):
    execution_directory = "processes_tests"
    execution_file = "solitary_wave_benchmark"
    need_scipy = True
    need_mpmath = True

class TestMeshMovingStrategy(ShallowWaterTestFactory):
    execution_directory = "nightly_tests"
    execution_file = "mesh_moving_strategy"
    need_numpy = True

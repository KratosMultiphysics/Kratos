
# Importing the Kratos Library
import KratosMultiphysics

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.MPMApplication.particle_mechanics_analysis import ParticleMechanicsAnalysis

class MPMTestFactory(KratosUnittest.TestCase):
    def setUp(self):
        # Within this location context:
        with KratosUnittest.WorkFolderScope(".", __file__):

            # Reading the ProjectParameters
            with open(self.file_name + "_parameters.json",'r') as parameter_file:
                ProjectParameters = KratosMultiphysics.Parameters(parameter_file.read())

            self.modify_parameters(ProjectParameters)

            # To avoid many prints
            if (ProjectParameters["problem_data"]["echo_level"].GetInt() == 0):
                KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)

            # Creating the test
            model = KratosMultiphysics.Model()
            self.test = ParticleMechanicsAnalysis(model, ProjectParameters)
            self.test.Initialize()

    def modify_parameters(self, project_parameters):
        """This function can be used in derived classes to modify existing parameters
        before the execution of the test (e.g. switch to MPI)
        """
        pass

    def test_execution(self):
        # Within this location context:
        with KratosUnittest.WorkFolderScope(".", __file__):
            self.test.RunSolutionLoop()

    def tearDown(self):
        # Within this location context:
        with KratosUnittest.WorkFolderScope(".", __file__):
            self.test.Finalize()

### Axis-Symmetric Tests
class AxisSymmetricCircularPlate2DTriTest(MPMTestFactory):
    file_name = "axisym_tests/circular_plate_axisym_test"

### Beam Tests
class BeamCantileverStaticLinearElasticPointLoad2DTriTest(MPMTestFactory):
    file_name = "beam_tests/cantilever_beam/static_point_load_2D_tri_test"

class BeamCantileverStaticLinearElasticLineLoad2DQuadTest(MPMTestFactory):
    file_name = "beam_tests/cantilever_beam/static_line_load_2D_quad_test"

class BeamCantileverStaticLinearElasticParticlePointLoad2DTriTest(MPMTestFactory):
    file_name = "beam_tests/cantilever_beam/particle_point_load_2D_tri_test"

class BeamCantileverStaticLinearElasticSurfaceLoad3DHexaTest(MPMTestFactory):
    file_name = "beam_tests/cantilever_beam/static_surface_load_3D_hexa_test"

class BeamCantileverStaticHyperelasticSelfWeightLoad2DQuadTest(MPMTestFactory):
    file_name = "beam_tests/hyperelastic_cantilever_beam/self_weight_load_2D_quad_test"

class BeamCantileverLinearStaticHyperelasticSelfWeightLoad2DQuadTest(MPMTestFactory):
    file_name = "beam_tests/hyperelastic_cantilever_beam/linear_self_weight_load_2D_quad_test"

class BeamCantileverDynamicConsistentMassTest(MPMTestFactory):
    file_name = "beam_tests/dynamic_cantilever/dynamic_cantilever_consistent_mass_test"

class BeamCantileverDynamicHyperelasticUPTest(MPMTestFactory):
    file_name = "beam_tests/dynamic_UP_hyperelastic_cantilever_test/dynamic_UP_hyperelastic_test"

### Cook's Membrane Tests
class CooksMembraneCompressibleTest(MPMTestFactory):
    file_name = "cooks_membrane_tests/compressible_cook_membrane_2D_test"

class CooksMembraneUPCompressibleTest(MPMTestFactory):
    file_name = "cooks_membrane_tests/UP_compressible_cook_membrane_2D_test"

class CooksMembraneUPIncompressibleTest(MPMTestFactory):
    file_name = "cooks_membrane_tests/UP_incompressible_cook_membrane_2D_test"

### Constitutive Law Tests
class CLLinearElastic3DQuadTest(MPMTestFactory):
    file_name = "cl_tests/solid_cl/linear_elastic_3D_hexa_test"

class CLDispNewtonianFluidTest(MPMTestFactory):
    file_name = "cl_tests/fluid_cl/newtonian_fluid_test"

### Gravity Application Tests
class GravityApplicationTest(MPMTestFactory):
    file_name = "gravity_tests/dynamic_gravity_application_test"

### Gravity Time Step Table Tests
class GravityTimeStepTableTest(MPMTestFactory):
    file_name = "gravity_tests/dynamic_gravity_time_step_table_test"

### Penalty Imposition Tests
class PenaltyImpositionBeamCantileverStaticHyperelasticSelfWeightLoad2DQuadTest(MPMTestFactory):
    file_name = "beam_tests/hyperelastic_cantilever_beam/penalty_self_weight_load_2D_quad_test"

### Slip Boundary Tests
class SlipBoundaryTest(MPMTestFactory):
    file_name = "slip_tests/slip_boundary_test"

### Explicit time integration tests
class ExplicitOscillatingPointUSLTest(MPMTestFactory):
    file_name = "explicit_tests/oscillating_point/usl_explicit_oscillating_point_test"

class ExplicitOscillatingPointUSFTest(MPMTestFactory):
    file_name = "explicit_tests/oscillating_point/usf_explicit_oscillating_point_test"

class ExplicitOscillatingPointMUSLTest(MPMTestFactory):
    file_name = "explicit_tests/oscillating_point/musl_explicit_oscillating_point_test"

class ExplicitOscillatingPointCentralDifferenceTest(MPMTestFactory):
    file_name = "explicit_tests/oscillating_point/central_difference_explicit_oscillating_point_test"

class ExplicitOscillatingPointYCompressibleTest(MPMTestFactory):
    file_name = "explicit_tests/oscillating_point/explicit_oscillating_point_Y_compressible_test"

class ExplicitOscillatingPointGravityTest(MPMTestFactory):
    file_name = "explicit_tests/oscillating_point/explicit_oscillating_point_gravity_test"

class ExplicitOscillatingPointTriTest(MPMTestFactory):
    file_name = "explicit_tests/oscillating_point/tri_explicit_oscillating_point_test"

class ExplicitAxisymDiskTriCompressibleTest(MPMTestFactory):
    file_name = "explicit_tests/axisymmetric_disk/tri_compressible_explicit_axisym_disk_test"

class ExplicitAxisymDiskQuadCompressibleTest(MPMTestFactory):
    file_name = "explicit_tests/axisymmetric_disk/quad_compressible_explicit_axisym_disk_test"

class Explicit3dHexCompressibleOscillatingPointTest(MPMTestFactory):
    file_name = "explicit_tests/oscillating_point_3d/3dhex_compressible_explicit_oscillating_point_test"

class Explicit3dTetCompressibleOscillatingPointTest(MPMTestFactory):
    file_name = "explicit_tests/oscillating_point_3d/3dtet_compressible_explicit_oscillating_point_test"

### PQMPM tests
class PQMPMExplicitQuadTest(MPMTestFactory):
    file_name = "pqmpm_tests/pqmpm_explicit_quad_test"

class PQMPMExplicitTriTest(MPMTestFactory):
    file_name = "pqmpm_tests/pqmpm_explicit_tri_test"

class PQMPMExplicitHexTest(MPMTestFactory):
    file_name = "pqmpm_tests/pqmpm_explicit_hex_test"

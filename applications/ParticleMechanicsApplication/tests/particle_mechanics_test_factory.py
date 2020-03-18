from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.ParticleMechanicsApplication.particle_mechanics_analysis import ParticleMechanicsAnalysis

class ParticleMechanicsTestFactory(KratosUnittest.TestCase):
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
class AxisSymmetricCircularPlate2DTriTest(ParticleMechanicsTestFactory):
    file_name = "axisym_tests/circular_plate_axisym_test"

### Beam Tests
class BeamCantileverStaticLinearElasticPointLoad2DTriTest(ParticleMechanicsTestFactory):
    file_name = "beam_tests/cantilever_beam/static_point_load_2D_tri_test"

class BeamCantileverStaticLinearElasticLineLoad2DQuadTest(ParticleMechanicsTestFactory):
    file_name = "beam_tests/cantilever_beam/static_line_load_2D_quad_test"

class BeamCantileverStaticLinearElasticParticlePointLoad2DTriTest(ParticleMechanicsTestFactory):
    file_name = "beam_tests/cantilever_beam/particle_point_load_2D_tri_test"

class BeamCantileverStaticLinearElasticSurfaceLoad3DHexaTest(ParticleMechanicsTestFactory):
    file_name = "beam_tests/cantilever_beam/static_surface_load_3D_hexa_test"

class BeamCantileverStaticHyperelasticSelfWeightLoad2DQuadTest(ParticleMechanicsTestFactory):
    file_name = "beam_tests/hyperelastic_cantilever_beam/self_weight_load_2D_quad_test"

class BeamCantileverLinearStaticHyperelasticSelfWeightLoad2DQuadTest(ParticleMechanicsTestFactory):
    file_name = "beam_tests/hyperelastic_cantilever_beam/linear_self_weight_load_2D_quad_test"

### Cook's Membrane Tests
class CooksMembraneCompressibleTest(ParticleMechanicsTestFactory):
    file_name = "cooks_membrane_tests/compressible_cook_membrane_2D_test"

class CooksMembraneUPCompressibleTest(ParticleMechanicsTestFactory):
    file_name = "cooks_membrane_tests/UP_compressible_cook_membrane_2D_test"

class CooksMembraneUPIncompressibleTest(ParticleMechanicsTestFactory):
    file_name = "cooks_membrane_tests/UP_incompressible_cook_membrane_2D_test"

### Constitutive Law Tests
class CLLinearElastic3DQuadTest(ParticleMechanicsTestFactory):
    file_name = "cl_tests/solid_cl/linear_elastic_3D_hexa_test"

### Gravity Application Tests
class GravityApplicationTest(ParticleMechanicsTestFactory):
    file_name = "gravity_tests/dynamic_gravity_application_test"

### Penalty Imposition Tests
class PenaltyImpositionBeamCantileverStaticHyperelasticSelfWeightLoad2DQuadTest(ParticleMechanicsTestFactory):
    file_name = "beam_tests/hyperelastic_cantilever_beam/penalty_self_weight_load_2D_quad_test"

### Slip Boundary Tests
class SlipBoundaryTest(ParticleMechanicsTestFactory):
    file_name = "slip_tests/slip_boundary_test"
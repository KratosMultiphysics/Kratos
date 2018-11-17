from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.ParticleMechanicsApplication

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as KratosUnittest
import particle_mechanics_analysis

# Other imports
import os

# This utility will control the execution scope in case we need to access files or we depend
# on specific relative locations of the files.

class controlledExecutionScope:
    def __init__(self, scope):
        self.currentPath = os.getcwd()
        self.scope = scope

    def __enter__(self):
        os.chdir(self.scope)

    def __exit__(self, type, value, traceback):
        os.chdir(self.currentPath)


class ParticleMechanicsTestFactory(KratosUnittest.TestCase):
    def setUp(self):
        # Within this location context:
        with controlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):

            # Reading the ProjectParameters
            with open(self.file_name + "_parameters.json",'r') as parameter_file:
                ProjectParameters = KratosMultiphysics.Parameters(parameter_file.read())

            self.modify_parameters(ProjectParameters)

            # To avoid many prints
            if (ProjectParameters["problem_data"]["echo_level"].GetInt() == 0):
                KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)

            # Creating the test
            model = KratosMultiphysics.Model()
            self.test = particle_mechanics_analysis.ParticleMechanicsAnalysis(model, ProjectParameters)
            self.test.Initialize()

    def modify_parameters(self, project_parameters):
        """This function can be used in derived classes to modify existing parameters
        before the execution of the test (e.g. switch to MPI)
        """
        pass

    def test_execution(self):
        # Within this location context:
        with controlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):
            self.test.RunSolutionLoop()

    def tearDown(self):
        # Within this location context:
        with controlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):
            self.test.Finalize()

### Beam Tests
class BeamCantileverStaticLinearElasticPointLoad2DTriTest(ParticleMechanicsTestFactory):
    file_name = "beam_tests/cantilever_beam/static_point_load_2D_tri_test"

class BeamCantileverStaticLinearElasticLineLoad2DQuadTest(ParticleMechanicsTestFactory):
    file_name = "beam_tests/cantilever_beam/static_line_load_2D_quad_test"

class BeamCantileverStaticLinearElasticSurfaceLoad3DHexaTest(ParticleMechanicsTestFactory):
    file_name = "beam_tests/cantilever_beam/static_surface_load_3D_hexa_test"

class BeamCantileverStaticHyperelasticSelfWeightLoad2DQuadTest(ParticleMechanicsTestFactory):
    file_name = "beam_tests/hyperelastic_cantilever_beam/self_weight_load_2D_quad_test"

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

### Slip Boundary Tests
class SlipBoundaryTest(ParticleMechanicsTestFactory):
    file_name = "slip_tests/slip_boundary_test"
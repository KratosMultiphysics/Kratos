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

class BeamCantileverLinearElasticPointLoad2DTriTest(ParticleMechanicsTestFactory):
    file_name = "beam_tests/cantilever_beam/linear_elastic_point_load_2d_tri_test"

class BeamCantileverLinearElasticPointLoad2DQuadTest(ParticleMechanicsTestFactory):
    file_name = "cl_tests/solid_cl/linear_elastic_2D_hexa_test"

class BeamCantileverLinearElasticLineLoad2DTriTest(ParticleMechanicsTestFactory):
    file_name = "cl_tests/solid_cl/linear_elastic_2D_hexa_test"

class BeamCantileverLinearElasticLineLoad2DQuadTest(ParticleMechanicsTestFactory):
    file_name = "cl_tests/solid_cl/linear_elastic_2D_hexa_test"

class BeamCantileverLinearElasticSurfaceLoad3DTetraTest(ParticleMechanicsTestFactory):
    file_name = "cl_tests/solid_cl/linear_elastic_2D_hexa_test"

class BeamCantileverLinearElasticSurfaceLoad3DHexaTest(ParticleMechanicsTestFactory):
    file_name = "cl_tests/solid_cl/linear_elastic_2D_hexa_test"

class CLLinearElastic2DQuadTest(ParticleMechanicsTestFactory):
    file_name = "cl_tests/solid_cl/linear_elastic_2D_hexa_test"
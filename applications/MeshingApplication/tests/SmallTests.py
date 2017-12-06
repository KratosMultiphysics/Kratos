import os

# Import Kratos
from KratosMultiphysics import *

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as KratosUnittest
import Kratos_Execute_Meshing_Test as Execute_Test

# This utiltiy will control the execution scope in case we need to acces files or we depend
# on specific relative locations of the files.

# TODO: Should we move this to KratosUnittest?
class controlledExecutionScope:
    def __init__(self, scope):
        self.currentPath = os.getcwd()
        self.scope = scope

    def __enter__(self):
        os.chdir(self.scope)

    def __exit__(self, type, value, traceback):
        os.chdir(self.currentPath)

class MeshingTestFactory(KratosUnittest.TestCase):

    def setUp(self):
        # Within this location context:
        with controlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):
            # Initialize GiD  I/O
            parameter_file = open(self.file_name + "_parameters.json", 'r')
            ProjectParameters = Parameters(parameter_file.read())

            # Creating the model part
            self.test = Execute_Test.Kratos_Execute_Test(ProjectParameters)

    def test_execution(self):
        # Within this location context:
        with controlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):
            self.test.Solve()

    def tearDown(self):
        pass
    
#class TwoDLevelSetTest(MeshingTestFactory):
    #file_name = "mmg_eulerian_test/2D_levelset_test"
    
#class ThreeDLevelSetTest(MeshingTestFactory):
    #file_name = "mmg_eulerian_test/3D_levelset_test"
    
class TwoDHessianTest(MeshingTestFactory):
    file_name = "mmg_eulerian_test/2D_hessian_test"
    
class ThreeDHessianTest(MeshingTestFactory):
    file_name = "mmg_eulerian_test/3D_hessian_test"
    
class TwoDCavityTest(MeshingTestFactory):
    file_name = "mmg_eulerian_test/2D_cavity_test"
    
class TwoDDynamicBeamTest(MeshingTestFactory):
    file_name = "mmg_lagrangian_test/beam2D_test"
    
class TwoDDynamicBeamLineLoadTest(MeshingTestFactory):
    file_name = "mmg_lagrangian_test/beam2D_line_load_test"
    
class ThreeDDynamicBeamTest(MeshingTestFactory):
    file_name = "mmg_lagrangian_test/beam3D_test"
    
class TwoDDynamicPlasticBeamTest(MeshingTestFactory):
    file_name = "mmg_lagrangian_test/beam2D_internal_variables_interpolation_test"
    

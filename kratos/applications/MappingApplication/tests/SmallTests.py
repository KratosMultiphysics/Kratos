import os

# Import Kratos
from KratosMultiphysics import *

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as KratosUnittest
from NearestNeighborMapperTest import NearestNeighborMapperTest

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


class NearestNeighborMapperTestFactory(KratosUnittest.TestCase):

    def setUp(self):
        # Within this location context:
        with controlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):
            output_post = True # set to "True" if GiD output is wanted
            self.nearest_neighbor_mapper_test = NearestNeighborMapperTest(output_post)

    def test_execution(self):
        # Within this location context:
        with controlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):
            # the numeric values are the output times for GiD
            self.nearest_neighbor_mapper_test.TestMapConstantScalarValues(1.0)
            self.nearest_neighbor_mapper_test.TestInverseMapConstantScalarValues(2.0)

            self.nearest_neighbor_mapper_test.TestMapConstantVectorValues(3.0)
            self.nearest_neighbor_mapper_test.TestInverseMapConstantVectorValues(4.0)

            self.nearest_neighbor_mapper_test.TestMapNonConstantScalarValues(5.0)
            self.nearest_neighbor_mapper_test.TestInverseMapNonConstantScalarValues(6.0)

            self.nearest_neighbor_mapper_test.TestMapNonConstantVectorValues(7.0)
            self.nearest_neighbor_mapper_test.TestInverseMapNonConstantVectorValues(8.0)

    def tearDown(self):
        pass


class NearestNeighborTest_1(NearestNeighborMapperTestFactory):
    file_name = "Mapper_Test_1/Mapper_Test_1"

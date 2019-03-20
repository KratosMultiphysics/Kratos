import os

# Import Kratos
from KratosMultiphysics import *

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as KratosUnittest
from NearestNeighborMapperTest import NearestNeighborMapperTest
from NearestElementMapperTest2D import NearestElementMapperTest2D
import KratosExecuteMapperTests as ExecuteMapperTests

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
        self.execution_path = os.path.dirname(os.path.realpath(__file__))
        # Within this location context:
        with controlledExecutionScope(self.execution_path):
            output_post = True # set to "True" if GiD output is wanted
            self.nearest_neighbor_mapper_test = NearestNeighborMapperTest(output_post)

    def test_execution(self):
        # Within this location context:
        with controlledExecutionScope(self.execution_path):
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



class NearestElementMapperTest2DFactory(KratosUnittest.TestCase):

    def setUp(self):
        self.execution_path = os.path.dirname(os.path.realpath(__file__))
        # Within this location context:
        with controlledExecutionScope(self.execution_path):
            output_post = True # set to "True" if GiD output is wanted
            self.nearest_element_mapper_test_2D = NearestElementMapperTest2D(output_post)

    def test_execution(self):
        # Within this location context:
        with controlledExecutionScope(self.execution_path):
            # the numeric values are the output times for GiD
            self.nearest_element_mapper_test_2D.TestMapConstantScalarValues(1.0)
            self.nearest_element_mapper_test_2D.TestInverseMapConstantScalarValues(2.0)

            self.nearest_element_mapper_test_2D.TestMapConstantVectorValues(3.0)
            self.nearest_element_mapper_test_2D.TestInverseMapConstantVectorValues(4.0)

            self.nearest_element_mapper_test_2D.TestMapNonConstantScalarValues(5.0)
            self.nearest_element_mapper_test_2D.TestInverseMapNonConstantScalarValues(6.0)

            self.nearest_element_mapper_test_2D.TestMapNonConstantVectorValues(7.0)
            self.nearest_element_mapper_test_2D.TestInverseMapNonConstantVectorValues(8.0)

    def tearDown(self):
        pass

class MapperTestsFactory(KratosUnittest.TestCase):
    def setUp(self):
        self.execution_path = os.path.dirname(os.path.realpath(__file__))
        # Within this location context:
        with controlledExecutionScope(self.execution_path):
            self.output_post = False # set to "True" if GiD output is wanted
            self.set_up_test_1 = False # set to "True" to print the coordinates and the prescribed values
            self.set_up_test_2 = False # set to "True" to print the mapped Values
            self.test_object = ExecuteMapperTests.KratosExecuteMapperTests(self.output_post,
                                                                           self.set_up_test_1,
                                                                           self.set_up_test_2)

    def test_execution(self):
        # Within this location context:
        with controlledExecutionScope(self.execution_path):
            # the numeric values are the output times for GiD
            for file_name in self.file_name_list:
                print("Testing Test \"" + file_name + "\" ... ", end='') # this is only printed in case the test fails
                self.test_object.SetUpMapper(file_name)

                self.test_object.TestMapConstantScalarValues(1.0)
                self.test_object.TestInverseMapConstantScalarValues(2.0)

                self.test_object.TestMapConstantVectorValues(3.0)
                self.test_object.TestInverseMapConstantVectorValues(4.0)

                self.test_object.TestMapNonConstantScalarValues(5.0)
                self.test_object.TestInverseMapNonConstantScalarValues(6.0)

                self.test_object.TestMapNonConstantVectorValues(7.0)
                self.test_object.TestInverseMapNonConstantVectorValues(8.0)

                print("succssful") # this is only printed in case the test fails

                if (self.output_post):
                    self.test_object.FinalizeGiD()

    def tearDown(self):
        if (self.set_up_test_1 or self.set_up_test_2):
            err # needed to get the output
        else:
            pass

class NearestNeighborTest_1(NearestNeighborMapperTestFactory):
    file_name = "Mapper_Test_1/Mapper_Test_1"

class NearestElementTest2D_1(NearestElementMapperTest2DFactory):
    file_name = "Mapper_Test_1/Mapper_Test_1"

class MapperTests(MapperTestsFactory):
    # Add new tests here
    file_name_1 = "NearestNeighbor_line"
    file_name_2 = "NearestNeighbor_surface"
    file_name_3 = "NearestNeighbor_volume"
    file_name_4 = "NearestElement_line"
    file_name_5 = "NearestElement_surface"
    file_name_6 = "NearestElement_volume"

    file_name_list = [file_name_1, file_name_2, file_name_3, file_name_4, file_name_5, file_name_6]
    file_name_list = [file_name_1, file_name_2, file_name_3, file_name_4, file_name_5]

from __future__ import print_function, absolute_import, division

import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest

class TestVariables(KratosUnittest.TestCase):
    def test_GetSourceVariable(self):
        component_variable = KratosMultiphysics.VELOCITY_X
        source_variable = KratosMultiphysics.VELOCITY
        false_source_variable = KratosMultiphysics.DISPLACEMENT
        self.assertEqual(component_variable.GetSourceVariable(), source_variable)
        self.assertNotEqual(component_variable.GetSourceVariable(), false_source_variable)

    def test_CreateVariable(self):
        NEW_STRING_VARIABLE = KratosMultiphysics.StringVariable("NEW_STRING_VARIABLE")
        NEW_TEST_BOOL_VARIABLE = KratosMultiphysics.BoolVariable("NEW_TEST_BOOL_VARIABLE")
        NEW_INT_VARIABLE = KratosMultiphysics.IntegerVariable("NEW_INT_VARIABLE")
        NEW_INT_VECTOR_VARIABLE = KratosMultiphysics.IntegerVectorVariable("NEW_INT_VECTOR_VARIABLE")
        NEW_DOUBLE_VARIABLE = KratosMultiphysics.DoubleVariable("NEW_DOUBLE_VARIABLE")
        NEW_DOUBLE_VECTOR_VARIABLE = KratosMultiphysics.VectorVariable("NEW_DOUBLE_VECTOR_VARIABLE")
        NEW_VARIABLE = KratosMultiphysics.Array1DVariable3("NEW_VARIABLE")

    def test_CreateDuplicateVariableFail(self):
        #This test is to show that the test should not fail when creating a variable with same name.
        NEW_STRING_VARIABLE_TWO = KratosMultiphysics.StringVariable("NEW_STRING_VARIABLE_TWO")
        NEW_STRING_VARIABLE_TWO = KratosMultiphysics.StringVariable("NEW_STRING_VARIABLE_TWO")

    def test_CreateDuplicateVariableComponentFail(self):
        #This test is to show that the test should not fail when creating a variable with same name.
        NEW_DOUBLE_VECTOR_VARIABLE_TWO = KratosMultiphysics.Array1DVariable3("NEW_DOUBLE_VECTOR_VARIABLE_TWO")
        NEW_DOUBLE_VECTOR_VARIABLE_TWO_X = KratosMultiphysics.Array1DComponentVariable("NEW_DOUBLE_VECTOR_VARIABLE_TWO_X", "NEW_DOUBLE_VECTOR_VARIABLE_TWO", 0)

    def test_CreateNodesWithCreatedVariable(self):
        #This test is to show that the test should not fail when using new variables in a ModelPart
        current_model = KratosMultiphysics.Model()
        model_part = current_model.CreateModelPart("for_test")

        new_double_variable = KratosMultiphysics.DoubleVariable("NEW_DOUBLE_VARIABLE")
        new_array_3_variable = KratosMultiphysics.Array1DVariable3("NEW_VARIABLE")

        model_part.AddNodalSolutionStepVariable(new_double_variable)
        model_part.AddNodalSolutionStepVariable(new_array_3_variable)

        new_array_3_variable = KratosMultiphysics.KratosGlobals.GetVariable("DISPLACEMENT")
        model_part.AddNodalSolutionStepVariable(new_array_3_variable)

        new_node = model_part.CreateNewNode(1, 0.0, 0.0, 0.0)

        new_node.SetSolutionStepValue(new_double_variable, 0, 1.5)
        new_node.SetSolutionStepValue(new_array_3_variable, 0, [1.50, 0.8, -9.147])


if __name__ == '__main__':
    KratosUnittest.main()

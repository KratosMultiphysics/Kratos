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

        NEW_Q_VARIABLE = KratosMultiphysics.DoubleQuaternionVariable("NEW_Q_VARIABLE")

        NEW_VARIABLE = KratosMultiphysics.Array1DVariable3("NEW_VARIABLE")


    @KratosUnittest.expectedFailure
    def test_CreateDuplicateVariableFail(self):
        NEW_STRING_VARIABLE_TWO = KratosMultiphysics.StringVariable("NEW_STRING_VARIABLE_TWO")
        NEW_STRING_VARIABLE_TWO = KratosMultiphysics.StringVariable("NEW_STRING_VARIABLE_TWO")

    @KratosUnittest.expectedFailure
    def test_CreateDuplicateVariableComponentFail(self):
        NEW_VARIABLE_TWO = KratosMultiphysics.Array1DVariable3("NEW_VARIABLE_TWO")
        NEW_VARIABLE_TWO_X = KratosMultiphysics.Array1DComponentVariable("NEW_VARIABLE_TWO_X", "NEW_VARIABLE_TWO", 1)



if __name__ == '__main__':
    KratosUnittest.main()

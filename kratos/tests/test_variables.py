from __future__ import print_function, absolute_import, division

import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest

class TestVariables(KratosUnittest.TestCase):
    def test_01_GetSourceVariable(self): # Tests with numbers to force the order
        component_variable = KratosMultiphysics.VELOCITY_X
        source_variable = KratosMultiphysics.VELOCITY
        false_source_variable = KratosMultiphysics.DISPLACEMENT
        self.assertEqual(component_variable.GetSourceVariable(), source_variable)
        self.assertNotEqual(component_variable.GetSourceVariable(), false_source_variable)

    def test_02_VariableComponentConstructor(self):
        NEW_VARIABLE = KratosMultiphysics.Array1DVariable3("NEW_VARIABLE")

        NEW_VARIABLE_X = KratosMultiphysics.Array1DComponentVariable("NEW_VARIABLE_X", "NEW_VARIABLE", 1)
        NEW_VARIABLE_Y = KratosMultiphysics.Array1DComponentVariable("NEW_VARIABLE_Y", "NEW_VARIABLE", 2)
        NEW_VARIABLE_Z = KratosMultiphysics.Array1DComponentVariable("NEW_VARIABLE_Z", "NEW_VARIABLE", 3)

    def test_03_CreateVariable(self):
        NEW_STRING_VARIABLE = KratosMultiphysics.StringVariable("NEW_STRING_VARIABLE")

        NEW_TEST_BOOL_VARIABLE = KratosMultiphysics.BoolVariable("NEW_TEST_BOOL_VARIABLE")

        NEW_INT_VARIABLE = KratosMultiphysics.IntegerVariable("NEW_INT_VARIABLE")

        NEW_INT_VECTOR_VARIABLE = KratosMultiphysics.IntegerVectorVariable("NEW_INT_VECTOR_VARIABLE")

        NEW_DOUBLE_VARIABLE = KratosMultiphysics.DoubleVariable("NEW_DOUBLE_VARIABLE")

        NEW_DOUBLE_VECTOR_VARIABLE = KratosMultiphysics.VectorVariable("NEW_DOUBLE_VECTOR_VARIABLE")

        NEW_Q_VARIABLE = KratosMultiphysics.DoubleQuaternionVariable("NEW_Q_VARIABLE")

    @KratosUnittest.expectedFailure
    def test_04_CreateDuplicateVariableFail(self):
        NEW_STRING_VARIABLE = KratosMultiphysics.StringVariable("NEW_STRING_VARIABLE")

    @KratosUnittest.expectedFailure
    def test_05_CreateDuplicateVariableComponentFail(self):
        NEW_VARIABLE_TWO = KratosMultiphysics.Array1DVariable3("NEW_VARIABLE_TWO")
        NEW_VARIABLE_TWO_X = KratosMultiphysics.Array1DComponentVariable("NEW_VARIABLE_TWO_X", "NEW_VARIABLE_TWO", 1)
        NEW_VARIABLE_TWO_Y = KratosMultiphysics.Array1DComponentVariable("NEW_VARIABLE_TWO_X", "NEW_VARIABLE_TWO", 1)


if __name__ == '__main__':
    KratosUnittest.main()

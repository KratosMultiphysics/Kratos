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

    def test_VariableComponentConstructor(self):
        new_model_part = KratosMultiphysics.ModelPart("NewModelPart")
        NEW_VARIABLE = KratosMultiphysics.Array1DVariable3("NEW_VARIABLE")
        new_model_part.AddNodalSolutionStepVariable(NEW_VARIABLE)

        NEW_VARIABLE_X = KratosMultiphysics.Array1DComponentVariable("NEW_VARIABLE_X", "NEW_VARIABLE", 1)
        NEW_VARIABLE_Y = KratosMultiphysics.Array1DComponentVariable("NEW_VARIABLE_Y", "NEW_VARIABLE", 2)
        NEW_VARIABLE_Z = KratosMultiphysics.Array1DComponentVariable("NEW_VARIABLE_Z", "NEW_VARIABLE", 3)

    def test_CreateVariable(self):
        new_model_part = KratosMultiphysics.ModelPart("NewModelPart")
        NEW_STRING_VARIABLE = KratosMultiphysics.StringVariable("NEW_STRING_VARIABLE")

        NEW_BOOL_VARIABLE = KratosMultiphysics.BoolVariable("NEW_BOOL_VARIABLE")
        new_model_part.AddNodalSolutionStepVariable(NEW_BOOL_VARIABLE)

        NEW_INT_VARIABLE = KratosMultiphysics.IntegerVariable("NEW_INT_VARIABLE")
        new_model_part.AddNodalSolutionStepVariable(NEW_INT_VARIABLE)

        NEW_INT_VECTOR_VARIABLE = KratosMultiphysics.IntegerVectorVariable("NEW_INT_VECTOR_VARIABLE")

        NEW_DOUBLE_VARIABLE = KratosMultiphysics.DoubleVariable("NEW_DOUBLE_VARIABLE")
        new_model_part.AddNodalSolutionStepVariable(NEW_DOUBLE_VARIABLE)

        NEW_DOUBLE_VECTOR_VARIABLE = KratosMultiphysics.VectorVariable("NEW_DOUBLE_VECTOR_VARIABLE")
        new_model_part.AddNodalSolutionStepVariable(NEW_DOUBLE_VECTOR_VARIABLE)

    def test_CreateDuplicateVariableFail(self):
        NEW_STRING_VARIABLE = KratosMultiphysics.StringVariable("NEW_STRING_VARIABLE")
        NEW_STRING_VARIABLE = KratosMultiphysics.StringVariable("NEW_STRING_VARIABLE")

    def test_CreateDuplicateVariableFail(self):
        NEW_VARIABLE = KratosMultiphysics.Array1DVariable3("NEW_VARIABLE")
        NEW_VARIABLE_X = KratosMultiphysics.Array1DComponentVariable("NEW_VARIABLE_X", "NEW_VARIABLE", 1)
        NEW_VARIABLE_Y = KratosMultiphysics.Array1DComponentVariable("NEW_VARIABLE_X", "NEW_VARIABLE", 1)


if __name__ == '__main__':
    KratosUnittest.main()

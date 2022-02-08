import sympy

import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.sympy_fe_utilities as KratosSympy


class TestSympyFEUtilities(KratosUnittest.TestCase):

    def test_ScalarOutput(self):
        x = sympy.var('x')
        y = sympy.var('y')
        z = sympy.var('z')

        expression = x+y/z**3
        code = KratosSympy.OutputScalar(expression, "f", "python", 2, False, " += ")

        expected_code =  "        f += x + y/z**3\n"

        self.assertEqual(code, expected_code)

    def test_VectorOutput(self):
        x = KratosSympy.DefineVector('x', 3)
        y = KratosSympy.DefineVector('y', 3)

        expression = x+y
        code = KratosSympy.OutputVector(expression, "f", "python", 1, False, " += ")

        expected_code = "\n".join([
            "    f[0] += x_0 + y_0",
            "    f[1] += x_1 + y_1",
            "    f[2] += x_2 + y_2",
            ""
        ])

        self.assertEqual(code, expected_code)

    def test_MatrixOutput(self):
        x = KratosSympy.DefineMatrix('x', 2, 3)
        y = KratosSympy.DefineMatrix('y', 2, 3)

        expression = x-y
        code = KratosSympy.OutputMatrix(expression, "f", "c", 0, False, " = ")

        expected_code = "\n".join([
            "f(0,0) = x_0_0 - y_0_0;",
            "f(0,1) = x_0_1 - y_0_1;",
            "f(0,2) = x_0_2 - y_0_2;",
            "f(1,0) = x_1_0 - y_1_0;",
            "f(1,1) = x_1_1 - y_1_1;",
            "f(1,2) = x_1_2 - y_1_2;",
            ""
        ])

        self.assertEqual(code, expected_code)

    def test_OutputSymbolicVariableDeclaration(self):
        x = sympy.var('x')
        y = sympy.var('y')

        expression = sympy.sqrt(x**2 + y**2)

        code = KratosSympy.OutputSymbolicVariableDeclaration(expression, "f", "c", 1, False)

        expected_code = "    const double f = sqrt(pow(x, 2) + pow(y, 2));\n"
        self.assertEqual(code, expected_code)


if __name__ == '__main__':
    KratosUnittest.main()

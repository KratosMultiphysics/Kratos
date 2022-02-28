import sympy

import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.sympy_fe_utilities as KratosSympy


class TestSympyFEUtilities(KratosUnittest.TestCase):

    def testScalarOutput(self):
        x = sympy.var('x')
        y = sympy.var('y')
        z = sympy.var('z')

        expression = x+y/z**3
        code = KratosSympy.OutputScalar(expression, "f", "python", 2, False, " += ")

        expected_code =  "        f += x + y/z**3\n"

        self.assertEqual(code, expected_code)

    def testVectorOutput(self):
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

    def testMatrixOutput(self):
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

    def testOutputSymbolicVariableDeclaration(self):
        x = sympy.var('x')
        y = sympy.var('y')

        expression = sympy.sqrt(x**2 + y**2)

        code = KratosSympy.OutputSymbolicVariableDeclaration(expression, "f", "c", 1, False)

        expected_code = "    const double f = sqrt(pow(x, 2) + pow(y, 2));\n"
        self.assertEqual(code, expected_code)

    def testOutputScalar_CollectingFactors(self):
        x = KratosSympy.DefineVector('x', 5)
        y = KratosSympy.DefineVector('y', 5)

        expression = (3*(x.T*(2*y + 3*x))*(y.T*y) + 6*y.T*x)[0, 0]
        code = KratosSympy.OutputScalar_CollectingFactors(expression, "myvariable", "c", 1, replace_indices=True)

        self.maxDiff = 2048
        expected_code = "\n".join([
            "    const double cmyvariable0 = 2*y[0];",
            "    const double cmyvariable1 = 2*y[1];",
            "    const double cmyvariable2 = 2*y[2];",
            "    const double cmyvariable3 = 2*y[3];",
            "    const double cmyvariable4 = 2*y[4];",
            "    myvariable=3*cmyvariable0*x[0] + 3*cmyvariable1*x[1] + 3*cmyvariable2*x[2] + 3*cmyvariable3*x[3] + 3*cmyvariable4*x[4] + 3*(pow(y[0], 2) + pow(y[1], 2) + pow(y[2], 2) + pow(y[3], 2) + pow(y[4], 2))*(x[0]*(cmyvariable0 + 3*x[0]) + x[1]*(cmyvariable1 + 3*x[1]) + x[2]*(cmyvariable2 + 3*x[2]) + x[3]*(cmyvariable3 + 3*x[3]) + x[4]*(cmyvariable4 + 3*x[4]));",
            ""
        ])
        self.assertEqual(code, expected_code)

    def testDefineSymmetricMatrix(self):
        nrows = 5
        ncols = 7
        A = KratosSympy.DefineSymmetricMatrix("A", nrows, ncols)

        square_region = min(nrows, ncols)
        A_square = A[:square_region, :square_region]

        self.assertEqual(A_square.T, A_square)

if __name__ == '__main__':
    KratosUnittest.main()

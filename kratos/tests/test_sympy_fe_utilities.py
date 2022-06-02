import sympy

import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.sympy_fe_utilities as KratosSympy


class TestSympyFEUtilities(KratosUnittest.TestCase):

    def testConvertVoigtMatrixToTensor2D(self):
        dim = 2
        C_voigt = sympy.Matrix([
            [1,2,3],
            [2,4,5],
            [3,5,6]])
        C = KratosSympy.ConvertVoigtMatrixToTensor(C_voigt)
        C_reference = [
            [[[1,3],[3,2]],[[3,6],[6,5]]],
            [[[3,6],[6,5]],[[2,5],[5,4]]]]
        for i in range(dim):
            for j in range(dim):
                for k in range(dim):
                    for l in range(dim):
                        self.assertEqual(C[i][j][k][l], C_reference[i][j][k][l])

    def testConvertVoigtMatrixToTensor3D(self):
        dim = 3
        C_voigt = sympy.Matrix([
            [1,2,3,4,5,6],
            [2,7,8,9,10,11],
            [3,8,12,13,14,15],
            [4,9,13,16,17,18],
            [5,10,14,17,19,20],
            [6,11,15,18,20,21]])
        C = KratosSympy.ConvertVoigtMatrixToTensor(C_voigt)
        C_reference = [
            [[[1,6,5],[6,2,4],[5,4,3]],[[6,21,20],[21,11,18],[20,18,15]],[[5,20,19],[20,10,17],[19,17,14]]],
            [[[6,21,20],[21,11,18],[20,18,15]],[[2,11,10],[11,7,9],[10,9,8]],[[4,18,17],[18,9,16],[17,16,13]]],
            [[[5,20,19],[20,10,17],[19,17,14]],[[4,18,17],[18,9,16],[17,16,13]],[[3,15,14],[15,8,13],[14,13,12]]]]
        for i in range(dim):
            for j in range(dim):
                for k in range(dim):
                    for l in range(dim):
                        self.assertEqual(C[i][j][k][l], C_reference[i][j][k][l])

    def testDoubleContraction(self):
        A_2nd_order = sympy.Matrix([
            [1,2],
            [2,3]])
        B_4th_order = sympy.MutableDenseNDimArray([
            [[[1,3],[3,2]],[[3,6],[6,5]]],
            [[[3,6],[6,5]],[[2,5],[5,4]]]], shape=(2,2,2,2))
        C_2nd_order = KratosSympy.DoubleContraction(A_2nd_order,B_4th_order)
        C_2nd_order_reference = [
            [19,42],
            [42,34]]
        for i in range(2):
            for j in range(2):
                self.assertEqual(C_2nd_order[i][j], C_2nd_order_reference[i][j])

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

    def testDefineSymmetricFourthOrderTensor(self):
        A = KratosSympy.DefineSymmetricFourthOrderTensor("A",3,3,3,3)
        self.assertEqual(A.rank(), 4)
        for i in range(3):
            for j in range(3):
                aux_i = min(i,j)
                aux_j = max(i,j)
                for k in range(3):
                    for l in range(3):
                        aux_k = min(k,l)
                        aux_l = max(k,l)
                        self.assertEqual(A[i,j,k,l].name, "A_{}_{}_{}_{}".format(aux_i,aux_j,aux_k,aux_l))

if __name__ == '__main__':
    KratosUnittest.main()

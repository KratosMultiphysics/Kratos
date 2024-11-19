import sympy

import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.sympy_fe_utilities as KratosSympy


class TestSympyFEUtilities(KratosUnittest.TestCase):

    def testDefineShapeFunctions1(self):
        n_dim = 2
        n_nodes = 3
        partition_of_unity = False
        N, DN = KratosSympy.DefineShapeFunctions(n_nodes, n_dim, partition_of_unity)

        for i in range(n_nodes):
            self.assertEqual(N[i], sympy.var(f"N_{i}"))
            for j in range(n_dim):
                self.assertEqual(DN[i,j], sympy.var(f"DN_{i}_{j}"))

    def testDefineTestFunctions2(self):
        n_dim = 2
        n_nodes = 3
        N, DN = KratosSympy.DefineShapeFunctions(n_nodes, n_dim, shape_functions_name="my_N", first_derivatives_name="my_DN")

        for i in range(n_nodes):
            self.assertEqual(N[i], sympy.var(f"my_N_{i}"))
            for j in range(n_dim):
                self.assertEqual(DN[i,j], sympy.var(f"my_DN_{i}_{j}"))

    def testDefineTestFunctions3(self):
        n_dim = 2
        n_nodes = 6
        N, DN, DDN = KratosSympy.DefineShapeFunctions(n_nodes, n_dim, shape_functions_name="my_N", first_derivatives_name="my_DN", second_derivatives_name="my_DDN")

        for i in range(n_nodes):
            self.assertEqual(N[i], sympy.var(f"my_N_{i}"))
            for j in range(n_dim):
                self.assertEqual(DN[i,j], sympy.var(f"my_DN_{i}_{j}"))
                for k in range(n_dim):
                    self.assertEqual(DDN[i][j,k], sympy.var(f"my_DDN_{i}_{j}_{k}"))

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
            [[[1,4,6],[4,2,5],[6,5,3]],[[4,16,18],[16,9,17],[18,17,13]],[[6,18,21],[18,11,20],[21,20,15]]],
            [[[4,16,18],[16,9,17],[18,17,13]],[[2,9,11],[9,7,10],[11,10,8]],[[5,17,20],[17,10,19],[20,19,14]]],
            [[[6,18,21],[18,11,20],[21,20,15]],[[5,17,20],[17,10,19],[20,19,14]],[[3,13,15],[13,8,14],[15,14,12]]]]
        for i in range(dim):
            for j in range(dim):
                for k in range(dim):
                    for l in range(dim):
                        self.assertEqual(C[i][j][k][l], C_reference[i][j][k][l])

    def testStrainToVoigt2D(self):
        strain = sympy.Matrix([
            [1.0, 2.0],
            [2.0, 3.0]])
        strain_voigt = KratosSympy.StrainToVoigt(strain)
        self._AuxCheckToVoigtResults(3, lambda i,j : 1.0 if i == j else 2.0, strain, strain_voigt)

    def testStrainToVoigt3D(self):
        strain = sympy.Matrix([
            [1.0, 2.0, 3.0],
            [2.0, 4.0, 5.0],
            [3.0, 5.0, 6.0]])
        strain_voigt = KratosSympy.StrainToVoigt(strain)
        self._AuxCheckToVoigtResults(6, lambda i,j : 1.0 if i == j else 2.0, strain, strain_voigt)

    def testMatrixToVoigt2D(self):
        matrix = sympy.Matrix([
            [1.0, 2.0],
            [2.0, 3.0]])
        matrix_voigt = KratosSympy.MatrixToVoigt(matrix)
        self._AuxCheckToVoigtResults(3, lambda i,j : 1.0, matrix, matrix_voigt)

    def testMatrixToVoigt3D(self):
        matrix = sympy.Matrix([
            [1.0, 2.0, 3.0],
            [2.0, 4.0, 5.0],
            [3.0, 5.0, 6.0]])
        matrix_voigt = KratosSympy.MatrixToVoigt(matrix)
        self._AuxCheckToVoigtResults(6, lambda i,j : 1.0, matrix, matrix_voigt)

    def testVoigtToMatrix2D(self):
        vector_voigt = sympy.Matrix([1.0,2.0,3.0])
        matrix = KratosSympy.VoigtToMatrix(vector_voigt)
        self._AuxCheckToMatrixResults(2, lambda i : 1.0, vector_voigt, matrix)

    def testVoigtToMatrix3D(self):
        vector_voigt = sympy.Matrix([1.0,2.0,3.0,4.0,5.0,6.0])
        matrix = KratosSympy.VoigtToMatrix(vector_voigt)
        self._AuxCheckToMatrixResults(3, lambda i : 1.0, vector_voigt, matrix)

    def testStrainToMatrix2D(self):
        strain_voigt = sympy.Matrix([1.0,2.0,3.0])
        strain_matrix = KratosSympy.StrainToMatrix(strain_voigt)
        self._AuxCheckToMatrixResults(2, lambda i : 1.0 if i < 2 else 0.5, strain_voigt, strain_matrix)

    def testStrainToMatrix3D(self):
        strain_voigt = sympy.Matrix([1.0,2.0,3.0,4.0,5.0,6.0])
        strain_matrix = KratosSympy.StrainToMatrix(strain_voigt)
        self._AuxCheckToMatrixResults(3, lambda i : 1.0 if i < 3 else 0.5, strain_voigt, strain_matrix)

    def _AuxCheckToVoigtResults(self, strain_size, factor_calculator, reference, voigt_output):
        dim = 2 if strain_size == 3 else 3
        aux_indices = KratosSympy._GetVoigtToTensorConversionIndices(strain_size)
        for i in range(dim):
            for j in range(dim):
                factor = factor_calculator(i,j)
                self.assertEqual(voigt_output[aux_indices[(i, j)]], factor * reference[i,j])

    def _AuxCheckToMatrixResults(self, dimension, factor_calculator, reference, matrix_output):
        strain_size = 3 if dimension == 2 else 6
        aux_indices = KratosSympy._GetTensorToVoigtConversionIndices(dimension)
        for i in range(strain_size):
            factor = factor_calculator(i)
            self.assertEqual(matrix_output[aux_indices[i][0],aux_indices[i][1]], factor * reference[i])

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

    def testListOfMatricesOutput(self):
        x = sympy.Matrix(1,2, lambda _, j : (sympy.Matrix(3,3, lambda m, n : sympy.var(f"x_{j}_{m}_{n}"))))

        code_0 = KratosSympy.OutputMatrix(x[0], "f0", "c", 0, False, " = ")
        code_1 = KratosSympy.OutputMatrix(x[1], "f1", "c", 0, False, " = ")

        expected_code_0 = "\n".join([
            "f0(0,0) = x_0_0_0;",
            "f0(0,1) = x_0_0_1;",
            "f0(0,2) = x_0_0_2;",
            "f0(1,0) = x_0_1_0;",
            "f0(1,1) = x_0_1_1;",
            "f0(1,2) = x_0_1_2;",
            "f0(2,0) = x_0_2_0;",
            "f0(2,1) = x_0_2_1;",
            "f0(2,2) = x_0_2_2;",
            ""
        ])
        expected_code_1 = "\n".join([
            "f1(0,0) = x_1_0_0;",
            "f1(0,1) = x_1_0_1;",
            "f1(0,2) = x_1_0_2;",
            "f1(1,0) = x_1_1_0;",
            "f1(1,1) = x_1_1_1;",
            "f1(1,2) = x_1_1_2;",
            "f1(2,0) = x_1_2_0;",
            "f1(2,1) = x_1_2_1;",
            "f1(2,2) = x_1_2_2;",
            ""
        ])

        self.assertEqual(code_0, expected_code_0)
        self.assertEqual(code_1, expected_code_1)

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

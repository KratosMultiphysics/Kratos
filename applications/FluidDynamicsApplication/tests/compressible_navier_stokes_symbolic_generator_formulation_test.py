import os
import sympy

import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnitTest

from KratosMultiphysics.FluidDynamicsApplication.symbolic_generation.compressible_navier_stokes \
    .compressible_navier_stokes_symbolic_generator import CompressibleNavierStokesSymbolicGenerator

from KratosMultiphysics.FluidDynamicsApplication.symbolic_generation.compressible_navier_stokes.src import generate_convective_flux
from KratosMultiphysics.FluidDynamicsApplication.symbolic_generation.compressible_navier_stokes.src import generate_diffusive_flux
from KratosMultiphysics.FluidDynamicsApplication.symbolic_generation.compressible_navier_stokes.src import generate_stabilization_matrix

from KratosMultiphysics.FluidDynamicsApplication.symbolic_generation.compressible_navier_stokes \
    .src.defines import CompressibleNavierStokesDefines as defs
from KratosMultiphysics.FluidDynamicsApplication.symbolic_generation.compressible_navier_stokes.src \
    .symbolic_parameters import FormulationParameters, PrimitiveMagnitudes, ShockCapturingParameters
from KratosMultiphysics.FluidDynamicsApplication.symbolic_generation.compressible_navier_stokes.src \
    .symbolic_geometry import GeometryDataFactory

import KratosMultiphysics.FluidDynamicsApplication

class CompressibleNavierStokesSymbolicGeneratorValidationTest(KratosUnitTest.TestCase):
    def setUp(self):
        self.files_to_remove = []

    def tearDown(self):
        for file_name in self.files_to_remove:
            KratosMultiphysics.kratos_utilities.DeleteFileIfExisting(file_name)
        self.files_to_remove = []

    @classmethod
    def _FormatVector(cls, vector, fmt = "{:>.6}"):
        """Mimics KratosMultiphysics.vector.__str__ but with a customizable format string."""
        output = "[{}]".format(vector.Size())
        output += "("
        output += ", ".join([fmt.format(x) for x in vector])
        output += ")"
        return output

    @classmethod
    def _GetGeneratorSettings(cls, geometry):
        """Returns the Kratos Parameters for the symbolic generator."""
        geometry_name = {
            "2D3N": "triangle",
            "2D4N": "quadrilateral",
            "3D4N": "tetrahedron"
        }[geometry]

        parameters = KratosMultiphysics.Parameters("""
        {
            "geometry": "PLEASE SPECIFY A GEOMETRY",
            "template_filename" : "PLEASE SPECIFY A TEMPLATE FILE",
            "output_filename"   : "PLEASE SPECIFY AN OUTPUT FILE",
            "language" : "python"
        }""")

        parameters["geometry"].SetString(geometry_name)
        parameters["template_filename"].SetString("test_{}.py_template".format(geometry))
        parameters["output_filename"].SetString("symbolic_test_{}.py".format(geometry))

        return parameters

    @classmethod
    def _Generate(cls, geometry):
        """Generates the code to be tested, and stores it in a file."""
        params = cls._GetGeneratorSettings(geometry)
        generator = CompressibleNavierStokesSymbolicGenerator(params)
        generator.Generate()
        generator.Write()

        return generator.output_filename

    @classmethod
    def _ImportSubTestSuite(cls, generated_file_name):
        """Imports the generated code as a sub-testsuite."""
        module_name = "compressible_symbolic_generation" + "." + generated_file_name[:-3]

        try:
            test_module = __import__(module_name, fromlist=module_name)
        except ModuleNotFoundError as err:
            raise RuntimeError("Failed to import generated file:", generated_file_name) from err
        return test_module.SubTestSuite()


    def _RunSubTestSuite(self, geometry, sub_testsuite, print_results):
        """Runs the generated code and compares it to the reference."""
        tests = [subtest_name for subtest_name in dir(sub_testsuite)
            if callable(getattr(sub_testsuite, subtest_name)) and subtest_name.startswith("test_")]

        for subtest_name in tests:
            with self.subTest(subtest_name):
                result, reference = getattr(sub_testsuite, subtest_name)()

                if print_results:
                    print("Result for {} -- {}: {}".format(geometry, subtest_name, self._FormatVector(result, "{:>10}")))

                self.assertVectorAlmostEqual(result, reference)

    def _RunTest(self, geometry, print_results=False, cleanup=True):
        """
        Runs the test.

        - geometry -- Choice of geometry. Format is xDyN, with x,y integers
        - print_results -- Set to `True` to print all results to console
        - cleanup -- Set to `True` in order to remove the generated code files
        """
        with KratosUnitTest.WorkFolderScope("compressible_symbolic_generation", __file__):
            generated_file = self._Generate(geometry)

            if cleanup:
                self.files_to_remove.append(os.path.abspath(generated_file))

            sub_testsuite = self._ImportSubTestSuite(generated_file)
            self._RunSubTestSuite(geometry, sub_testsuite, print_results)

    def testSymbolicTriangle(self):
        self._RunTest("2D3N")

    def testSymbolicQuadrilateral(self):
        self._RunTest("2D4N")

    def testSymbolicTetrahedron(self):
        self._RunTest("3D4N")


class CompressibleNavierStokesSymbolicGeneratorUnitTest(KratosUnitTest.TestCase):
    def _assertSympyMatrixEqual(self, first, second, msg=None):
        """Asserts that two sympy matrices are equivalent."""
        class LazyMsg:
            def __init__(self, fmt, *args):
                self.fmt = fmt
                self.args = args

            def __str__(self):
                if self.args:
                    return self.fmt.format(*self.args)
                return self.fmt

        if not msg:
            msg = ""

        self.assertEqual(first.shape, second.shape, msg=LazyMsg("\nMismatching dimensions. {}x{} != {}x{} {}",
                first.shape[0], first.shape[1], second.shape[0], second.shape[1], msg))

        for i in range(first.shape[0]):
            for j in range(first.shape[1]):
                self.assertTrue(first[i,j].equals(second[i,j]),
                    msg=LazyMsg("\nMissmatching entry in position [{},{}]:\n >  first: {}\n > second: {}\n{}",
                        i, j, first[i,j], second[i,j], msg))

    class _DummyGeometry:
        ndims = 2
        nnodes = 3
        blocksize = ndims+2

    def testComputeEulerJacobianMatrix(self):
        defs.SetFormat("python")

        g = self._DummyGeometry()
        U = defs.Vector("U", g.ndims+2)
        DU = defs.Matrix("DU", g.ndims+2, g.ndims)
        params = FormulationParameters(g, "python")
        primitives = PrimitiveMagnitudes(g, params, U, DU)
        A = generate_convective_flux.ComputeEulerJacobianMatrix(U, params, primitives)

        A0_expected = sympy.Matrix([
            [0, 1, 0, 0],
            [(-U[1]**2 + 0.5*(U[1]**2 + U[2]**2)*(params.gamma - 1))/U[0]**2, U[1]*(3.0 - 1.0*params.gamma)/U[0], 1.0*U[2]*(1 - params.gamma)/U[0], params.gamma - 1],
            [-U[1]*U[2]/U[0]**2, U[2]/U[0], U[1]/U[0], 0],
            [1.0*U[1]*(-U[0]*U[3]*params.gamma + U[1]**2*params.gamma - U[1]**2 + U[2]**2*params.gamma - U[2]**2)/U[0]**3, (U[0]*U[3] + 1.0*U[1]**2*(1 - params.gamma) - (params.gamma - 1)*(-U[0]*U[3] + 0.5*U[1]**2 + 0.5*U[2]**2))/U[0]**2, 1.0*U[1]*U[2]*(1 - params.gamma)/U[0]**2, U[1]*params.gamma/U[0]]
        ])

        self._assertSympyMatrixEqual(A[0],  A0_expected)

    def testComputeDiffusiveFlux(self):
        defs.SetFormat("python")

        g = self._DummyGeometry()
        params = FormulationParameters(g, "python")
        U = defs.Vector("U", g.ndims+2)
        DU = defs.Matrix("DU", g.ndims+2, g.ndims)
        primitives = PrimitiveMagnitudes(g, params, U, DU)
        G = generate_diffusive_flux.ComputeDiffusiveFlux(primitives, params)

        G_expected = sympy.Matrix([
            [0, 0],
            [params.mu*(8/3*DU[0,0]*U[1] + 2/3*DU[0,1]*U[2] - 8/3*DU[1,0]*U[0] - 2/3*DU[2,1]*U[0])/U[0]**2, params.mu*(DU[0,0]*U[2] + DU[0,1]*U[1] - DU[1,1]*U[0] - DU[2,0]*U[0])/U[0]**2],
            [params.mu*(DU[0,0]*U[2] + DU[0,1]*U[1] - DU[1,1]*U[0] - DU[2,0]*U[0])/U[0]**2, params.mu*(2/3*DU[0,0]*U[1] + 8/3*DU[0,1]*U[2] - 2/3*DU[1,0]*U[0] - 8/3*DU[2,1]*U[0])/U[0]**2],
            [(params.c_v*params.mu*(U[1]*(8/3*DU[0,0]*U[1] + 2/3*DU[0,1]*U[2] - 8/3*DU[1,0]*U[0] - 2/3*DU[2,1]*U[0]) + U[2]*(DU[0,0]*U[2] + DU[0,1]*U[1] - DU[1,1]*U[0] - DU[2,0]*U[0])) + 1.0*params.lamb*(DU[0,0]*U[0]*U[3] - DU[0,0]*U[1]**2 - DU[0,0]*U[2]**2 + DU[1,0]*U[0]*U[1] + DU[2,0]*U[0]*U[2] - DU[3,0]*U[0]**2))/(U[0]**3*params.c_v), (params.c_v*params.mu*(U[1]*(DU[0,0]*U[2] + DU[0,1]*U[1] - DU[1,1]*U[0] - DU[2,0]*U[0]) + U[2]*(2/3*DU[0,0]*U[1] + 8/3*DU[0,1]*U[2] - 2/3*DU[1,0]*U[0] - 8/3*DU[2,1]*U[0])) + 1.0*params.lamb*(DU[0,1]*U[0]*U[3] - DU[0,1]*U[1]**2 - DU[0,1]*U[2]**2 + DU[1,1]*U[0]*U[1] + DU[2,1]*U[0]*U[2] - DU[3,1]*U[0]**2))/(U[0]**3*params.c_v)]
        ])

        self._assertSympyMatrixEqual(G,  G_expected)

    def testComputeDiffusiveFluxWithShockCapturing(self):
        defs.SetFormat("python")

        g = self._DummyGeometry()
        params = FormulationParameters(g, "python")
        sc_params = ShockCapturingParameters()
        U = defs.Vector("U", g.ndims+2)
        DU = defs.Matrix("DU", g.ndims+2, g.ndims)
        primitives = PrimitiveMagnitudes(g, params, U, DU)
        G = generate_diffusive_flux.ComputeDiffusiveFluxWithShockCapturing(DU, primitives, params, sc_params)

        G_expected = sympy.Matrix([
            [DU[0,0]*sc_params.alpha, DU[0,1]*sc_params.alpha],
            [(2*(params.mu + sc_params.mu)*(DU[0,0]*U[1] - DU[1,0]*U[0]) + (-sc_params.beta + 2/3*params.mu + 2/3*sc_params.mu)*(DU[0,0]*U[1] + DU[0,1]*U[2] - DU[1,0]*U[0] - DU[2,1]*U[0]))/U[0]**2, (params.mu + sc_params.mu)*(DU[0,0]*U[2] + DU[0,1]*U[1] - DU[1,1]*U[0] - DU[2,0]*U[0])/U[0]**2],
            [(params.mu + sc_params.mu)*(DU[0,0]*U[2] + DU[0,1]*U[1] - DU[1,1]*U[0] - DU[2,0]*U[0])/U[0]**2, (2*(params.mu + sc_params.mu)*(DU[0,1]*U[2] - DU[2,1]*U[0]) + (-sc_params.beta + 2/3*params.mu + 2/3*sc_params.mu)*(DU[0,0]*U[1] + DU[0,1]*U[2] - DU[1,0]*U[0] - DU[2,1]*U[0]))/U[0]**2],
            [(params.c_v*(U[1]*(2*(params.mu + sc_params.mu)*(DU[0,0]*U[1] - DU[1,0]*U[0]) + (-sc_params.beta + 2/3*params.mu + 2/3*sc_params.mu)*(DU[0,0]*U[1] + DU[0,1]*U[2] - DU[1,0]*U[0] - DU[2,1]*U[0])) + U[2]*(params.mu + sc_params.mu)*(DU[0,0]*U[2] + DU[0,1]*U[1] - DU[1,1]*U[0] - DU[2,0]*U[0])) + 1.0*(sc_params.lamb + params.lamb)*(DU[0,0]*U[0]*U[3] - DU[0,0]*U[1]**2 - DU[0,0]*U[2]**2 + DU[1,0]*U[0]*U[1] + DU[2,0]*U[0]*U[2] - DU[3,0]*U[0]**2))/(U[0]**3*params.c_v), (params.c_v*(U[1]*(params.mu + sc_params.mu)*(DU[0,0]*U[2] + DU[0,1]*U[1] - DU[1,1]*U[0] - DU[2,0]*U[0]) + U[2]*(2*(params.mu + sc_params.mu)*(DU[0,1]*U[2] - DU[2,1]*U[0]) + (-sc_params.beta + 2/3*params.mu + 2/3*sc_params.mu)*(DU[0,0]*U[1] + DU[0,1]*U[2] - DU[1,0]*U[0] - DU[2,1]*U[0]))) + 1.0*(sc_params.lamb + params.lamb)*(DU[0,1]*U[0]*U[3] - DU[0,1]*U[1]**2 - DU[0,1]*U[2]**2 + DU[1,1]*U[0]*U[1] + DU[2,1]*U[0]*U[2] - DU[3,1]*U[0]**2))/(U[0]**3*params.c_v)]
        ])

        G.simplify()
        G_expected.simplify()

        self._assertSympyMatrixEqual(G,  G_expected)

    def testComputeStabilizationMatrix(self):
        g = self._DummyGeometry()
        params = FormulationParameters(g, "python")
        tau = generate_stabilization_matrix.ComputeStabilizationMatrix(params)

        tau_expected = sympy.Matrix([
            [sympy.Symbol('tau1'),                    0,                    0,                    0],
            [                   0, sympy.Symbol('tau2'),                    0,                    0],
            [                   0,                    0, sympy.Symbol('tau2'),                    0],
            [                   0,                    0,                    0, sympy.Symbol('tau3')]
        ])

        self._assertSympyMatrixEqual(tau,  tau_expected)

    class _DummyGenerator(CompressibleNavierStokesSymbolicGenerator):
        def __init__(self, geometry_class):
            self.is_explicit = True
            self.echo_level = 0
            self.geometry = geometry_class
            # self.write_language = "python"
            # self.shock_capturing = True
            # self.primitive_interpolation = "nodal"
            # self.outstring = None

        def ComputeNonLinearOperator(self, A, H, S, Ug):
            return super()._ComputeNonLinearOperator(A, H, S, Ug)

        def ComputeNonLinearAdjointOperator(self, A, H, Q, S, Ug, V):
            return super()._ComputeNonLinearAdjointOperator(A, H, Q, S, Ug, V)

        def ComputeVariationalFormulation(self, A, acc, G, H, L_adj, Q, S, Ug, V):
            return super()._ComputeVariationalFormulation(A, acc, G, H, L_adj, Q, S, Ug, V)

        def ComputeOSSProjectionsAtGaussPoint(self, acc, bdf, dUdt, f, forcing_terms, H, i_gauss, mg, params, primitives, projections, res, rg, U, Ug, Un, Unn):
            return super()._ComputeOSSProjectionsAtGaussPoint(acc, bdf, dUdt, f, forcing_terms, H, i_gauss, mg, params, primitives, projections, res, rg, U, Ug, Un, Unn)


    def testComputeNonLinearOperator(self):
        dim = 2
        blocksize = dim+2
        A = [defs.Matrix('A[{}]'.format(d), blocksize, blocksize) for d in range(dim)]
        H = defs.Matrix('H', blocksize, dim)
        S = defs.Matrix('S', blocksize, blocksize)
        U = defs.Vector('U', blocksize)

        dummy_geneator = self._DummyGenerator(self._DummyGeometry)
        L = dummy_geneator.ComputeNonLinearOperator(A, H, S, U)

        L_expected = sympy.Matrix([
            [A[0][0,0]*H[0,0] + A[0][0,1]*H[1,0] + A[0][0,2]*H[2,0] + A[0][0,3]*H[3,0] + A[1][0,0]*H[0,1] + A[1][0,1]*H[1,1] + A[1][0,2]*H[2,1] + A[1][0,3]*H[3,1] - S[0,0]*U[0] - S[0,1]*U[1] - S[0,2]*U[2] - S[0,3]*U[3]],
            [A[0][1,0]*H[0,0] + A[0][1,1]*H[1,0] + A[0][1,2]*H[2,0] + A[0][1,3]*H[3,0] + A[1][1,0]*H[0,1] + A[1][1,1]*H[1,1] + A[1][1,2]*H[2,1] + A[1][1,3]*H[3,1] - S[1,0]*U[0] - S[1,1]*U[1] - S[1,2]*U[2] - S[1,3]*U[3]],
            [A[0][2,0]*H[0,0] + A[0][2,1]*H[1,0] + A[0][2,2]*H[2,0] + A[0][2,3]*H[3,0] + A[1][2,0]*H[0,1] + A[1][2,1]*H[1,1] + A[1][2,2]*H[2,1] + A[1][2,3]*H[3,1] - S[2,0]*U[0] - S[2,1]*U[1] - S[2,2]*U[2] - S[2,3]*U[3]],
            [A[0][3,0]*H[0,0] + A[0][3,1]*H[1,0] + A[0][3,2]*H[2,0] + A[0][3,3]*H[3,0] + A[1][3,0]*H[0,1] + A[1][3,1]*H[1,1] + A[1][3,2]*H[2,1] + A[1][3,3]*H[3,1] - S[3,0]*U[0] - S[3,1]*U[1] - S[3,2]*U[2] - S[3,3]*U[3]]
        ])

        self._assertSympyMatrixEqual(L, L_expected)

    def testComputeNonLinearAdjointOperator(self):
        dim = 2
        blocksize = dim+2
        A = [defs.Matrix('A[{}]'.format(d), blocksize, blocksize) for d in range(dim)]
        H = defs.Matrix('H', blocksize, dim)
        Q = defs.Matrix('Q', blocksize, dim)
        S = defs.Matrix('S', blocksize, blocksize)
        U = defs.Vector('U', blocksize)
        V = defs.Vector('V', blocksize)

        dummy_geneator = self._DummyGenerator(self._DummyGeometry)
        Ladj = dummy_geneator.ComputeNonLinearAdjointOperator(A, H, Q, S, U, V)

        Ladj_expected = sympy.Matrix([
            [A[0][0,0]*Q[0,0] + A[0][1,0]*Q[1,0] + A[0][2,0]*Q[2,0] + A[0][3,0]*Q[3,0] + A[1][0,0]*Q[0,1] + A[1][1,0]*Q[1,1] + A[1][2,0]*Q[2,1] + A[1][3,0]*Q[3,1] + S[0,0]*V[0] + S[1,0]*V[1] + S[2,0]*V[2] + S[3,0]*V[3]],
            [A[0][0,1]*Q[0,0] + A[0][1,1]*Q[1,0] + A[0][2,1]*Q[2,0] + A[0][3,1]*Q[3,0] + A[1][0,1]*Q[0,1] + A[1][1,1]*Q[1,1] + A[1][2,1]*Q[2,1] + A[1][3,1]*Q[3,1] + S[0,1]*V[0] + S[1,1]*V[1] + S[2,1]*V[2] + S[3,1]*V[3]],
            [A[0][0,2]*Q[0,0] + A[0][1,2]*Q[1,0] + A[0][2,2]*Q[2,0] + A[0][3,2]*Q[3,0] + A[1][0,2]*Q[0,1] + A[1][1,2]*Q[1,1] + A[1][2,2]*Q[2,1] + A[1][3,2]*Q[3,1] + S[0,2]*V[0] + S[1,2]*V[1] + S[2,2]*V[2] + S[3,2]*V[3]],
            [A[0][0,3]*Q[0,0] + A[0][1,3]*Q[1,0] + A[0][2,3]*Q[2,0] + A[0][3,3]*Q[3,0] + A[1][0,3]*Q[0,1] + A[1][1,3]*Q[1,1] + A[1][2,3]*Q[2,1] + A[1][3,3]*Q[3,1] + S[0,3]*V[0] + S[1,3]*V[1] + S[2,3]*V[2] + S[3,3]*V[3]]
        ])

        self._assertSympyMatrixEqual(Ladj, Ladj_expected)

    def testComputeVariationalFormulation(self):
        dim = 2
        blocksize = dim+2
        A = [defs.Matrix('A[{}]'.format(d), blocksize, blocksize) for d in range(dim)]
        H = defs.Matrix('H', blocksize, dim)
        G = defs.Matrix('G', blocksize, dim)
        Ladj = defs.Vector('Ladj', blocksize)
        Q = defs.Matrix('Q', blocksize, dim)
        S = defs.Matrix('S', blocksize, blocksize)
        U = defs.Vector('U', blocksize)
        V = defs.Vector('V', blocksize)
        acc = defs.Vector('acc', blocksize)

        dummy_geneator = self._DummyGenerator(self._DummyGeometry)
        rv, subscales = dummy_geneator.ComputeVariationalFormulation(A, acc, G, H, Ladj, Q, S, U, V)

        rv_expected = sympy.Matrix([
            G[0,0]*Q[0,0] + G[0,1]*Q[0,1] + G[1,0]*Q[1,0] + G[1,1]*Q[1,1] + G[2,0]*Q[2,0] + G[2,1]*Q[2,1] + G[3,0]*Q[3,0] + G[3,1]*Q[3,1] + Ladj[0]*subscales[0] + Ladj[1]*subscales[1] + Ladj[2]*subscales[2] + Ladj[3]*subscales[3] + V[0]*(S[0,0]*U[0] + S[0,1]*U[1] + S[0,2]*U[2] + S[0,3]*U[3]) - V[0]*(A[0][0,0]*H[0,0] + A[0][0,1]*H[1,0] + A[0][0,2]*H[2,0] + A[0][0,3]*H[3,0] + A[1][0,0]*H[0,1] + A[1][0,1]*H[1,1] + A[1][0,2]*H[2,1] + A[1][0,3]*H[3,1]) + V[1]*(S[1,0]*U[0] + S[1,1]*U[1] + S[1,2]*U[2] + S[1,3]*U[3]) - V[1]*(A[0][1,0]*H[0,0] + A[0][1,1]*H[1,0] + A[0][1,2]*H[2,0] + A[0][1,3]*H[3,0] + A[1][1,0]*H[0,1] + A[1][1,1]*H[1,1] + A[1][1,2]*H[2,1] + A[1][1,3]*H[3,1]) + V[2]*(S[2,0]*U[0] + S[2,1]*U[1] + S[2,2]*U[2] + S[2,3]*U[3]) - V[2]*(A[0][2,0]*H[0,0] + A[0][2,1]*H[1,0] + A[0][2,2]*H[2,0] + A[0][2,3]*H[3,0] + A[1][2,0]*H[0,1] + A[1][2,1]*H[1,1] + A[1][2,2]*H[2,1] + A[1][2,3]*H[3,1]) + V[3]*(S[3,0]*U[0] + S[3,1]*U[1] + S[3,2]*U[2] + S[3,3]*U[3]) - V[3]*(A[0][3,0]*H[0,0] + A[0][3,1]*H[1,0] + A[0][3,2]*H[2,0] + A[0][3,3]*H[3,0] + A[1][3,0]*H[0,1] + A[1][3,1]*H[1,1] + A[1][3,2]*H[2,1] + A[1][3,3]*H[3,1])
        ])

        subscales_expected = defs.Vector('subscales', blocksize)

        self._assertSympyMatrixEqual(rv, rv_expected)
        self._assertSympyMatrixEqual(subscales, subscales_expected)

    def testComputeOSSProjectionsAtGaussPoint(self):
        geometry = GeometryDataFactory("quadrilateral")
        dim = geometry.ndims
        blocksize = geometry.blocksize
        numnodes = geometry.nnodes

        params = FormulationParameters(geometry, "python")
        H = defs.Matrix('H', blocksize, dim)
        Ug = defs.Vector('Ug', blocksize)
        primitives = PrimitiveMagnitudes(geometry, params, Ug, H)
        U = defs.Matrix('U', numnodes, blocksize)
        mg  = sympy.Symbol('mg')
        f = defs.Vector('f', dim)
        acc = defs.Vector('acc', blocksize)
        rg  = sympy.Symbol('rg')

        bdf = None
        dUdt = defs.Matrix('dUdt', numnodes, blocksize)
        Un = None
        Unn = None

        forcing_terms = {
            "mass":    defs.Vector('forcing_terms["mass"]', numnodes),
            "thermal": defs.Vector('forcing_terms["thermal"]', numnodes),
            "force":   defs.Matrix('forcing_terms["force"]', numnodes, dim)
        }

        projections = {
            "rho"      : defs.ZeroVector(numnodes),
            "momentum" : defs.ZeroVector(numnodes*dim),
            "energy"   : defs.ZeroVector(numnodes)
        }

        # Just an example nonsense set of expressions
        inner_prod = lambda U,V: sum([u*v for u,v in zip(U,V)])
        res = sympy.Matrix([
            acc[0] + inner_prod(primitives.V, H[0,:]) - mg,
            Ug[0]*sympy.Matrix(acc[1:-1]) + Ug[0]*primitives.grad_V*primitives.V - f,
            acc[-1] + inner_prod(primitives.V, primitives.grad_T) - rg - inner_prod(primitives.V, f)
        ])

        dummy_geneator = self._DummyGenerator(geometry)
        dummy_geneator.ComputeOSSProjectionsAtGaussPoint(acc, bdf, dUdt, f, forcing_terms, H, 0, mg, params, primitives, projections, res, rg, U, Ug, Un, Unn)

        N = geometry.N_gauss(0)
        DN_DX = geometry.DN()

        expected_rho_proj = sympy.Matrix([
            [N[0]*(-N[0]*forcing_terms["mass"][0] + N[0]*dUdt[0,0] - N[1]*forcing_terms["mass"][1] + N[1]*dUdt[1,0] - N[2]*forcing_terms["mass"][2] + N[2]*dUdt[2,0] - N[3]*forcing_terms["mass"][3] + N[3]*dUdt[3,0] + (DN_DX[0,0]*U[0,0] + DN_DX[1,0]*U[1,0] + DN_DX[2,0]*U[2,0] + DN_DX[3,0]*U[3,0])*(N[0]*U[0,1] + N[1]*U[1,1] + N[2]*U[2,1] + N[3]*U[3,1])/(N[0]*U[0,0] + N[1]*U[1,0] + N[2]*U[2,0] + N[3]*U[3,0]) + (DN_DX[0,1]*U[0,0] + DN_DX[1,1]*U[1,0] + DN_DX[2,1]*U[2,0] + DN_DX[3,1]*U[3,0])*(N[0]*U[0,2] + N[1]*U[1,2] + N[2]*U[2,2] + N[3]*U[3,2])/(N[0]*U[0,0] + N[1]*U[1,0] + N[2]*U[2,0] + N[3]*U[3,0]))],
            [N[1]*(-N[0]*forcing_terms["mass"][0] + N[0]*dUdt[0,0] - N[1]*forcing_terms["mass"][1] + N[1]*dUdt[1,0] - N[2]*forcing_terms["mass"][2] + N[2]*dUdt[2,0] - N[3]*forcing_terms["mass"][3] + N[3]*dUdt[3,0] + (DN_DX[0,0]*U[0,0] + DN_DX[1,0]*U[1,0] + DN_DX[2,0]*U[2,0] + DN_DX[3,0]*U[3,0])*(N[0]*U[0,1] + N[1]*U[1,1] + N[2]*U[2,1] + N[3]*U[3,1])/(N[0]*U[0,0] + N[1]*U[1,0] + N[2]*U[2,0] + N[3]*U[3,0]) + (DN_DX[0,1]*U[0,0] + DN_DX[1,1]*U[1,0] + DN_DX[2,1]*U[2,0] + DN_DX[3,1]*U[3,0])*(N[0]*U[0,2] + N[1]*U[1,2] + N[2]*U[2,2] + N[3]*U[3,2])/(N[0]*U[0,0] + N[1]*U[1,0] + N[2]*U[2,0] + N[3]*U[3,0]))],
            [N[2]*(-N[0]*forcing_terms["mass"][0] + N[0]*dUdt[0,0] - N[1]*forcing_terms["mass"][1] + N[1]*dUdt[1,0] - N[2]*forcing_terms["mass"][2] + N[2]*dUdt[2,0] - N[3]*forcing_terms["mass"][3] + N[3]*dUdt[3,0] + (DN_DX[0,0]*U[0,0] + DN_DX[1,0]*U[1,0] + DN_DX[2,0]*U[2,0] + DN_DX[3,0]*U[3,0])*(N[0]*U[0,1] + N[1]*U[1,1] + N[2]*U[2,1] + N[3]*U[3,1])/(N[0]*U[0,0] + N[1]*U[1,0] + N[2]*U[2,0] + N[3]*U[3,0]) + (DN_DX[0,1]*U[0,0] + DN_DX[1,1]*U[1,0] + DN_DX[2,1]*U[2,0] + DN_DX[3,1]*U[3,0])*(N[0]*U[0,2] + N[1]*U[1,2] + N[2]*U[2,2] + N[3]*U[3,2])/(N[0]*U[0,0] + N[1]*U[1,0] + N[2]*U[2,0] + N[3]*U[3,0]))],
            [N[3]*(-N[0]*forcing_terms["mass"][0] + N[0]*dUdt[0,0] - N[1]*forcing_terms["mass"][1] + N[1]*dUdt[1,0] - N[2]*forcing_terms["mass"][2] + N[2]*dUdt[2,0] - N[3]*forcing_terms["mass"][3] + N[3]*dUdt[3,0] + (DN_DX[0,0]*U[0,0] + DN_DX[1,0]*U[1,0] + DN_DX[2,0]*U[2,0] + DN_DX[3,0]*U[3,0])*(N[0]*U[0,1] + N[1]*U[1,1] + N[2]*U[2,1] + N[3]*U[3,1])/(N[0]*U[0,0] + N[1]*U[1,0] + N[2]*U[2,0] + N[3]*U[3,0]) + (DN_DX[0,1]*U[0,0] + DN_DX[1,1]*U[1,0] + DN_DX[2,1]*U[2,0] + DN_DX[3,1]*U[3,0])*(N[0]*U[0,2] + N[1]*U[1,2] + N[2]*U[2,2] + N[3]*U[3,2])/(N[0]*U[0,0] + N[1]*U[1,0] + N[2]*U[2,0] + N[3]*U[3,0]))]
        ])

        expected_mom_proj = sympy.Matrix([
            [N[0]*(-N[0]*forcing_terms["force"][0,0] - N[1]*forcing_terms["force"][1,0] - N[2]*forcing_terms["force"][2,0] - N[3]*forcing_terms["force"][3,0] + (-(DN_DX[0,0]*U[0,0] + DN_DX[1,0]*U[1,0] + DN_DX[2,0]*U[2,0] + DN_DX[3,0]*U[3,0])*(N[0]*U[0,1] + N[1]*U[1,1] + N[2]*U[2,1] + N[3]*U[3,1]) + (DN_DX[0,0]*U[0,1] + DN_DX[1,0]*U[1,1] + DN_DX[2,0]*U[2,1] + DN_DX[3,0]*U[3,1])*(N[0]*U[0,0] + N[1]*U[1,0] + N[2]*U[2,0] + N[3]*U[3,0]))*(N[0]*U[0,1] + N[1]*U[1,1] + N[2]*U[2,1] + N[3]*U[3,1])/(N[0]*U[0,0] + N[1]*U[1,0] + N[2]*U[2,0] + N[3]*U[3,0])**2 + (-(DN_DX[0,0]*U[0,0] + DN_DX[1,0]*U[1,0] + DN_DX[2,0]*U[2,0] + DN_DX[3,0]*U[3,0])*(N[0]*U[0,2] + N[1]*U[1,2] + N[2]*U[2,2] + N[3]*U[3,2]) + (DN_DX[0,0]*U[0,2] + DN_DX[1,0]*U[1,2] + DN_DX[2,0]*U[2,2] + DN_DX[3,0]*U[3,2])*(N[0]*U[0,0] + N[1]*U[1,0] + N[2]*U[2,0] + N[3]*U[3,0]))*(N[0]*U[0,2] + N[1]*U[1,2] + N[2]*U[2,2] + N[3]*U[3,2])/(N[0]*U[0,0] + N[1]*U[1,0] + N[2]*U[2,0] + N[3]*U[3,0])**2 + (N[0]*U[0,0] + N[1]*U[1,0] + N[2]*U[2,0] + N[3]*U[3,0])*(N[0]*dUdt[0,1] + N[1]*dUdt[1,1] + N[2]*dUdt[2,1] + N[3]*dUdt[3,1]))],
            [N[0]*(-N[0]*forcing_terms["force"][0,1] - N[1]*forcing_terms["force"][1,1] - N[2]*forcing_terms["force"][2,1] - N[3]*forcing_terms["force"][3,1] + (-(DN_DX[0,1]*U[0,0] + DN_DX[1,1]*U[1,0] + DN_DX[2,1]*U[2,0] + DN_DX[3,1]*U[3,0])*(N[0]*U[0,1] + N[1]*U[1,1] + N[2]*U[2,1] + N[3]*U[3,1]) + (DN_DX[0,1]*U[0,1] + DN_DX[1,1]*U[1,1] + DN_DX[2,1]*U[2,1] + DN_DX[3,1]*U[3,1])*(N[0]*U[0,0] + N[1]*U[1,0] + N[2]*U[2,0] + N[3]*U[3,0]))*(N[0]*U[0,1] + N[1]*U[1,1] + N[2]*U[2,1] + N[3]*U[3,1])/(N[0]*U[0,0] + N[1]*U[1,0] + N[2]*U[2,0] + N[3]*U[3,0])**2 + (-(DN_DX[0,1]*U[0,0] + DN_DX[1,1]*U[1,0] + DN_DX[2,1]*U[2,0] + DN_DX[3,1]*U[3,0])*(N[0]*U[0,2] + N[1]*U[1,2] + N[2]*U[2,2] + N[3]*U[3,2]) + (DN_DX[0,1]*U[0,2] + DN_DX[1,1]*U[1,2] + DN_DX[2,1]*U[2,2] + DN_DX[3,1]*U[3,2])*(N[0]*U[0,0] + N[1]*U[1,0] + N[2]*U[2,0] + N[3]*U[3,0]))*(N[0]*U[0,2] + N[1]*U[1,2] + N[2]*U[2,2] + N[3]*U[3,2])/(N[0]*U[0,0] + N[1]*U[1,0] + N[2]*U[2,0] + N[3]*U[3,0])**2 + (N[0]*U[0,0] + N[1]*U[1,0] + N[2]*U[2,0] + N[3]*U[3,0])*(N[0]*dUdt[0,2] + N[1]*dUdt[1,2] + N[2]*dUdt[2,2] + N[3]*dUdt[3,2]))],
            [N[1]*(-N[0]*forcing_terms["force"][0,0] - N[1]*forcing_terms["force"][1,0] - N[2]*forcing_terms["force"][2,0] - N[3]*forcing_terms["force"][3,0] + (-(DN_DX[0,0]*U[0,0] + DN_DX[1,0]*U[1,0] + DN_DX[2,0]*U[2,0] + DN_DX[3,0]*U[3,0])*(N[0]*U[0,1] + N[1]*U[1,1] + N[2]*U[2,1] + N[3]*U[3,1]) + (DN_DX[0,0]*U[0,1] + DN_DX[1,0]*U[1,1] + DN_DX[2,0]*U[2,1] + DN_DX[3,0]*U[3,1])*(N[0]*U[0,0] + N[1]*U[1,0] + N[2]*U[2,0] + N[3]*U[3,0]))*(N[0]*U[0,1] + N[1]*U[1,1] + N[2]*U[2,1] + N[3]*U[3,1])/(N[0]*U[0,0] + N[1]*U[1,0] + N[2]*U[2,0] + N[3]*U[3,0])**2 + (-(DN_DX[0,0]*U[0,0] + DN_DX[1,0]*U[1,0] + DN_DX[2,0]*U[2,0] + DN_DX[3,0]*U[3,0])*(N[0]*U[0,2] + N[1]*U[1,2] + N[2]*U[2,2] + N[3]*U[3,2]) + (DN_DX[0,0]*U[0,2] + DN_DX[1,0]*U[1,2] + DN_DX[2,0]*U[2,2] + DN_DX[3,0]*U[3,2])*(N[0]*U[0,0] + N[1]*U[1,0] + N[2]*U[2,0] + N[3]*U[3,0]))*(N[0]*U[0,2] + N[1]*U[1,2] + N[2]*U[2,2] + N[3]*U[3,2])/(N[0]*U[0,0] + N[1]*U[1,0] + N[2]*U[2,0] + N[3]*U[3,0])**2 + (N[0]*U[0,0] + N[1]*U[1,0] + N[2]*U[2,0] + N[3]*U[3,0])*(N[0]*dUdt[0,1] + N[1]*dUdt[1,1] + N[2]*dUdt[2,1] + N[3]*dUdt[3,1]))],
            [N[1]*(-N[0]*forcing_terms["force"][0,1] - N[1]*forcing_terms["force"][1,1] - N[2]*forcing_terms["force"][2,1] - N[3]*forcing_terms["force"][3,1] + (-(DN_DX[0,1]*U[0,0] + DN_DX[1,1]*U[1,0] + DN_DX[2,1]*U[2,0] + DN_DX[3,1]*U[3,0])*(N[0]*U[0,1] + N[1]*U[1,1] + N[2]*U[2,1] + N[3]*U[3,1]) + (DN_DX[0,1]*U[0,1] + DN_DX[1,1]*U[1,1] + DN_DX[2,1]*U[2,1] + DN_DX[3,1]*U[3,1])*(N[0]*U[0,0] + N[1]*U[1,0] + N[2]*U[2,0] + N[3]*U[3,0]))*(N[0]*U[0,1] + N[1]*U[1,1] + N[2]*U[2,1] + N[3]*U[3,1])/(N[0]*U[0,0] + N[1]*U[1,0] + N[2]*U[2,0] + N[3]*U[3,0])**2 + (-(DN_DX[0,1]*U[0,0] + DN_DX[1,1]*U[1,0] + DN_DX[2,1]*U[2,0] + DN_DX[3,1]*U[3,0])*(N[0]*U[0,2] + N[1]*U[1,2] + N[2]*U[2,2] + N[3]*U[3,2]) + (DN_DX[0,1]*U[0,2] + DN_DX[1,1]*U[1,2] + DN_DX[2,1]*U[2,2] + DN_DX[3,1]*U[3,2])*(N[0]*U[0,0] + N[1]*U[1,0] + N[2]*U[2,0] + N[3]*U[3,0]))*(N[0]*U[0,2] + N[1]*U[1,2] + N[2]*U[2,2] + N[3]*U[3,2])/(N[0]*U[0,0] + N[1]*U[1,0] + N[2]*U[2,0] + N[3]*U[3,0])**2 + (N[0]*U[0,0] + N[1]*U[1,0] + N[2]*U[2,0] + N[3]*U[3,0])*(N[0]*dUdt[0,2] + N[1]*dUdt[1,2] + N[2]*dUdt[2,2] + N[3]*dUdt[3,2]))],
            [N[2]*(-N[0]*forcing_terms["force"][0,0] - N[1]*forcing_terms["force"][1,0] - N[2]*forcing_terms["force"][2,0] - N[3]*forcing_terms["force"][3,0] + (-(DN_DX[0,0]*U[0,0] + DN_DX[1,0]*U[1,0] + DN_DX[2,0]*U[2,0] + DN_DX[3,0]*U[3,0])*(N[0]*U[0,1] + N[1]*U[1,1] + N[2]*U[2,1] + N[3]*U[3,1]) + (DN_DX[0,0]*U[0,1] + DN_DX[1,0]*U[1,1] + DN_DX[2,0]*U[2,1] + DN_DX[3,0]*U[3,1])*(N[0]*U[0,0] + N[1]*U[1,0] + N[2]*U[2,0] + N[3]*U[3,0]))*(N[0]*U[0,1] + N[1]*U[1,1] + N[2]*U[2,1] + N[3]*U[3,1])/(N[0]*U[0,0] + N[1]*U[1,0] + N[2]*U[2,0] + N[3]*U[3,0])**2 + (-(DN_DX[0,0]*U[0,0] + DN_DX[1,0]*U[1,0] + DN_DX[2,0]*U[2,0] + DN_DX[3,0]*U[3,0])*(N[0]*U[0,2] + N[1]*U[1,2] + N[2]*U[2,2] + N[3]*U[3,2]) + (DN_DX[0,0]*U[0,2] + DN_DX[1,0]*U[1,2] + DN_DX[2,0]*U[2,2] + DN_DX[3,0]*U[3,2])*(N[0]*U[0,0] + N[1]*U[1,0] + N[2]*U[2,0] + N[3]*U[3,0]))*(N[0]*U[0,2] + N[1]*U[1,2] + N[2]*U[2,2] + N[3]*U[3,2])/(N[0]*U[0,0] + N[1]*U[1,0] + N[2]*U[2,0] + N[3]*U[3,0])**2 + (N[0]*U[0,0] + N[1]*U[1,0] + N[2]*U[2,0] + N[3]*U[3,0])*(N[0]*dUdt[0,1] + N[1]*dUdt[1,1] + N[2]*dUdt[2,1] + N[3]*dUdt[3,1]))],
            [N[2]*(-N[0]*forcing_terms["force"][0,1] - N[1]*forcing_terms["force"][1,1] - N[2]*forcing_terms["force"][2,1] - N[3]*forcing_terms["force"][3,1] + (-(DN_DX[0,1]*U[0,0] + DN_DX[1,1]*U[1,0] + DN_DX[2,1]*U[2,0] + DN_DX[3,1]*U[3,0])*(N[0]*U[0,1] + N[1]*U[1,1] + N[2]*U[2,1] + N[3]*U[3,1]) + (DN_DX[0,1]*U[0,1] + DN_DX[1,1]*U[1,1] + DN_DX[2,1]*U[2,1] + DN_DX[3,1]*U[3,1])*(N[0]*U[0,0] + N[1]*U[1,0] + N[2]*U[2,0] + N[3]*U[3,0]))*(N[0]*U[0,1] + N[1]*U[1,1] + N[2]*U[2,1] + N[3]*U[3,1])/(N[0]*U[0,0] + N[1]*U[1,0] + N[2]*U[2,0] + N[3]*U[3,0])**2 + (-(DN_DX[0,1]*U[0,0] + DN_DX[1,1]*U[1,0] + DN_DX[2,1]*U[2,0] + DN_DX[3,1]*U[3,0])*(N[0]*U[0,2] + N[1]*U[1,2] + N[2]*U[2,2] + N[3]*U[3,2]) + (DN_DX[0,1]*U[0,2] + DN_DX[1,1]*U[1,2] + DN_DX[2,1]*U[2,2] + DN_DX[3,1]*U[3,2])*(N[0]*U[0,0] + N[1]*U[1,0] + N[2]*U[2,0] + N[3]*U[3,0]))*(N[0]*U[0,2] + N[1]*U[1,2] + N[2]*U[2,2] + N[3]*U[3,2])/(N[0]*U[0,0] + N[1]*U[1,0] + N[2]*U[2,0] + N[3]*U[3,0])**2 + (N[0]*U[0,0] + N[1]*U[1,0] + N[2]*U[2,0] + N[3]*U[3,0])*(N[0]*dUdt[0,2] + N[1]*dUdt[1,2] + N[2]*dUdt[2,2] + N[3]*dUdt[3,2]))],
            [N[3]*(-N[0]*forcing_terms["force"][0,0] - N[1]*forcing_terms["force"][1,0] - N[2]*forcing_terms["force"][2,0] - N[3]*forcing_terms["force"][3,0] + (-(DN_DX[0,0]*U[0,0] + DN_DX[1,0]*U[1,0] + DN_DX[2,0]*U[2,0] + DN_DX[3,0]*U[3,0])*(N[0]*U[0,1] + N[1]*U[1,1] + N[2]*U[2,1] + N[3]*U[3,1]) + (DN_DX[0,0]*U[0,1] + DN_DX[1,0]*U[1,1] + DN_DX[2,0]*U[2,1] + DN_DX[3,0]*U[3,1])*(N[0]*U[0,0] + N[1]*U[1,0] + N[2]*U[2,0] + N[3]*U[3,0]))*(N[0]*U[0,1] + N[1]*U[1,1] + N[2]*U[2,1] + N[3]*U[3,1])/(N[0]*U[0,0] + N[1]*U[1,0] + N[2]*U[2,0] + N[3]*U[3,0])**2 + (-(DN_DX[0,0]*U[0,0] + DN_DX[1,0]*U[1,0] + DN_DX[2,0]*U[2,0] + DN_DX[3,0]*U[3,0])*(N[0]*U[0,2] + N[1]*U[1,2] + N[2]*U[2,2] + N[3]*U[3,2]) + (DN_DX[0,0]*U[0,2] + DN_DX[1,0]*U[1,2] + DN_DX[2,0]*U[2,2] + DN_DX[3,0]*U[3,2])*(N[0]*U[0,0] + N[1]*U[1,0] + N[2]*U[2,0] + N[3]*U[3,0]))*(N[0]*U[0,2] + N[1]*U[1,2] + N[2]*U[2,2] + N[3]*U[3,2])/(N[0]*U[0,0] + N[1]*U[1,0] + N[2]*U[2,0] + N[3]*U[3,0])**2 + (N[0]*U[0,0] + N[1]*U[1,0] + N[2]*U[2,0] + N[3]*U[3,0])*(N[0]*dUdt[0,1] + N[1]*dUdt[1,1] + N[2]*dUdt[2,1] + N[3]*dUdt[3,1]))],
            [N[3]*(-N[0]*forcing_terms["force"][0,1] - N[1]*forcing_terms["force"][1,1] - N[2]*forcing_terms["force"][2,1] - N[3]*forcing_terms["force"][3,1] + (-(DN_DX[0,1]*U[0,0] + DN_DX[1,1]*U[1,0] + DN_DX[2,1]*U[2,0] + DN_DX[3,1]*U[3,0])*(N[0]*U[0,1] + N[1]*U[1,1] + N[2]*U[2,1] + N[3]*U[3,1]) + (DN_DX[0,1]*U[0,1] + DN_DX[1,1]*U[1,1] + DN_DX[2,1]*U[2,1] + DN_DX[3,1]*U[3,1])*(N[0]*U[0,0] + N[1]*U[1,0] + N[2]*U[2,0] + N[3]*U[3,0]))*(N[0]*U[0,1] + N[1]*U[1,1] + N[2]*U[2,1] + N[3]*U[3,1])/(N[0]*U[0,0] + N[1]*U[1,0] + N[2]*U[2,0] + N[3]*U[3,0])**2 + (-(DN_DX[0,1]*U[0,0] + DN_DX[1,1]*U[1,0] + DN_DX[2,1]*U[2,0] + DN_DX[3,1]*U[3,0])*(N[0]*U[0,2] + N[1]*U[1,2] + N[2]*U[2,2] + N[3]*U[3,2]) + (DN_DX[0,1]*U[0,2] + DN_DX[1,1]*U[1,2] + DN_DX[2,1]*U[2,2] + DN_DX[3,1]*U[3,2])*(N[0]*U[0,0] + N[1]*U[1,0] + N[2]*U[2,0] + N[3]*U[3,0]))*(N[0]*U[0,2] + N[1]*U[1,2] + N[2]*U[2,2] + N[3]*U[3,2])/(N[0]*U[0,0] + N[1]*U[1,0] + N[2]*U[2,0] + N[3]*U[3,0])**2 + (N[0]*U[0,0] + N[1]*U[1,0] + N[2]*U[2,0] + N[3]*U[3,0])*(N[0]*dUdt[0,2] + N[1]*dUdt[1,2] + N[2]*dUdt[2,2] + N[3]*dUdt[3,2]))]
        ])

        expected_etot_proj = sympy.Matrix([
            [N[0]*(N[0]*dUdt[0,3] - N[0]*forcing_terms["thermal"][0] + N[1]*dUdt[1,3] - N[1]*forcing_terms["thermal"][1] + N[2]*dUdt[2,3] - N[2]*forcing_terms["thermal"][2] + N[3]*dUdt[3,3] - N[3]*forcing_terms["thermal"][3] - (N[0]*U[0,1] + N[1]*U[1,1] + N[2]*U[2,1] + N[3]*U[3,1])*(N[0]*forcing_terms["force"][0,0] + N[1]*forcing_terms["force"][1,0] + N[2]*forcing_terms["force"][2,0] + N[3]*forcing_terms["force"][3,0])/(N[0]*U[0,0] + N[1]*U[1,0] + N[2]*U[2,0] + N[3]*U[3,0]) - (N[0]*U[0,2] + N[1]*U[1,2] + N[2]*U[2,2] + N[3]*U[3,2])*(N[0]*forcing_terms["force"][0,1] + N[1]*forcing_terms["force"][1,1] + N[2]*forcing_terms["force"][2,1] + N[3]*forcing_terms["force"][3,1])/(N[0]*U[0,0] + N[1]*U[1,0] + N[2]*U[2,0] + N[3]*U[3,0]) + 1.0*(N[0]*U[0,1] + N[1]*U[1,1] + N[2]*U[2,1] + N[3]*U[3,1])*(-(DN_DX[0,0]*U[0,0] + DN_DX[1,0]*U[1,0] + DN_DX[2,0]*U[2,0] + DN_DX[3,0]*U[3,0])*(N[0]*U[0,0] + N[1]*U[1,0] + N[2]*U[2,0] + N[3]*U[3,0])*(N[0]*U[0,3] + N[1]*U[1,3] + N[2]*U[2,3] + N[3]*U[3,3]) + (DN_DX[0,0]*U[0,0] + DN_DX[1,0]*U[1,0] + DN_DX[2,0]*U[2,0] + DN_DX[3,0]*U[3,0])*(N[0]*U[0,1] + N[1]*U[1,1] + N[2]*U[2,1] + N[3]*U[3,1])**2 + (DN_DX[0,0]*U[0,0] + DN_DX[1,0]*U[1,0] + DN_DX[2,0]*U[2,0] + DN_DX[3,0]*U[3,0])*(N[0]*U[0,2] + N[1]*U[1,2] + N[2]*U[2,2] + N[3]*U[3,2])**2 - (DN_DX[0,0]*U[0,1] + DN_DX[1,0]*U[1,1] + DN_DX[2,0]*U[2,1] + DN_DX[3,0]*U[3,1])*(N[0]*U[0,0] + N[1]*U[1,0] + N[2]*U[2,0] + N[3]*U[3,0])*(N[0]*U[0,1] + N[1]*U[1,1] + N[2]*U[2,1] + N[3]*U[3,1]) - (DN_DX[0,0]*U[0,2] + DN_DX[1,0]*U[1,2] + DN_DX[2,0]*U[2,2] + DN_DX[3,0]*U[3,2])*(N[0]*U[0,0] + N[1]*U[1,0] + N[2]*U[2,0] + N[3]*U[3,0])*(N[0]*U[0,2] + N[1]*U[1,2] + N[2]*U[2,2] + N[3]*U[3,2]) + (DN_DX[0,0]*U[0,3] + DN_DX[1,0]*U[1,3] + DN_DX[2,0]*U[2,3] + DN_DX[3,0]*U[3,3])*(N[0]*U[0,0] + N[1]*U[1,0] + N[2]*U[2,0] + N[3]*U[3,0])**2)/(params.c_v*(N[0]*U[0,0] + N[1]*U[1,0] + N[2]*U[2,0] + N[3]*U[3,0])**4) + 1.0*(N[0]*U[0,2] + N[1]*U[1,2] + N[2]*U[2,2] + N[3]*U[3,2])*(-(DN_DX[0,1]*U[0,0] + DN_DX[1,1]*U[1,0] + DN_DX[2,1]*U[2,0] + DN_DX[3,1]*U[3,0])*(N[0]*U[0,0] + N[1]*U[1,0] + N[2]*U[2,0] + N[3]*U[3,0])*(N[0]*U[0,3] + N[1]*U[1,3] + N[2]*U[2,3] + N[3]*U[3,3]) + (DN_DX[0,1]*U[0,0] + DN_DX[1,1]*U[1,0] + DN_DX[2,1]*U[2,0] + DN_DX[3,1]*U[3,0])*(N[0]*U[0,1] + N[1]*U[1,1] + N[2]*U[2,1] + N[3]*U[3,1])**2 + (DN_DX[0,1]*U[0,0] + DN_DX[1,1]*U[1,0] + DN_DX[2,1]*U[2,0] + DN_DX[3,1]*U[3,0])*(N[0]*U[0,2] + N[1]*U[1,2] + N[2]*U[2,2] + N[3]*U[3,2])**2 - (DN_DX[0,1]*U[0,1] + DN_DX[1,1]*U[1,1] + DN_DX[2,1]*U[2,1] + DN_DX[3,1]*U[3,1])*(N[0]*U[0,0] + N[1]*U[1,0] + N[2]*U[2,0] + N[3]*U[3,0])*(N[0]*U[0,1] + N[1]*U[1,1] + N[2]*U[2,1] + N[3]*U[3,1]) - (DN_DX[0,1]*U[0,2] + DN_DX[1,1]*U[1,2] + DN_DX[2,1]*U[2,2] + DN_DX[3,1]*U[3,2])*(N[0]*U[0,0] + N[1]*U[1,0] + N[2]*U[2,0] + N[3]*U[3,0])*(N[0]*U[0,2] + N[1]*U[1,2] + N[2]*U[2,2] + N[3]*U[3,2]) + (DN_DX[0,1]*U[0,3] + DN_DX[1,1]*U[1,3] + DN_DX[2,1]*U[2,3] + DN_DX[3,1]*U[3,3])*(N[0]*U[0,0] + N[1]*U[1,0] + N[2]*U[2,0] + N[3]*U[3,0])**2)/(params.c_v*(N[0]*U[0,0] + N[1]*U[1,0] + N[2]*U[2,0] + N[3]*U[3,0])**4))],
            [N[1]*(N[0]*dUdt[0,3] - N[0]*forcing_terms["thermal"][0] + N[1]*dUdt[1,3] - N[1]*forcing_terms["thermal"][1] + N[2]*dUdt[2,3] - N[2]*forcing_terms["thermal"][2] + N[3]*dUdt[3,3] - N[3]*forcing_terms["thermal"][3] - (N[0]*U[0,1] + N[1]*U[1,1] + N[2]*U[2,1] + N[3]*U[3,1])*(N[0]*forcing_terms["force"][0,0] + N[1]*forcing_terms["force"][1,0] + N[2]*forcing_terms["force"][2,0] + N[3]*forcing_terms["force"][3,0])/(N[0]*U[0,0] + N[1]*U[1,0] + N[2]*U[2,0] + N[3]*U[3,0]) - (N[0]*U[0,2] + N[1]*U[1,2] + N[2]*U[2,2] + N[3]*U[3,2])*(N[0]*forcing_terms["force"][0,1] + N[1]*forcing_terms["force"][1,1] + N[2]*forcing_terms["force"][2,1] + N[3]*forcing_terms["force"][3,1])/(N[0]*U[0,0] + N[1]*U[1,0] + N[2]*U[2,0] + N[3]*U[3,0]) + 1.0*(N[0]*U[0,1] + N[1]*U[1,1] + N[2]*U[2,1] + N[3]*U[3,1])*(-(DN_DX[0,0]*U[0,0] + DN_DX[1,0]*U[1,0] + DN_DX[2,0]*U[2,0] + DN_DX[3,0]*U[3,0])*(N[0]*U[0,0] + N[1]*U[1,0] + N[2]*U[2,0] + N[3]*U[3,0])*(N[0]*U[0,3] + N[1]*U[1,3] + N[2]*U[2,3] + N[3]*U[3,3]) + (DN_DX[0,0]*U[0,0] + DN_DX[1,0]*U[1,0] + DN_DX[2,0]*U[2,0] + DN_DX[3,0]*U[3,0])*(N[0]*U[0,1] + N[1]*U[1,1] + N[2]*U[2,1] + N[3]*U[3,1])**2 + (DN_DX[0,0]*U[0,0] + DN_DX[1,0]*U[1,0] + DN_DX[2,0]*U[2,0] + DN_DX[3,0]*U[3,0])*(N[0]*U[0,2] + N[1]*U[1,2] + N[2]*U[2,2] + N[3]*U[3,2])**2 - (DN_DX[0,0]*U[0,1] + DN_DX[1,0]*U[1,1] + DN_DX[2,0]*U[2,1] + DN_DX[3,0]*U[3,1])*(N[0]*U[0,0] + N[1]*U[1,0] + N[2]*U[2,0] + N[3]*U[3,0])*(N[0]*U[0,1] + N[1]*U[1,1] + N[2]*U[2,1] + N[3]*U[3,1]) - (DN_DX[0,0]*U[0,2] + DN_DX[1,0]*U[1,2] + DN_DX[2,0]*U[2,2] + DN_DX[3,0]*U[3,2])*(N[0]*U[0,0] + N[1]*U[1,0] + N[2]*U[2,0] + N[3]*U[3,0])*(N[0]*U[0,2] + N[1]*U[1,2] + N[2]*U[2,2] + N[3]*U[3,2]) + (DN_DX[0,0]*U[0,3] + DN_DX[1,0]*U[1,3] + DN_DX[2,0]*U[2,3] + DN_DX[3,0]*U[3,3])*(N[0]*U[0,0] + N[1]*U[1,0] + N[2]*U[2,0] + N[3]*U[3,0])**2)/(params.c_v*(N[0]*U[0,0] + N[1]*U[1,0] + N[2]*U[2,0] + N[3]*U[3,0])**4) + 1.0*(N[0]*U[0,2] + N[1]*U[1,2] + N[2]*U[2,2] + N[3]*U[3,2])*(-(DN_DX[0,1]*U[0,0] + DN_DX[1,1]*U[1,0] + DN_DX[2,1]*U[2,0] + DN_DX[3,1]*U[3,0])*(N[0]*U[0,0] + N[1]*U[1,0] + N[2]*U[2,0] + N[3]*U[3,0])*(N[0]*U[0,3] + N[1]*U[1,3] + N[2]*U[2,3] + N[3]*U[3,3]) + (DN_DX[0,1]*U[0,0] + DN_DX[1,1]*U[1,0] + DN_DX[2,1]*U[2,0] + DN_DX[3,1]*U[3,0])*(N[0]*U[0,1] + N[1]*U[1,1] + N[2]*U[2,1] + N[3]*U[3,1])**2 + (DN_DX[0,1]*U[0,0] + DN_DX[1,1]*U[1,0] + DN_DX[2,1]*U[2,0] + DN_DX[3,1]*U[3,0])*(N[0]*U[0,2] + N[1]*U[1,2] + N[2]*U[2,2] + N[3]*U[3,2])**2 - (DN_DX[0,1]*U[0,1] + DN_DX[1,1]*U[1,1] + DN_DX[2,1]*U[2,1] + DN_DX[3,1]*U[3,1])*(N[0]*U[0,0] + N[1]*U[1,0] + N[2]*U[2,0] + N[3]*U[3,0])*(N[0]*U[0,1] + N[1]*U[1,1] + N[2]*U[2,1] + N[3]*U[3,1]) - (DN_DX[0,1]*U[0,2] + DN_DX[1,1]*U[1,2] + DN_DX[2,1]*U[2,2] + DN_DX[3,1]*U[3,2])*(N[0]*U[0,0] + N[1]*U[1,0] + N[2]*U[2,0] + N[3]*U[3,0])*(N[0]*U[0,2] + N[1]*U[1,2] + N[2]*U[2,2] + N[3]*U[3,2]) + (DN_DX[0,1]*U[0,3] + DN_DX[1,1]*U[1,3] + DN_DX[2,1]*U[2,3] + DN_DX[3,1]*U[3,3])*(N[0]*U[0,0] + N[1]*U[1,0] + N[2]*U[2,0] + N[3]*U[3,0])**2)/(params.c_v*(N[0]*U[0,0] + N[1]*U[1,0] + N[2]*U[2,0] + N[3]*U[3,0])**4))],
            [N[2]*(N[0]*dUdt[0,3] - N[0]*forcing_terms["thermal"][0] + N[1]*dUdt[1,3] - N[1]*forcing_terms["thermal"][1] + N[2]*dUdt[2,3] - N[2]*forcing_terms["thermal"][2] + N[3]*dUdt[3,3] - N[3]*forcing_terms["thermal"][3] - (N[0]*U[0,1] + N[1]*U[1,1] + N[2]*U[2,1] + N[3]*U[3,1])*(N[0]*forcing_terms["force"][0,0] + N[1]*forcing_terms["force"][1,0] + N[2]*forcing_terms["force"][2,0] + N[3]*forcing_terms["force"][3,0])/(N[0]*U[0,0] + N[1]*U[1,0] + N[2]*U[2,0] + N[3]*U[3,0]) - (N[0]*U[0,2] + N[1]*U[1,2] + N[2]*U[2,2] + N[3]*U[3,2])*(N[0]*forcing_terms["force"][0,1] + N[1]*forcing_terms["force"][1,1] + N[2]*forcing_terms["force"][2,1] + N[3]*forcing_terms["force"][3,1])/(N[0]*U[0,0] + N[1]*U[1,0] + N[2]*U[2,0] + N[3]*U[3,0]) + 1.0*(N[0]*U[0,1] + N[1]*U[1,1] + N[2]*U[2,1] + N[3]*U[3,1])*(-(DN_DX[0,0]*U[0,0] + DN_DX[1,0]*U[1,0] + DN_DX[2,0]*U[2,0] + DN_DX[3,0]*U[3,0])*(N[0]*U[0,0] + N[1]*U[1,0] + N[2]*U[2,0] + N[3]*U[3,0])*(N[0]*U[0,3] + N[1]*U[1,3] + N[2]*U[2,3] + N[3]*U[3,3]) + (DN_DX[0,0]*U[0,0] + DN_DX[1,0]*U[1,0] + DN_DX[2,0]*U[2,0] + DN_DX[3,0]*U[3,0])*(N[0]*U[0,1] + N[1]*U[1,1] + N[2]*U[2,1] + N[3]*U[3,1])**2 + (DN_DX[0,0]*U[0,0] + DN_DX[1,0]*U[1,0] + DN_DX[2,0]*U[2,0] + DN_DX[3,0]*U[3,0])*(N[0]*U[0,2] + N[1]*U[1,2] + N[2]*U[2,2] + N[3]*U[3,2])**2 - (DN_DX[0,0]*U[0,1] + DN_DX[1,0]*U[1,1] + DN_DX[2,0]*U[2,1] + DN_DX[3,0]*U[3,1])*(N[0]*U[0,0] + N[1]*U[1,0] + N[2]*U[2,0] + N[3]*U[3,0])*(N[0]*U[0,1] + N[1]*U[1,1] + N[2]*U[2,1] + N[3]*U[3,1]) - (DN_DX[0,0]*U[0,2] + DN_DX[1,0]*U[1,2] + DN_DX[2,0]*U[2,2] + DN_DX[3,0]*U[3,2])*(N[0]*U[0,0] + N[1]*U[1,0] + N[2]*U[2,0] + N[3]*U[3,0])*(N[0]*U[0,2] + N[1]*U[1,2] + N[2]*U[2,2] + N[3]*U[3,2]) + (DN_DX[0,0]*U[0,3] + DN_DX[1,0]*U[1,3] + DN_DX[2,0]*U[2,3] + DN_DX[3,0]*U[3,3])*(N[0]*U[0,0] + N[1]*U[1,0] + N[2]*U[2,0] + N[3]*U[3,0])**2)/(params.c_v*(N[0]*U[0,0] + N[1]*U[1,0] + N[2]*U[2,0] + N[3]*U[3,0])**4) + 1.0*(N[0]*U[0,2] + N[1]*U[1,2] + N[2]*U[2,2] + N[3]*U[3,2])*(-(DN_DX[0,1]*U[0,0] + DN_DX[1,1]*U[1,0] + DN_DX[2,1]*U[2,0] + DN_DX[3,1]*U[3,0])*(N[0]*U[0,0] + N[1]*U[1,0] + N[2]*U[2,0] + N[3]*U[3,0])*(N[0]*U[0,3] + N[1]*U[1,3] + N[2]*U[2,3] + N[3]*U[3,3]) + (DN_DX[0,1]*U[0,0] + DN_DX[1,1]*U[1,0] + DN_DX[2,1]*U[2,0] + DN_DX[3,1]*U[3,0])*(N[0]*U[0,1] + N[1]*U[1,1] + N[2]*U[2,1] + N[3]*U[3,1])**2 + (DN_DX[0,1]*U[0,0] + DN_DX[1,1]*U[1,0] + DN_DX[2,1]*U[2,0] + DN_DX[3,1]*U[3,0])*(N[0]*U[0,2] + N[1]*U[1,2] + N[2]*U[2,2] + N[3]*U[3,2])**2 - (DN_DX[0,1]*U[0,1] + DN_DX[1,1]*U[1,1] + DN_DX[2,1]*U[2,1] + DN_DX[3,1]*U[3,1])*(N[0]*U[0,0] + N[1]*U[1,0] + N[2]*U[2,0] + N[3]*U[3,0])*(N[0]*U[0,1] + N[1]*U[1,1] + N[2]*U[2,1] + N[3]*U[3,1]) - (DN_DX[0,1]*U[0,2] + DN_DX[1,1]*U[1,2] + DN_DX[2,1]*U[2,2] + DN_DX[3,1]*U[3,2])*(N[0]*U[0,0] + N[1]*U[1,0] + N[2]*U[2,0] + N[3]*U[3,0])*(N[0]*U[0,2] + N[1]*U[1,2] + N[2]*U[2,2] + N[3]*U[3,2]) + (DN_DX[0,1]*U[0,3] + DN_DX[1,1]*U[1,3] + DN_DX[2,1]*U[2,3] + DN_DX[3,1]*U[3,3])*(N[0]*U[0,0] + N[1]*U[1,0] + N[2]*U[2,0] + N[3]*U[3,0])**2)/(params.c_v*(N[0]*U[0,0] + N[1]*U[1,0] + N[2]*U[2,0] + N[3]*U[3,0])**4))],
            [N[3]*(N[0]*dUdt[0,3] - N[0]*forcing_terms["thermal"][0] + N[1]*dUdt[1,3] - N[1]*forcing_terms["thermal"][1] + N[2]*dUdt[2,3] - N[2]*forcing_terms["thermal"][2] + N[3]*dUdt[3,3] - N[3]*forcing_terms["thermal"][3] - (N[0]*U[0,1] + N[1]*U[1,1] + N[2]*U[2,1] + N[3]*U[3,1])*(N[0]*forcing_terms["force"][0,0] + N[1]*forcing_terms["force"][1,0] + N[2]*forcing_terms["force"][2,0] + N[3]*forcing_terms["force"][3,0])/(N[0]*U[0,0] + N[1]*U[1,0] + N[2]*U[2,0] + N[3]*U[3,0]) - (N[0]*U[0,2] + N[1]*U[1,2] + N[2]*U[2,2] + N[3]*U[3,2])*(N[0]*forcing_terms["force"][0,1] + N[1]*forcing_terms["force"][1,1] + N[2]*forcing_terms["force"][2,1] + N[3]*forcing_terms["force"][3,1])/(N[0]*U[0,0] + N[1]*U[1,0] + N[2]*U[2,0] + N[3]*U[3,0]) + 1.0*(N[0]*U[0,1] + N[1]*U[1,1] + N[2]*U[2,1] + N[3]*U[3,1])*(-(DN_DX[0,0]*U[0,0] + DN_DX[1,0]*U[1,0] + DN_DX[2,0]*U[2,0] + DN_DX[3,0]*U[3,0])*(N[0]*U[0,0] + N[1]*U[1,0] + N[2]*U[2,0] + N[3]*U[3,0])*(N[0]*U[0,3] + N[1]*U[1,3] + N[2]*U[2,3] + N[3]*U[3,3]) + (DN_DX[0,0]*U[0,0] + DN_DX[1,0]*U[1,0] + DN_DX[2,0]*U[2,0] + DN_DX[3,0]*U[3,0])*(N[0]*U[0,1] + N[1]*U[1,1] + N[2]*U[2,1] + N[3]*U[3,1])**2 + (DN_DX[0,0]*U[0,0] + DN_DX[1,0]*U[1,0] + DN_DX[2,0]*U[2,0] + DN_DX[3,0]*U[3,0])*(N[0]*U[0,2] + N[1]*U[1,2] + N[2]*U[2,2] + N[3]*U[3,2])**2 - (DN_DX[0,0]*U[0,1] + DN_DX[1,0]*U[1,1] + DN_DX[2,0]*U[2,1] + DN_DX[3,0]*U[3,1])*(N[0]*U[0,0] + N[1]*U[1,0] + N[2]*U[2,0] + N[3]*U[3,0])*(N[0]*U[0,1] + N[1]*U[1,1] + N[2]*U[2,1] + N[3]*U[3,1]) - (DN_DX[0,0]*U[0,2] + DN_DX[1,0]*U[1,2] + DN_DX[2,0]*U[2,2] + DN_DX[3,0]*U[3,2])*(N[0]*U[0,0] + N[1]*U[1,0] + N[2]*U[2,0] + N[3]*U[3,0])*(N[0]*U[0,2] + N[1]*U[1,2] + N[2]*U[2,2] + N[3]*U[3,2]) + (DN_DX[0,0]*U[0,3] + DN_DX[1,0]*U[1,3] + DN_DX[2,0]*U[2,3] + DN_DX[3,0]*U[3,3])*(N[0]*U[0,0] + N[1]*U[1,0] + N[2]*U[2,0] + N[3]*U[3,0])**2)/(params.c_v*(N[0]*U[0,0] + N[1]*U[1,0] + N[2]*U[2,0] + N[3]*U[3,0])**4) + 1.0*(N[0]*U[0,2] + N[1]*U[1,2] + N[2]*U[2,2] + N[3]*U[3,2])*(-(DN_DX[0,1]*U[0,0] + DN_DX[1,1]*U[1,0] + DN_DX[2,1]*U[2,0] + DN_DX[3,1]*U[3,0])*(N[0]*U[0,0] + N[1]*U[1,0] + N[2]*U[2,0] + N[3]*U[3,0])*(N[0]*U[0,3] + N[1]*U[1,3] + N[2]*U[2,3] + N[3]*U[3,3]) + (DN_DX[0,1]*U[0,0] + DN_DX[1,1]*U[1,0] + DN_DX[2,1]*U[2,0] + DN_DX[3,1]*U[3,0])*(N[0]*U[0,1] + N[1]*U[1,1] + N[2]*U[2,1] + N[3]*U[3,1])**2 + (DN_DX[0,1]*U[0,0] + DN_DX[1,1]*U[1,0] + DN_DX[2,1]*U[2,0] + DN_DX[3,1]*U[3,0])*(N[0]*U[0,2] + N[1]*U[1,2] + N[2]*U[2,2] + N[3]*U[3,2])**2 - (DN_DX[0,1]*U[0,1] + DN_DX[1,1]*U[1,1] + DN_DX[2,1]*U[2,1] + DN_DX[3,1]*U[3,1])*(N[0]*U[0,0] + N[1]*U[1,0] + N[2]*U[2,0] + N[3]*U[3,0])*(N[0]*U[0,1] + N[1]*U[1,1] + N[2]*U[2,1] + N[3]*U[3,1]) - (DN_DX[0,1]*U[0,2] + DN_DX[1,1]*U[1,2] + DN_DX[2,1]*U[2,2] + DN_DX[3,1]*U[3,2])*(N[0]*U[0,0] + N[1]*U[1,0] + N[2]*U[2,0] + N[3]*U[3,0])*(N[0]*U[0,2] + N[1]*U[1,2] + N[2]*U[2,2] + N[3]*U[3,2]) + (DN_DX[0,1]*U[0,3] + DN_DX[1,1]*U[1,3] + DN_DX[2,1]*U[2,3] + DN_DX[3,1]*U[3,3])*(N[0]*U[0,0] + N[1]*U[1,0] + N[2]*U[2,0] + N[3]*U[3,0])**2)/(params.c_v*(N[0]*U[0,0] + N[1]*U[1,0] + N[2]*U[2,0] + N[3]*U[3,0])**4))]
        ])

        self.assertTrue(projections["rho"].equals(expected_rho_proj))
        self.assertTrue(projections["momentum"].equals(expected_mom_proj))
        self.assertTrue(projections["energy"].equals(expected_etot_proj))

if __name__ == '__main__':
    suites = KratosUnitTest.KratosSuites

    suites["small"].addTests(KratosUnitTest.TestLoader().loadTestsFromTestCases([CompressibleNavierStokesSymbolicGeneratorUnitTest]))
    suites["validation"].addTests(KratosUnitTest.TestLoader().loadTestsFromTestCases([CompressibleNavierStokesSymbolicGeneratorValidationTest]))

    suites["all"].addTests(suites["small"])
    suites["all"].addTests(suites["validation"])

    KratosUnitTest.runTests(suites)

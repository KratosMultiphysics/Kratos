import os

import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnitTest

from KratosMultiphysics.FluidDynamicsApplication.automatic_differentiation.compressible_navier_stokes \
    .compressible_navier_stokes_symbolic_generator import CompressibleNavierStokesSymbolicGenerator

import KratosMultiphysics.FluidDynamicsApplication

class CompressibleNavierStokesSymbolicGeneratorFormulationTest(KratosUnitTest.TestCase):
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
        module_name = "compressible_automatic_differentiation" + "." + generated_file_name[:-3]

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
        with KratosUnitTest.WorkFolderScope("compressible_automatic_differentiation", __file__):
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

if __name__ == '__main__':
    KratosUnitTest.main()
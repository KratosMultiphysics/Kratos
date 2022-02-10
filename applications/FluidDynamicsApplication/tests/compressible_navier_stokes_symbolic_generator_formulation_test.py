import os

import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnitTest

from KratosMultiphysics.FluidDynamicsApplication.symbolic_generation.compressible_navier_stokes \
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
    def _GetGeneratorSettings(cls, geometry):
        """
        Returns the Kratos Parameters for the symbolic generator
        """
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
        params = cls._GetGeneratorSettings(geometry)
        generator = CompressibleNavierStokesSymbolicGenerator(params)
        generator.Generate()
        generator.Write()

        return generator.output_filename

    @classmethod
    def _ImportSubTestSuite(cls, generated_file_name):
        module_name = "compressible_symbolic_generation" + "." + generated_file_name[:-3]

        try:
            test_module = __import__(module_name, fromlist=module_name)
        except ModuleNotFoundError as err:
            raise RuntimeError("Failed to import generated file:", generated_file_name) from err
        return test_module.SubTestSuite()


    def _RunSubTestSuite(self, sub_testsuite,print_results):
        tests = [subtest_name for subtest_name in dir(sub_testsuite)
            if callable(getattr(sub_testsuite, subtest_name)) and subtest_name.startswith("test_")]

        for subtest_name in tests:
            try:
                result, reference = getattr(sub_testsuite, subtest_name)()
            except Exception as e:
                raise RuntimeError("Error in sub-test " + subtest_name) from e

            if print_results:
                print(result)

            self.assertVectorAlmostEqual(result, reference, msg="Failure in sub-test " + subtest_name)

    def _RunTest(self, geometry, print_results=False, cleanup=True):
        """
        Runs the test.

        kwargs
        ------
        - cleanup: Instructs the test whether to remove the generated code
        - geometry: Choice of geometry. Format is xDyN, with x,y integers
        - print_results: Whether to print all the results even if they are correct
        """

        generated_file = self._Generate(geometry)

        if cleanup:
            self.files_to_remove.append(os.path.abspath(generated_file))

        sub_testsuite = self._ImportSubTestSuite(generated_file)
        self._RunSubTestSuite(sub_testsuite, print_results)

    def test_SymbolicTriangle(self):
        with KratosUnitTest.WorkFolderScope("compressible_symbolic_generation", __file__):
            self._RunTest("2D3N")

    def test_SymbolicQuadrilateral(self):
        with KratosUnitTest.WorkFolderScope("compressible_symbolic_generation", __file__):
            self._RunTest("2D4N")

    def test_SymbolicTetrahedron(self):
        with KratosUnitTest.WorkFolderScope("compressible_symbolic_generation", __file__):
            self._RunTest("3D4N")

if __name__ == '__main__':
    KratosUnitTest.main()

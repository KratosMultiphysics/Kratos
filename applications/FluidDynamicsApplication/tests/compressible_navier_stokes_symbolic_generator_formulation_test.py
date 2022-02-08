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
        module_name = "compressible_symbolic_generation"

        try:
            test_module = __import__(module_name + "." + generated_file_name[:-3], fromlist=module_name)
        except ModuleNotFoundError as err:
            print("Failed to import generated file:", module_name)
            raise err
        return test_module.SubTestSuite()


    def _RunSubTestSuite(self, sub_testsuite,print_results):
        tests = [testname for testname in dir(sub_testsuite)
            if callable(getattr(sub_testsuite, testname)) and testname.startswith("test_")]

        for testname in tests:
            try:
                result, reference = getattr(sub_testsuite, testname)()
            except Exception as e:
                raise RuntimeError("Error in test " + testname) from e

            if print_results:
                print(result)

            try:
                self.assertVectorAlmostEqual(result, reference)
            except AssertionError as e:
                raise AssertionError("AssertionError in " + testname) from e

    def _RunTest(self, **kwargs):
        """
        Runs the test.

        kwargs
        ------

        - regenerate: Instructs the test whether to call the symbolic generator
            or use a pre-existing source file.
        - cleanup: Instructs the test whether to remove generated code and
            compiled binary
        - compiler: Choice of compiler. Must accept GCC-like commands (e.g -00 -g)
        - geometry: Choice of geometry. Format is xDyN, with x,y integers
        - print_results: Whether to print all the results even if they are correct
        """
        regenerate = True if "regenerate" not in kwargs else kwargs["regenerate"]
        cleanup    = True if "cleanup" not in kwargs else kwargs["cleanup"]
        geometry   = kwargs["geometry"]
        print_results = False if "print_results" not in kwargs else kwargs["print_results"]


        if regenerate:
            generated_file = self._Generate(geometry)
        else:
            generated_file = self._GetGeneratorSettings(geometry)["output_filename"].GetString()

        if cleanup:
            self.files_to_remove.append(os.path.abspath(generated_file))

        sub_testsuite = self._ImportSubTestSuite(generated_file)
        self._RunSubTestSuite(sub_testsuite, print_results)

    def test_SymbolicTrinagle(self):
        args = {
            "regenerate": True,
            "cleanup": True,
            "print_results": False,
            "geometry": "2D3N"
        }

        with KratosUnitTest.WorkFolderScope("compressible_symbolic_generation", __file__):
            self._RunTest(**args)

    def test_SymbolicQuadrilateral(self):
        args = {
            "regenerate": True,
            "cleanup": True,
            "print_results": False,
            "geometry": "2D4N"
        }

        with KratosUnitTest.WorkFolderScope("compressible_symbolic_generation", __file__):
            self._RunTest(**args)

    def test_SymbolicTetrahedron(self):
        args = {
            "regenerate": True,
            "cleanup": True,
            "print_results": False,
            "geometry": "3D4N"
        }

        with KratosUnitTest.WorkFolderScope("compressible_symbolic_generation", __file__):
            self._RunTest(**args)

if __name__ == '__main__':
    KratosUnitTest.main()

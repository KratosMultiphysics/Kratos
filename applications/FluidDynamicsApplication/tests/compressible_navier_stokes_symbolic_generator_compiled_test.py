import os
import sys
import unittest

import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnitTest

from KratosMultiphysics.FluidDynamicsApplication.symbolic_generation.compressible_navier_stokes \
    .compressible_navier_stokes_symbolic_generator import CompressibleNavierStokesSymbolicGenerator


class CompressibleNavierStokesSymbolicGeneratorCompilationTest(KratosUnitTest.TestCase):
    def setUp(self):
        self.files_to_remove = []

    def tearDown(self):
        for file_name in self.files_to_remove:
            try:
                os.remove(file_name)
            except FileNotFoundError:
                pass
        self.files_to_remove = []

    @classmethod
    def _GetGeneratorSettings(cls, geometry):
        """
        Returns the Kratos Parameters for the symbolic generator
        """
        geometry_name = {
            "2D3N": "triangle",
            "2D4N": "quadrilateral",
            "3D4N": "hexahedron"
        }[geometry]

        parameters = KratosMultiphysics.Parameters("""
        {
            "geometry": "PLEASE SPECIFY A GEOMETRY",
            "template_filename" : "PLEASE SPECIFY A TEMPLATE FILE",
            "output_filename"   : "PLEASE SPECIFY AN OUTPUT FILE"
        }""")

        parameters["geometry"].SetString(geometry_name)
        parameters["template_filename"].SetString("compilable_{}.template".format(geometry))
        parameters["output_filename"].SetString("symbolic_test_{}.cpp".format(geometry))

        return parameters

    @classmethod
    def _ReadTestStdErr(cls):
        try:
            with open("test_stderr.log", "r") as f:
                return f.read()
        except Exception as exc:
            return "Failed to parse stderr:\n" + str(exc)

    def _GenerateCompileAndRun(self, **kwargs):
        """
        Runs the test.

        kwargs
        ------

        - regenerate: Instructs the test whether to call the symbolic generator
            or use a pre-existing source file.
        - recompile: Instructs the test whether to call the compiler
            or use a pre-existing binary file. If regenerate is enabled,
            recompile will be overriden to True.
        - cleanup: Instructs the test whether to remove generated code and
            compiled binary
        - compiler: Choice of compiler. Must accept GCC-like commands (e.g -00 -g)
        - geometry: Choice of geometry. Format is xDyN, with x,y integers
        - print_results: Whether to print all the results even if they are correct
        """
        regenerate = True if "regenerate" not in kwargs else kwargs["regenerate"]
        recompile  = True if "recompile" not in kwargs else kwargs["recompile"]
        cleanup    = True if "cleanup" not in kwargs else kwargs["cleanup"]
        compiler   = "g++" if "compiler" not in kwargs else kwargs["compiler"]
        geometry   = "2D3N" if "geometry" not in kwargs else kwargs["geometry"]
        print_results = False if "print_results" not in kwargs else kwargs["print_results"]

        recompile = recompile or regenerate
        parameters = self._GetGeneratorSettings(geometry)

        # Adding files to be removed
        if cleanup:
            self.files_to_remove.append(os.path.abspath(parameters["output_filename"].GetString()))
            self.files_to_remove.append(os.path.abspath("test.out"))
            self.files_to_remove.append(os.path.abspath("test_stderr.log"))

        # Checking Compiler exists
        if recompile:
            print("\n----------------Checking compiler exists----------------")
            errcode = os.system("{} --version 2> test_stderr.log".format(compiler))
            self.assertEqual(errcode, 0, "Compiler {} not available.\nStderr:\n{}".format(compiler, self._ReadTestStdErr()))

        # Generation
        if regenerate:
            print("\n--------------------Generating--------------------------")
            generator = CompressibleNavierStokesSymbolicGenerator(parameters)
            generator.Generate()
            generator.Write()

        # Compilation
        if recompile:
            outfile = parameters["output_filename"].GetString()

            print("\n---------------------Compiling---------------------------")
            flags = "-g -O0 -Wall -Wextra -fsanitize=address -fno-omit-frame-pointer -fsanitize=undefined"
            errcode = os.system("{} {} -o test.out {} 2> test_stderr.log".format(compiler, outfile, flags))
            self.assertEqual(errcode, 0, "Compiler returned error code {}\nStderr:\n{}".format(errcode, self._ReadTestStdErr()))

        # Testing
        print("\n----------------------Testing----------------------------")
        command = "./test.out"
        if print_results:
            command = command + " --print-results"
        command = command + " 2> test_stderr.log" # stderr redirected in order to print it later

        errcode = os.system(command)
        self.assertEqual(errcode, 0, "Tests returned error code {}.\nStderr:\n{}\n".format(errcode, self._ReadTestStdErr()))

        return 0

    def _GetCompilerFromEnvironment(self):
        """
        Gets the compiler from environment variables (Useful in continuous integration runs).
        """
        if "CXX" not in os.environ:
            return None

        compiler = os.environ["CXX"]

        if len(compiler) == 0:
            return None

        return compiler

    @unittest.skipIf("linux" not in sys.platform, "This test is only available on linux")
    def test_Quadrilateral(self):
        args = {
            "regenerate": True,
            "recompile": True,
            "cleanup": True,
            "print_results": False,
            "compiler": "g++",
            "geometry": "2D4N"
        }

        ci_compiler = self._GetCompilerFromEnvironment()
        if ci_compiler is not None:
            args["compiler"] = ci_compiler

        with KratosUnitTest.WorkFolderScope("compressible_symbolic_generation", __file__):
            self._GenerateCompileAndRun(**args)

    @unittest.skipIf("linux" not in sys.platform, "This test is only available on linux")
    def test_Triangle(self):
        args = {
            "regenerate": True,
            "recompile": True,
            "cleanup": True,
            "print_results": False,
            "compiler": "g++",
            "geometry": "2D3N"
        }

        ci_compiler = self._GetCompilerFromEnvironment()
        if ci_compiler is not None:
            args["compiler"] = ci_compiler

        with KratosUnitTest.WorkFolderScope("compressible_symbolic_generation", __file__):
            self._GenerateCompileAndRun(**args)


if __name__ == '__main__':
    KratosUnitTest.main()

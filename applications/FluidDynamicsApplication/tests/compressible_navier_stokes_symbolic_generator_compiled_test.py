import os
import sys
import subprocess
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

    def _GetGeneratorSettings(geometry):
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
        - dump_values: Whether to print all the results or not
        """
        regenerate = True if "regenerate" not in kwargs else kwargs["regenerate"]
        recompile = True if "recompile" not in kwargs else kwargs["recompile"]
        cleanup = True if "cleanup" not in kwargs else kwargs["cleanup"]
        compiler = "g++" if "compiler" not in kwargs else kwargs["compiler"]
        geometry = "2D3N" if "geometry" not in kwargs else kwargs["geometry"]
        dump_values = False if "dump_values" not in kwargs else kwargs["dump_values"]

        recompile = recompile or regenerate
        parameters = self._GetGeneratorSettings(geometry)

        # Adding files to be removed
        if cleanup:
            self.files_to_remove.append(os.path.abspath(parameters["output_filename"].GetString()))
            self.files_to_remove.append(os.path.abspath("test.out"))

        # Checking Compiler exists
        if recompile:
            print("\n----------------Checking compiler exists----------------")
            errcode = os.system("{} --version".format(compiler))
            self.assertEqual(errcode, 0, "Compiler {} not available".format(compiler))

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
            errcode = os.system("{} {} -o test.out -g -O0 -Wall -Wextra -fsanitize=address -fno-omit-frame-pointer -fsanitize=undefined".format(compiler, outfile))
            self.assertEqual(errcode, 0, "Compiler returned error code {}\n".format(errcode))

        # Testing
        print("\n----------------------Testing----------------------------")
        errcode = os.system("./test.out" + (" --dump" if dump_values else ""))
        self.assertEqual(errcode, 0, "Tests returned error code {}\n".format(errcode))

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
            "dump_values": False,
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
            "dump_values": False,
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

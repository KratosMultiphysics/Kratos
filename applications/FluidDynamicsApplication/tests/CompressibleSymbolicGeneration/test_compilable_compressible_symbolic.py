import os
import sys

import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnitTest

from KratosMultiphysics.FluidDynamicsApplication.symbolic_generation.compressible_navier_stokes \
    .compressible_navier_stokes_symbolic_generator import CompressibleNavierStokesSymbolicGenerator

class CompressibleNavierStokesSymbolicGeneratorCompilationTest(KratosUnitTest.TestCase):
    def generate_compile_run(self, **kwargs):
        regenerate = True if "regenerate" not in kwargs else kwargs["regenerate"]
        recompile = True if "recompile" not in kwargs else kwargs["recompile"]
        compiler = "g++" if "compiler" not in kwargs else kwargs["compiler"]
        geometry = "2D3N" if "geometry" not in kwargs else kwargs["geometry"]
        dump_values = False if "dump_values" not in kwargs else kwargs["dump_values"]

        # Generation
        geometry_name = {
            "2D3N": "triangle",
            "2D4N": "quadrilateral",
            "3D4N": "hexahedron"
        }[geometry]

        parameters = KratosMultiphysics.Parameters("""
        {
            "geometry": "__geometry_name__",
            "template_filename" : "compilable___xDyN__.template",
            "output_filename"   : "symbolic_test___xDyN__.cpp"
        }""")

        parameters["geometry"].SetString(geometry_name)
        parameters["template_filename"].SetString("compilable_{}.template".format(geometry))
        parameters["output_filename"].SetString("symbolic_test_{}.cpp".format(geometry))

        if regenerate:
            print("\n--------------------Generating--------------------------")
            generator = CompressibleNavierStokesSymbolicGenerator(parameters)
            generator.Generate()
            generator.Write()

        # Compilation
        if recompile or regenerate:
            print("\n----------------Checking compiler exists----------------")
            errcode = os.system("{} --version".format(compiler))
            self.assertEqual(errcode, 0, "Compiler {} not available".format(compiler))

            outfile = parameters["output_filename"].GetString()

            print("\n---------------------Compiling---------------------------")
            errcode = os.system("{} {} -o test.out -g -O0 -Wall -Wextra -fsanitize=address -fno-omit-frame-pointer -fsanitize=undefined".format(compiler, outfile))
            self.assertEqual(errcode, 0, "Compiler returned error code {}\n".format(errcode))

        # Testing
        print("\n----------------------Testing----------------------------")
        errcode = os.system("./test.out" + (" --dump" if dump_values else ""))
        self.assertEqual(errcode, 0, "Tests returned error code {}\n".format(errcode))

        return 0

    def test_Quadrilateral(self):
        args = {
            "regenerate": True,
            "recompile": True,
            "compiler": "g++",
            "geometry": "2D4N",
            "dump_values": False
        }
        with KratosUnitTest.WorkFolderScope(".", __file__):
            self.generate_compile_run(**args)

    def test_Triangle(self):
        args = {
            "regenerate": True,
            "recompile": True,
            "compiler": "g++",
            "geometry": "2D4N",
            "dump_values": False
        }
        with KratosUnitTest.WorkFolderScope(".", __file__):
            self.generate_compile_run(**args)


if __name__ == '__main__':
    KratosUnitTest.main()

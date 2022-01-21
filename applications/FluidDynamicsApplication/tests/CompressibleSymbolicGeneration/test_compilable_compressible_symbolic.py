import pathlib
import os

import KratosMultiphysics
from KratosMultiphysics.FluidDynamicsApplication.symbolic_generation.compressible_navier_stokes \
    .compressible_navier_stokes_symbolic_generator import CompressibleNavierStokesSymbolicGenerator

# DATA ENTRY
regenerate = True
recompile = True
compiler = "g++"
geometry = "2D3N"
dump_values = False

# Generation
geometry_name = {
    "2D3N": "triangle",
    "2D4N": "quadrilateral",
}[geometry]

parameters = KratosMultiphysics.Parameters("""
{
    "geometry": "@geometry_name",
    "template_filename" : "compilable_@geometry_xDyN.template",
    "output_filename"   : "symbolic_test_@geometry_xDyN.cpp"
}
""".replace("@geometry_xDyN", geometry).replace("@geometry_name", geometry_name))

path = pathlib.Path(__file__).parent
os.chdir(path)

if regenerate:
    print("\n--------------------Generating--------------------------")
    generator = CompressibleNavierStokesSymbolicGenerator(parameters)
    generator.Generate()
    generator.Write()

# Compilation
if recompile or regenerate:
    print("\n----------------Checking compiler exists----------------")
    errcode = os.system("{} --version".format(compiler))
    if errcode != 0:
        print("\n Compiler {} not available".format(compiler))
        exit(errcode)

    outfile = parameters["output_filename"].GetString()

    print("\n---------------------Compiling---------------------------")
    errcode = os.system("{} {} -o test.out -g -O0 -Wall -Wextra -fsanitize=address -fno-omit-frame-pointer -fsanitize=undefined".format(compiler, outfile))
    if errcode != 0:
        print("\nCompiler returned error code {}\n".format(errcode))
        exit(errcode)

# Testing
print("\n----------------------Testing----------------------------")
errcode = os.system("./test.out" + (" --dump" if dump_values else ""))
if errcode != 0:
    print("\nTests returned error code {}\n".format(errcode))
    exit(errcode)

exit(0)

import pathlib
import os

import KratosMultiphysics
from KratosMultiphysics.FluidDynamicsApplication.symbolic_generation.compressible_navier_stokes \
    .compressible_navier_stokes_symbolic_generator import CompressibleNavierStokesSymbolicGenerator

# DATA ENTRY
regenerate = True
recompile = True
compiler = "g++"
test_args = "sod"

parameters = KratosMultiphysics.Parameters("""
{
    "geometry": "quadrilateral",
    "template_filename" : "compilable_2D4N.template",
    "output_filename"   : "symbolic_test.cpp"
}
""")

# Generation
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
        print("\nGCC returned error code {}\n".format(errcode))
        exit(errcode)

# Testing
print("\n----------------------Testing----------------------------")
errcode = os.system("./test.out {}".format(test_args))
if errcode != 0:
    print("\nTests returned error code {}\n".format(errcode))
    exit(errcode)

exit(0)

import pathlib
import os

import KratosMultiphysics
from KratosMultiphysics.FluidDynamicsApplication.symbolic_generation.compressible_navier_stokes \
    .compressible_navier_stokes_symbolic_generator import CompressibleNavierStokesSymbolicGenerator

parameters = KratosMultiphysics.Parameters("""
{
    "geometry": "quadrilateral",
    "template_filename" : "compilable_template.cpp",
    "output_filename"   : "test_2D4N_symbolic.cpp"
}
""")

path = pathlib.Path(__file__).parent
os.chdir(path)

generator = CompressibleNavierStokesSymbolicGenerator(parameters)
generator.Generate()
generator.Write()

# Compilation
outfile = parameters["output_filename"].GetString()

errcode = os.system("g++ {} -o test.out".format(outfile))
if errcode != 0:
    print("\nGCC returned error code {}\n".format(errcode))
    exit(errcode)

errcode = os.system("./test.out")
exit(errcode)

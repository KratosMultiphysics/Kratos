import pathlib
import os

import KratosMultiphysics
from KratosMultiphysics.FluidDynamicsApplication.symbolic_generation.compressible_navier_stokes. \
    compressible_navier_stokes_symbolic_generator import CompressibleNavierStokesSymbolicGenerator

parameters = KratosMultiphysics.Parameters("""
{
    "geometry": "quadrilateral",
    "template_filename" : "templates/compressible_navier_stokes_explicit_cpp_template_with_integration.cpp",
    "output_filename"   : "templates/compressible_explicit_navier_stokes_quad.cpp"
}
""")

path = pathlib.Path(__file__).parent
os.chdir(path)

generator = CompressibleNavierStokesSymbolicGenerator(parameters)
generator.Generate()
generator.Write()

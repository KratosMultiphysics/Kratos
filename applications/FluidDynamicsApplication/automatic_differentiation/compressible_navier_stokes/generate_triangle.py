import pathlib
import os

import KratosMultiphysics
from KratosMultiphysics.FluidDynamicsApplication.compressible_navier_stokes \
    .compressible_navier_stokes_symbolic_generator import CompressibleNavierStokesSymbolicGenerator

parameters = KratosMultiphysics.Parameters("""
{
    "geometry": "triangle",
    "template_filename" : "templates/compressible_navier_stokes_explicit_template_2D3N.cpp",
    "output_filename"   : "compressible_navier_stokes_explicit_2D3N.cpp"
}
""")

path = pathlib.Path(__file__).parent
os.chdir(path)

generator = CompressibleNavierStokesSymbolicGenerator(parameters)
generator.Generate()
generator.Write()
import pathlib
import os

import KratosMultiphysics
from KratosMultiphysics.FluidDynamicsApplication.symbolic_generation.compressible_navier_stokes \
    .compressible_navier_stokes_condition_symbolic_generator import CompressibleNavierStokesConditionSymbolicGenerator


parameters = KratosMultiphysics.Parameters("""
{
    "geometry": "line",
    "template_filename" : "templates/compressible_navier_stokes_explicit_condition_template_2D2N.cpp",
    "output_filename"   : "compressible_navier_stokes_explicit_condition_2D2N.cpp"
}
""")

path = pathlib.Path(__file__).parent
os.chdir(path)

generator = CompressibleNavierStokesConditionSymbolicGenerator(parameters)
generator.Generate()
generator.Write()

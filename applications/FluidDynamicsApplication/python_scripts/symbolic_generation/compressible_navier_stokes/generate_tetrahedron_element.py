import pathlib
import os

import KratosMultiphysics
from KratosMultiphysics.FluidDynamicsApplication.symbolic_generation.compressible_navier_stokes \
    .compressible_navier_stokes_element_symbolic_generator import CompressibleNavierStokesElementSymbolicGenerator


parameters = KratosMultiphysics.Parameters("""
{
    "geometry": "tetrahedron",
    "template_filename" : "templates/compressible_navier_stokes_explicit_template_3D4N.cpp",
    "output_filename"   : "compressible_navier_stokes_explicit_3D4N.cpp"
}
""")

path = pathlib.Path(__file__).parent
os.chdir(path)

generator = CompressibleNavierStokesElementSymbolicGenerator(parameters)
generator.Generate()
generator.Write()

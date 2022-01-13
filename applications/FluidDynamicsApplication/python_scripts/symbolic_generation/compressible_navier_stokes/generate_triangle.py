import pathlib
import os

import KratosMultiphysics
from KratosMultiphysics.FluidDynamicsApplication.symbolic_generation.compressible_navier_stokes. \
    compressible_navier_stokes_symbolic_generator import CompressibleNavierStokesSymbolicGenerator


parameters = KratosMultiphysics.Parameters("""
{
    "2D" :
    {
        "geometry": "triangle",
        "template_filename" : "templates/compressible_navier_stokes_explicit_cpp_template_with_integration.cpp",
        "output_filename"   : "templates/compressible_explicit_navier_stokes.cpp"
    },
    "3D" :
    {
        "geometry": "tetrahedron",
        "template_filename" : "templates/compressible_explicit_navier_stokes.cpp",
        "output_filename"   : "templates/compressible_explicit_navier_stokes.cpp"
    }
}
""")

path = pathlib.Path(__file__).parent
os.chdir(path)

generator_2d = CompressibleNavierStokesSymbolicGenerator(parameters["2D"])
generator_2d.Generate()
generator_2d.Write()

generator_3d = CompressibleNavierStokesSymbolicGenerator(parameters["3D"])
generator_3d.Generate()
generator_3d.Write()
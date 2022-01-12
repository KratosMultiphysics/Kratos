import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.FluidDynamicsApplication.symbolic_generation.compressible_navier_stokes.compressible_navier_stokes_symbolic_generator import SymbolicGenerator


class TestSymbolicGeneration(KratosUnittest.TestCase):
    def SetUp(self):
        self.scope = KratosUnittest.WorkFolderScope()
        self.scope.__enter__()

    def tearDown(self):
        self.scope.__exit__()

    def test_Generator(self):
        parameters = KratosMultiphysics.Parameters("""
        {
            "2D" :
            {
                "geometry": "triangle",
                "template_filename" : "test_compressible_symbolic.template",
                "output_filename"   : "test_compressible_symbolic.result"
            },
            "3D" :
            {
                "geometry": "tetrahedron",
                "template_filename" : "test_compressible_symbolic.result",
                "output_filename"   : "test_compressible_symbolic.result"
            }
        }
        """)

        generator_2d = SymbolicGenerator(parameters["2D"])
        generator_2d.Generate()
        generator_2d.Write()

        generator_3d = SymbolicGenerator(parameters["3D"])
        generator_3d.Generate()
        generator_3d.Write()


if __name__ == '__main__':
    KratosUnittest.main()

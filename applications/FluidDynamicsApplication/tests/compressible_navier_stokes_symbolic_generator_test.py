import os
import hashlib

import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest

from KratosMultiphysics.FluidDynamicsApplication.symbolic_generation.compressible_navier_stokes \
    .compressible_navier_stokes_symbolic_generator import CompressibleNavierStokesSymbolicGenerator


class CompressibleNavierStokesSymbolicGeneratorTest(KratosUnittest.TestCase):
    @classmethod
    def _GetParameters(cls, geometry_name):
        return KratosMultiphysics.Parameters("""
        {
            "Parameters" : {
                "geometry": "{geometry}",
                "template_filename" : "CompressibleSymbolicGeneration/test_compressible_symbolic_{geometry}.template",
                "output_filename"   : "CompressibleSymbolicGeneration/test_compressible_symbolic_{geometry}.result",
                "echo_level" : 1
            },
            "reference": "CompressibleSymbolicGeneration/test_compressible_symbolic_{geometry}.reference"
        }
        """.replace("{geometry}", geometry_name)
        )

    @classmethod
    def RemoveTestOutput(cls, parameters):
        try:
            os.remove(parameters["Parameters"]["output_filename"].GetString())
        except FileNotFoundError:
            pass

    def _RunAndHashResult(self, geometry_name):
        """Runs the symbolic generator and computes the sha512 hash of the output file"""
        parameters = self._GetParameters(geometry_name)

        with KratosUnittest.WorkFolderScope(".", __file__):
            self.RemoveTestOutput(parameters)

            generator = CompressibleNavierStokesSymbolicGenerator(parameters["Parameters"])
            generator.Generate()
            generator.Write()

            with open(parameters["Parameters"]["output_filename"].GetString(), "r") as f:
                result = hashlib.sha512(f.read().encode("utf-8"))

            self.RemoveTestOutput(parameters)

        return result.hexdigest()


    def testGeneratorTetra(self):
        reference_hash = "9fcd17fc370c71a668278d7c2e0753e465562a4fbfebd180a8dbfda0c0422948bec06d23c9dbf7d956201a499dbdeb4d33305620d499a8f5b3ad614bfcf64d52"
        obtained_hash  = self._RunAndHashResult("tetrahedron")
        self.assertEqual(reference_hash, obtained_hash)

    def testGeneratorTriangle(self):
        reference_hash = "6beb5916b633dfd6a4a8ac82ddc0ab9900410059d9df8f9fbd09b2c013d52f1894988804e3c1f77bc06e5989e626d10d3abd54e24d04b0166f084a1626a40a6d"
        obtained_hash  = self._RunAndHashResult("triangle")
        self.assertEqual(reference_hash, obtained_hash)

    def testGeneratorQuad(self):
        reference_hash = "4c8f0284d7b0f40ffd49424cf174f31b21c91e1e6fd5ea98870db38ad0f04c76f53d5c5bef7fa01fcdb8118132e248ba10388ff0d747800277e0e9bac4ca2fab"
        obtained_hash  = self._RunAndHashResult("quadrilateral")
        self.assertEqual(reference_hash, obtained_hash)



if __name__ == '__main__':
    KratosUnittest.main()

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
        reference_hash = "63c2a914c7127f4272bb653739aef97507ddf396f21a21bf4a31198640878fa93b8a3f7028f12cdef1ca8deda8b1bf9fa8331211b026ece3680d0dc40cc733df"
        obtained_hash  = self._RunAndHashResult("tetrahedron")
        self.assertEqual(reference_hash, obtained_hash)

    def testGeneratorTriangle(self):
        reference_hash = "8fd1983a5cd95f1c5c1504f3a5096672d1e8e714c0bde1c12f42d8ca2f04f58e6d01530d66a1cdb0832b09404663f295941614c5ee1f47aa28b8ee1c97f7d96f"
        obtained_hash  = self._RunAndHashResult("triangle")
        self.assertEqual(reference_hash, obtained_hash)

    def testGeneratorQuad(self):
        reference_hash = "86d19a9d4b6ac20f0a0fcdabb79cecec4d312577d98d4ec4cedc3bba56c139656b56212aa40b81ee46c820980af2c5be87b087511ec7fe475acac4c0f9bf6bea"
        obtained_hash  = self._RunAndHashResult("quadrilateral")
        self.assertEqual(reference_hash, obtained_hash)



if __name__ == '__main__':
    KratosUnittest.main()

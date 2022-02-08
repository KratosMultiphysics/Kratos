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
        reference_hash = "8aa5ef176f59f717c02920e35b6cb301d80e4881339d5294671e7233417681df0d0e96e641a6c920b01770abc88aef070cfa894ef9012fdfb2778ea005be6537"
        obtained_hash  = self._RunAndHashResult("tetrahedron")
        self.assertEqual(reference_hash, obtained_hash)

    def testGeneratorTriangle(self):
        reference_hash = "ea5c4837899d6064d4b080b5dbb16ee348e078248d88df4814f8791868b066afcadd576ad3a80946e48e835d1fcdcd4a69ef6304e727102ea557ef03ef9993f9"
        obtained_hash  = self._RunAndHashResult("triangle")
        self.assertEqual(reference_hash, obtained_hash)

    def testGeneratorQuad(self):
        reference_hash = "94dad96b43402bcd2175b88d0338483c838b7ab2f6fa3f820e9ce6374ff0bed31eb23cd1d645c46a9e1583cbced92544326bf5e80de9b2b96659f56cba8cdb34"
        obtained_hash  = self._RunAndHashResult("quadrilateral")
        self.assertEqual(reference_hash, obtained_hash)



if __name__ == '__main__':
    KratosUnittest.main()

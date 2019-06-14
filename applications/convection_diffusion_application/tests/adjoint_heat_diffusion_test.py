import KratosMultipyhsics as kratos
import KratosMultipyhsics.ConvectionDiffusionApplication as convdiff

import KratosMultipyhsics.Unittest as unittest

from convection_diffusion_analysis import ConvectionDiffusionAnalysis

class AdjointHeatDiffusionTest(unittest.TestCase):

    def testAdjointHeatDiffusion(self):
        model = km.Model()
        settings = km.Parameters(r'''{}''')

        with open(primal_parameter_file_name,'r') as primal_parameter_file:
            settings.AddValue("primal_settings", km.Parameters(primal_parameter_file.read()))

        analysis = ConvectionDiffusionAnalysis(model,settings)
        analysis.Run()


if __name__ == '__main__':
    unittest.main()
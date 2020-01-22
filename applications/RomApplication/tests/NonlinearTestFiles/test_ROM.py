from __future__ import print_function, absolute_import, division
import KratosMultiphysics
import numpy as np

import KratosMultiphysics.KratosUnittest as KratosUnittest
from .MainKratosROM import ConvDiffROM


class ROMLinearTest(KratosUnittest.TestCase):
#########################################################################################

    def test_ConvDiff_ROM_2D(self):

        with KratosUnittest.WorkFolderScope(".", __file__):
            with open("ProjectParametersROM.json",'r') as parameter_file:
                parameters = KratosMultiphysics.Parameters(parameter_file.read())

            model = KratosMultiphysics.Model()
            simulation = ConvDiffROM(model,parameters)
            simulation.Run()
            ObtainedOutput = simulation.EvaluateQuantityOfInterest()
            ExpectedOutput = [-12.148084271927585, 14.819025485595288, -155.5281654573427, -9.353101257284601, 10.72554942206833, -124.59170838744552, 86.38273570629585, 16.913999272914044, 111.60080675691219, 37.69212146428647, 90.0, 90.0, 232.55768788233212, 253.10649201246736, 90.0, 90.0]
            NodalArea = [0.03703703702962964, 0.05555555555555557, 0.055555555555555566, 0.11111111112222223, 0.05555555555000001, 0.05555555555, 0.11111111112222222, 0.11111111112222222, 0.11111111112222222, 0.01851851851481482, 0.018518518514814813, 0.055555555549999996, 0.055555555549999996, 0.05555555555555555, 0.05555555555555555, 0.03703703702962963]
            UP=0
            DOWN=0
            for i in range (16):
                UP += (NodalArea[i]*(    (ExpectedOutput[i] - ObtainedOutput[i]   )**2)  )
                DOWN +=  NodalArea[i]
            L2 = np.sqrt(UP/DOWN)
            self.assertLess(L2, 1.0e-4)


##########################################################################################

if __name__ == '__main__':
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
    KratosUnittest.main()

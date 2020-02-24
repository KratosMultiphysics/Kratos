from __future__ import print_function, absolute_import, division
import KratosMultiphysics
import numpy as np

from .MainKratosROM import TestConvectionDiffusionStationaryROM
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utilities



class ROMStationaryConvDiff(KratosUnittest.TestCase):
#########################################################################################

    def test_ConvDiff_Stationary_ROM_2D(self):

        with KratosUnittest.WorkFolderScope(".", __file__):
            with open("ProjectParametersROM.json",'r') as parameter_file:
                parameters = KratosMultiphysics.Parameters(parameter_file.read())

            model = KratosMultiphysics.Model()
            simulation = TestConvectionDiffusionStationaryROM(model,parameters)
            simulation.Run()
            ObtainedOutput = simulation.EvaluateQuantityOfInterest()
            ExpectedOutput = np.array([-12.148084271927585, 14.819025485595288, -155.5281654573427, -9.353101257284601, 10.72554942206833, -124.59170838744552, 86.38273570629585, 16.913999272914044, 111.60080675691219, 37.69212146428647, 90.0, 90.0, 232.55768788233212, 253.10649201246736, 90.0, 90.0])
            NodalArea = np.array([0.03703703702962964, 0.05555555555555557, 0.055555555555555566, 0.11111111112222223, 0.05555555555000001, 0.05555555555, 0.11111111112222222, 0.11111111112222222, 0.11111111112222222, 0.01851851851481482, 0.018518518514814813, 0.055555555549999996, 0.055555555549999996, 0.05555555555555555, 0.05555555555555555, 0.03703703702962963])
            L2 = np.sqrt(      (sum(NodalArea*((ExpectedOutput/ExpectedOutput - ObtainedOutput/ExpectedOutput )**2)))  /     (sum(NodalArea))      )*100
            self.assertLess(L2, 0.1) #percent
            # Cleaning
            kratos_utilities.DeleteDirectoryIfExisting("__pycache__")

##########################################################################################

if __name__ == '__main__':
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
    KratosUnittest.main()

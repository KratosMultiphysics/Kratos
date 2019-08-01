from __future__ import print_function, absolute_import, division
import KratosMultiphysics
import numpy as np

import KratosMultiphysics.KratosUnittest as KratosUnittest


class ROMLinearTest(KratosUnittest.TestCase):
#########################################################################################

    def test_Linear_ROM_2D(self):
        from Main_To_Lauch_RomSolverRADIATOR import ConvectionDiffusionAnalysisWithFlush

        with open("ProjectParameters.json",'r') as parameter_file:
            parameters = KratosMultiphysics.Parameters(parameter_file.read())

        model = KratosMultiphysics.Model()
        simulation = ConvectionDiffusionAnalysisWithFlush(model,parameters)
        simulation.Run()
        QoI = simulation.EvaluateQuantityOfInterest()
        ExpectedOutput = [37.69212146428647, 232.5576878823321, 10.725549422068331, 86.38273570629583, 253.10649201246736, 14.819025485595288, 111.60080675691219, -9.353101257284605, 16.91399927291404, 90.0, -12.148084271927583, 90.0, -155.5281654573427, 90.0, -124.59170838744552, 90.0]
        Error = 0
        for i in range(16):
            Error += (QoI[i] - ExpectedOutput[i])**2
        Error = np.sqrt(Error)
        self.assertLess(Error, 1.0e-6)    


##########################################################################################        

if __name__ == '__main__':
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
    KratosUnittest.main()

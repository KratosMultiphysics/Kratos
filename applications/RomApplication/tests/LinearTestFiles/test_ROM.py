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
        ExpectedOutput = [1.9679718937478738, 1.4213130344171971, 1.4213130344171971, 1.0933177188253513, 1.0933177187269525, 1.0933177187269525, 0.7653224031351064, 0.7653224031351064, 0.4373270874448616, 0.9839859468739369, 0.9839859468739369, 0.4373270875432602, 0.4373270875432602, 0.10933177185301551, 0.10933177185301551, 0.0]
        Error = 0
        for i in range(16):
            Error += (QoI[i] - ExpectedOutput[i])**2
        Error = np.sqrt(Error)
        self.assertLess(Error, 1.0e-6)    


##########################################################################################        

if __name__ == '__main__':
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
    KratosUnittest.main()

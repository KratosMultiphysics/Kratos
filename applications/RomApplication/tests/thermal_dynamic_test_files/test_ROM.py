from __future__ import print_function, absolute_import, division
import KratosMultiphysics
import numpy as np

import KratosMultiphysics.KratosUnittest as KratosUnittest
from .MainKratosROM import TestConvectionDiffusionTransientROM
import KratosMultiphysics.kratos_utilities as kratos_utilities

class ROMDynamicConvDiff(KratosUnittest.TestCase):
#########################################################################################

    def test_ConvDiff_Dynamic_ROM_2D(self):

        with KratosUnittest.WorkFolderScope(".", __file__):
            with open("ProjectParametersROM.json",'r') as parameter_file:
                parameters = KratosMultiphysics.Parameters(parameter_file.read())
            model = KratosMultiphysics.Model()      
            Simulation = TestConvectionDiffusionTransientROM(model,parameters)
            Simulation.Run()
            ObtainedOutput = Simulation.EvaluateQuantityOfInterest()
            ExpectedOutput = np.load('ExpectedOutput.npy')
            NodalArea = Simulation.EvaluateQuantityOfInterest2()
            for i in range (np.shape(ObtainedOutput)[1]):
                L2 = np.sqrt(      (sum(NodalArea*((ExpectedOutput[:,i]/ExpectedOutput[:,i] - ObtainedOutput[:,i]/ExpectedOutput[:,i] )**2)))  /     (sum(NodalArea))      )*100
                self.assertLess(L2, 0.1) #percent            
            # Cleaning
            kratos_utilities.DeleteDirectoryIfExisting("__pycache__")



##########################################################################################

if __name__ == '__main__':
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
    KratosUnittest.main()

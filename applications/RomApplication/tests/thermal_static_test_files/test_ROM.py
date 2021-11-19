from __future__ import print_function, absolute_import, division
import KratosMultiphysics

try:
    import numpy as np
    numpy_available = True
    from .MainKratosROM import TestConvectionDiffusionStationaryROM
except:
    numpy_available = False

import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utilities


@KratosUnittest.skipIfApplicationsNotAvailable("ConvectionDiffusionApplication")
class ROMStationaryConvDiff(KratosUnittest.TestCase):
#########################################################################################

    @KratosUnittest.skipUnless(numpy_available, "numpy is required for RomApplication")
    def test_ConvDiff_Stationary_ROM_2D(self):

        with KratosUnittest.WorkFolderScope(".", __file__):
            with open("ProjectParameters.json",'r') as parameter_file:
                parameters = KratosMultiphysics.Parameters(parameter_file.read())

            model = KratosMultiphysics.Model()
            simulation = TestConvectionDiffusionStationaryROM(model,parameters)
            simulation.Run()
            ObtainedOutput = simulation.EvaluateQuantityOfInterest()
            ExpectedOutput = np.load('ExpectedOutput.npy')
            NodalArea = simulation.EvaluateQuantityOfInterest2()
            l2 = np.sqrt(      (sum(NodalArea*((1 - ObtainedOutput/ExpectedOutput )**2)))  /     (sum(NodalArea))      )*100
            self.assertLess(l2, 1e-12) #percent
            # Cleaning
            kratos_utilities.DeleteDirectoryIfExisting("__pycache__")

##########################################################################################

if __name__ == '__main__':
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
    KratosUnittest.main()

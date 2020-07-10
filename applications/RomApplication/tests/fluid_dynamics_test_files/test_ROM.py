from __future__ import print_function, absolute_import, division
import KratosMultiphysics
import numpy as np

from .MainKratosROM import TestFluidDynamicsROM
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utilities



class ROMFluidDynamics(KratosUnittest.TestCase):
#########################################################################################

    def test_Fluid_Dynamic_ROM_2D(self):

        with KratosUnittest.WorkFolderScope(".", __file__):
            with open("ProjectParameters.json",'r') as parameter_file:
                parameters = KratosMultiphysics.Parameters(parameter_file.read())
            model = KratosMultiphysics.Model()
            Simulation = TestFluidDynamicsROM(model,parameters)
            Simulation.Run()
            ObtainedOutput = Simulation.EvaluateQuantityOfInterest()
            ExpectedOutput = np.load('ExpectedOutput.npy')

            for i in range (np.shape(ObtainedOutput)[1]):
                UP=0
                DOWN=0
                UP = sum((ExpectedOutput[:,i] - ObtainedOutput[:,i])**2)
                DOWN = sum((ExpectedOutput[:,i])**2)
                L2 = np.sqrt(UP/DOWN)
                self.assertLess(L2, 1e-12)
            # Cleaning
            kratos_utilities.DeleteDirectoryIfExisting("__pycache__")

##########################################################################################

if __name__ == '__main__':
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
    KratosUnittest.main()

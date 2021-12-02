from __future__ import print_function, absolute_import, division
import KratosMultiphysics

try:
    import numpy as np
    from .MainKratosROM import TestFluidDynamicsROM
    numpy_available = True
except:
    numpy_available = False

import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utilities


@KratosUnittest.skipIfApplicationsNotAvailable("FluidDynamicsApplication")
class ROMFluidDynamics(KratosUnittest.TestCase):
#########################################################################################

    @KratosUnittest.skipUnless(numpy_available, "numpy is required for RomApplication")
    def test_Fluid_Dynamics_ROM_2D(self):

        with KratosUnittest.WorkFolderScope(".", __file__):
            with open("ProjectParameters.json",'r') as parameter_file:
                parameters = KratosMultiphysics.Parameters(parameter_file.read())
            model = KratosMultiphysics.Model()
            Simulation = TestFluidDynamicsROM(model,parameters)
            Simulation.Run()
            ObtainedOutput = Simulation.EvaluateQuantityOfInterest()
            ExpectedOutput = np.load('ExpectedOutput.npy')

            for i in range (np.shape(ObtainedOutput)[1]):
                up = sum((ExpectedOutput[:,i] - ObtainedOutput[:,i])**2)
                down = sum((ExpectedOutput[:,i])**2)
                l2 = np.sqrt(up/down)
                self.assertLess(l2, 1e-12)
            # Cleaning
            kratos_utilities.DeleteDirectoryIfExisting("__pycache__")

##########################################################################################

if __name__ == '__main__':
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
    KratosUnittest.main()

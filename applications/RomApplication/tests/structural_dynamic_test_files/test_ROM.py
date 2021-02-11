from __future__ import print_function, absolute_import, division
import KratosMultiphysics

try:
    import numpy as np
    from .MainKratosROM import TestStructuralMechanicsDynamicROM
    numpy_available = True
except:
    numpy_available = False

import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utilities



class ROMDynamicStruct(KratosUnittest.TestCase):
#########################################################################################

    @KratosUnittest.skipIf(numpy_available == False, "numpy is required for RomApplication")
    def test_Struct_Dynamic_ROM_2D(self):

        with KratosUnittest.WorkFolderScope(".", __file__):
            with open("ProjectParameters.json",'r') as parameter_file:
                parameters = KratosMultiphysics.Parameters(parameter_file.read())
            model = KratosMultiphysics.Model()
            Simulation = TestStructuralMechanicsDynamicROM(model,parameters)
            Simulation.Run()
            ObtainedOutput = Simulation.EvaluateQuantityOfInterest()
            ExpectedOutput = np.load('ExpectedOutput.npy')
            NodalArea = Simulation.EvaluateQuantityOfInterest2()

            for i in range (np.shape(ObtainedOutput)[1]):
                UP = sum((ExpectedOutput[:,i] - ObtainedOutput[:,i])**2)
                DOWN = sum((ExpectedOutput[:,i])**2)
                L2 = np.sqrt(UP/DOWN)
                self.assertLess(L2, 1e-10)
            # Cleaning
            kratos_utilities.DeleteDirectoryIfExisting("__pycache__")

##########################################################################################

if __name__ == '__main__':
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
    KratosUnittest.main()

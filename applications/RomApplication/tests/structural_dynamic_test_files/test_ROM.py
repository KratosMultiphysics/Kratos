from __future__ import print_function, absolute_import, division
import KratosMultiphysics
import numpy as np

from .MainKratosROM import TestStructuralMechanicsDynamicROM
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utilities



class ROMDynamicStruct(KratosUnittest.TestCase):
#########################################################################################

    def test_Struct_Dynamic_ROM_2D(self):

        with KratosUnittest.WorkFolderScope(".", __file__):
            with open("ProjectParametersROM.json",'r') as parameter_file:
                parameters = KratosMultiphysics.Parameters(parameter_file.read())
            model = KratosMultiphysics.Model()      
            Simulation = TestStructuralMechanicsDynamicROM(model,parameters)
            Simulation.Run()
            ObtainedOutput = Simulation.EvaluateQuantityOfInterest()
            ExpectedOutput = np.load('ExpectedOutput.npy')
            NodalArea = Simulation.EvaluateQuantityOfInterest2()

            for i in range (np.shape(ObtainedOutput)[1]):
                UP=0
                DOWN=0
                for j in range((np.shape(ObtainedOutput)[0])):
                    UP += (NodalArea[j]*(    (ExpectedOutput[j,i] - ObtainedOutput[j,i]   )**2)  )
                    DOWN +=  NodalArea[j]
                L2 = np.sqrt(UP/DOWN)
                self.assertLess(L2, 1.0e-4)
            # Cleaning
            kratos_utilities.DeleteDirectoryIfExisting("__pycache__")

##########################################################################################

if __name__ == '__main__':
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
    KratosUnittest.main()

from __future__ import print_function, absolute_import, division
import KratosMultiphysics
import numpy as np

from .MainKratosROM import TestStructuralMechanicsStaticROM
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utilities


class ROMStaticStruct(KratosUnittest.TestCase):
#########################################################################################

    def test_Struct_Static_ROM_2D(self):

        with KratosUnittest.WorkFolderScope(".", __file__):
            with open("ProjectParametersROM.json",'r') as parameter_file:
                parameters = KratosMultiphysics.Parameters(parameter_file.read())
            model = KratosMultiphysics.Model()      
            Simulation = TestStructuralMechanicsStaticROM(model,parameters)
            Simulation.Run()
            ObtainedOutput = Simulation.EvaluateQuantityOfInterest()
            ExpectedOutput = np.load('ExpectedOutput.npy')
            NodalArea = Simulation.EvaluateQuantityOfInterest2()

            UP=0
            DOWN=0
            for i in range (len(ObtainedOutput)):
                UP += (NodalArea[i]*(    (ExpectedOutput[i] - ObtainedOutput[i]   )**2)  )
                DOWN +=  NodalArea[i]
            L2 = np.sqrt(UP/DOWN)
            self.assertLess(L2, 1.0e-4)
            # Cleaning
            kratos_utilities.DeleteDirectoryIfExisting("__pycache__")
            

##########################################################################################

if __name__ == '__main__':
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
    KratosUnittest.main()

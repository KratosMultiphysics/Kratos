from __future__ import print_function, absolute_import, division
import KratosMultiphysics
import numpy as np

import KratosMultiphysics.KratosUnittest as KratosUnittest
from .MainKratosROM import StructStaticROM


class ROMStaticStruct(KratosUnittest.TestCase):
#########################################################################################

    def test_Struct_Static_ROM_2D(self):

        with KratosUnittest.WorkFolderScope(".", __file__):
            with open("ProjectParametersROM.json",'r') as parameter_file:
                parameters = KratosMultiphysics.Parameters(parameter_file.read())
            model = KratosMultiphysics.Model()      
            Simulation = StructStaticROM(model,parameters)
            Simulation.Run()
            ObtainedOutput = Simulation.EvaluateQuantityOfInterest()
            ExpectedOutput = np.load('ExpectedOutput.npy')
            NodalArea = Simulation.EvaluateQuantityOfInterest2()

            UP=0
            DOWN=0
            for i in range (16):
                UP += (NodalArea[i]*(    (ExpectedOutput[i] - ObtainedOutput[i]   )**2)  )
                DOWN +=  NodalArea[i]
            L2 = np.sqrt(UP/DOWN)
            self.assertLess(L2, 1.0e-4)


##########################################################################################

if __name__ == '__main__':
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
    KratosUnittest.main()

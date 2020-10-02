from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics
from KratosMultiphysics import Logger
import KratosMultiphysics.FemToDemApplication as KratosFemDem
# import KratosMultiphysics.FemToDemApplication.MainCouplingFemDem as MainCouplingFemDem
import main_coupling_for_testing
import KratosMultiphysics.KratosUnittest as KratosUnittest
import os
import shutil

def Wait():
    input("Press Something")

def KratosPrintInfo(message):
    KratosMultiphysics.Logger.Print(message, label="")
    KratosMultiphysics.Logger.Flush()

#============================================================================================================================
class MainCouplingFemDemForTestingSolution(main_coupling_for_testing.MainCouplingFemDemForTestingSolution):
#============================================================================================================================

#============================================================================================================================
    def CheckControlValuesForTesting(self):  # KratosPrintInfo(str(dy))
        

        # Here we check the vertical displacement of a node
        node = self.FEM_Solution.main_model_part.GetNode(39)
        dy = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y)
        if self.FEM_Solution.step == 20:
            if dy != -0.046660044881251264:
                raise ValueError('The computed displacement at step = 20 is not correct')
        elif self.FEM_Solution.step == 50:
            if dy != -0.3004937948812545:
                raise ValueError('The computed displacement at step = 50 is not correct')
        elif self.FEM_Solution.step == 90:
            if dy != -0.9778132951993863:
                raise ValueError('The computed displacement at step = 90 is not correct')
        elif self.FEM_Solution.step == 140:
            ref = 0.3262364577734523
            if (dy - ref) / ref > 1e-5:
                raise ValueError('The computed displacement at step = 140 is not correct')
        
        vy = node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_Y)
        if self.FEM_Solution.step == 20:
            if vy != -0.9564750000003245:
                raise ValueError('The computed velocity at step = 20 is not correct')
        elif self.FEM_Solution.step == 50:
            if vy != -2.4279749999996674:
                raise ValueError('The computed velocity at step = 50 is not correct')
        elif self.FEM_Solution.step == 90:
            if vy != -3.3520719444389666:
                raise ValueError('The computed velocity at step = 90 is not correct')
        elif self.FEM_Solution.step == 140:
            ref = 4.7178528813310345
            if (vy - ref) / ref > 1e-5:
                raise ValueError('The computed displacement at step = 140 is not correct')

class TestAnalytics(KratosUnittest.TestCase):
    
    def setUp(self):
        pass

    @classmethod
    def free_fall(self):
        model = KratosMultiphysics.Model()
        MainCouplingFemDemForTestingSolution(model, "small_tests/free_fall_contact/").Run()



if __name__ == "__main__":
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(Logger.Severity.WARNING)
    KratosUnittest.main()
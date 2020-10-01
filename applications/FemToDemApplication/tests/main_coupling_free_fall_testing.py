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
        node = self.FEM_Solution.main_model_part.GetNode(1)
        dy = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y)
        if self.FEM_Solution.step == 26:
            KratosPrintInfo(str(dy))
            # if dy != 1.972140108200514e-05:
            #     raise ValueError('The computed displacement at step = 26 is not correct')
        elif self.FEM_Solution.step == 36:
            KratosPrintInfo(str(dy))
            # if dy != 2.7306779366995683e-05:
            #     raise ValueError('The computed displacement at step = 36 is not correct')
        elif self.FEM_Solution.step == 46:
            KratosPrintInfo(str(dy))
            # if dy != 2.5780167460400344e-05:
            #     raise ValueError('The computed displacement at step = 46 is not correct')
        elif self.FEM_Solution.step == 61:
            KratosPrintInfo(str(dy))
            # if dy != 1.4423580956402616e-05:
            #     raise ValueError('The computed displacement at step = 61 is not correct')
        


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
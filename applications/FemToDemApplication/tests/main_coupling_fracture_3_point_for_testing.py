
import KratosMultiphysics
from KratosMultiphysics import Logger
import KratosMultiphysics.FemToDemApplication as KratosFemDem
import main_coupling_for_testing
import KratosMultiphysics.KratosUnittest as KratosUnittest
import os

def Wait():
    input("Press Something")

def KratosPrintInfo(message):
    KratosMultiphysics.Logger.Print(message, label="")
    KratosMultiphysics.Logger.Flush()

#============================================================================================================================
class MainCouplingFemDemForTestingSolution(main_coupling_for_testing.MainCouplingFemDemForTestingSolution):
#============================================================================================================================

#============================================================================================================================
    def CheckControlValuesForTesting(self):  # KratosPrintInfo(str(damage))

        # Here we check the damage obtained at each FE
        element = self.FEM_Solution.main_model_part.GetElement(301)
        damage = element.CalculateOnIntegrationPoints(KratosFemDem.DAMAGE_ELEMENT, self.FEM_Solution.main_model_part.ProcessInfo)[0]
        tol = 1e-8
        if self.FEM_Solution.step == 15:
            ref_value = 0.8342902589556865
            if abs((damage - ref_value) / ref_value) > tol:
                raise ValueError('The computed damage at step = 15 is not correct')
        elif self.FEM_Solution.step == 20:
            ref_value = 0.8342902589556865
            if abs((damage - ref_value) / ref_value) > tol:
                raise ValueError('The computed damage at step = 20 is not correct')
        elif self.FEM_Solution.step == 25:
            ref_value = 0.8342902589556865
            if abs((damage - ref_value) / ref_value) > tol:
                raise ValueError('The computed damage at step = 25 is not correct')
        elif self.FEM_Solution.step == 28:
            ref_value = 0.8342902589556865
            if abs((damage - ref_value) / ref_value) > tol:
                raise ValueError('The computed damage at step = 28 is not correct')

        # Here we check the vertical displacement of a node
        node = self.FEM_Solution.main_model_part.GetNode(125)
        Ry = node.GetSolutionStepValue(KratosMultiphysics.REACTION_Y)
        if self.FEM_Solution.step == 15:
            ref = -14404.407899045087
            if abs((Ry - ref) / ref) > tol:
                raise ValueError('The computed displacement at step = 15 is not correct')
        elif self.FEM_Solution.step == 20:
            ref = -7449.905801561186
            if abs((Ry - ref) / ref) > tol:
                raise ValueError('The computed displacement at step = 36 is not correct')
        elif self.FEM_Solution.step == 25:
            ref = -5271.126026352192
            if abs((Ry - ref) / ref) > tol:
                raise ValueError('The computed displacement at step = 46 is not correct')
        elif self.FEM_Solution.step == 28:
            ref = -4207.087723385579
            if abs((Ry - ref) / ref) > tol:
                raise ValueError('The computed displacement at step = 61 is not correct')


class TestAnalytics(KratosUnittest.TestCase):

    def setUp(self):
        pass

    @classmethod
    def fracture_3_point(self):
        model = KratosMultiphysics.Model()
        MainCouplingFemDemForTestingSolution(model, os.path.join(os.path.dirname(os.path.realpath(__file__)), "fracture_tests", "three_point_bending")).Run()


if __name__ == "__main__":
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(Logger.Severity.WARNING)
    KratosUnittest.main()
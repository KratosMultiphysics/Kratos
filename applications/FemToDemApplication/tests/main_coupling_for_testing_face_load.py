
import KratosMultiphysics
from KratosMultiphysics import Logger
import KratosMultiphysics.FemToDemApplication as KratosFemDem
# import KratosMultiphysics.FemToDemApplication.MainCouplingFemDem as MainCouplingFemDem
import main_coupling_for_testing
import KratosMultiphysics.KratosUnittest as KratosUnittest

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
        tol = 1e-9
        node = self.FEM_Solution.main_model_part.GetNode(54)
        dy = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y)
        if self.FEM_Solution.step == 1:
            ref = 6.888472011951311e-06
            if abs((dy - ref) / ref) > tol:
                raise ValueError('The computed displacement at step = 1 is not correct ' + str(dy))
        elif self.FEM_Solution.step == 2:
            ref = 6.88849294197138e-06
            if abs((dy - ref) / ref) > tol:
                raise ValueError('The computed displacement at step = 2 is not correct ' + str(dy))
        elif self.FEM_Solution.step == 3:
            ref = 6.888492942036134e-06
            if abs((dy - ref) / ref) > tol:
                raise ValueError('The computed displacement at step = 3 is not correct ' + str(dy))
        elif self.FEM_Solution.step == 4:
            ref = 6.888492942036134e-06
            if abs((dy - ref) / ref) > tol:
                raise ValueError('The computed displacement at step = 4 is not correct ' + str(dy))

        # Here we check the stresses and strains at one FE
        element = self.FEM_Solution.main_model_part.GetElement(1)
        Sx = element.CalculateOnIntegrationPoints(KratosFemDem.STRESS_VECTOR_INTEGRATED,           self.FEM_Solution.main_model_part.ProcessInfo)[0][0]
        Ex = element.CalculateOnIntegrationPoints(KratosMultiphysics.GREEN_LAGRANGE_STRAIN_VECTOR, self.FEM_Solution.main_model_part.ProcessInfo)[0][0]

        if self.FEM_Solution.step == 1:
            refSx = 88714.82305985004
            refEx = 1.976729635484453e-06
            if  abs((Sx - refSx) / refSx) > tol or abs((Ex - refEx) / refEx) > tol:
                raise ValueError('The computed stress/strain at step = 1 is not correct')
        elif self.FEM_Solution.step == 2:
            refSx = 88714.59315099663
            refEx = 1.9767219795336755e-06
            if  abs((Sx - refSx) / refSx) > tol or abs((Ex - refEx) / refEx) > tol:
                raise ValueError('The computed stress/strain at step = 2 is not correct')
        elif self.FEM_Solution.step == 3:
            refSx = 88714.59314896743
            refEx = 1.9767219794759574e-06
            if  abs((Sx - refSx) / refSx) > tol or abs((Ex - refEx) / refEx) > tol:
                raise ValueError('The computed stress/strain at step = 3 is not correct')
        elif self.FEM_Solution.step == 4:
            refSx = 88714.59314896737
            refEx = 1.9767219794759574e-06
            if  abs((Sx - refSx) / refSx) > tol or abs((Ex - refEx) / refEx) > tol:
                raise ValueError('The computed stress/strain at step = 4 is not correct')


class TestAnalytics(KratosUnittest.TestCase):
    
    def setUp(self):
        pass

    @classmethod
    def face_load(self):
        model = KratosMultiphysics.Model()
        MainCouplingFemDemForTestingSolution(model, "small_tests/pressure_face_load/").Run()

if __name__ == "__main__":
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(Logger.Severity.WARNING)
    KratosUnittest.main()
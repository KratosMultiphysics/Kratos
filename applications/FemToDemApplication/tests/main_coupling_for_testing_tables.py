
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

def IsOk(value, ref, tol):
    if abs(value - ref) / ref > tol:
        return False
    return True

#============================================================================================================================
class MainCouplingFemDemForTestingSolution(main_coupling_for_testing.MainCouplingFemDemForTestingSolution):
#============================================================================================================================

#============================================================================================================================
    def CheckControlValuesForTesting(self):  # KratosPrintInfo(str(dy))

        # Here we check the damage obtained at each FE
        tol = 1e-6
        element = self.FEM_Solution.main_model_part.GetElement(1)
        damage = element.CalculateOnIntegrationPoints(KratosFemDem.DAMAGE_ELEMENT, self.FEM_Solution.main_model_part.ProcessInfo)[0]

        if self.FEM_Solution.step == 26:
            if not IsOk(damage, 0.11526580049026725, tol):
                raise ValueError('The computed damage at step = 26 is not correct')
        elif self.FEM_Solution.step == 36:
            if not IsOk(damage, 0.4722648310044538, tol):
                raise ValueError('The computed damage at step = 36 is not correct')
        elif self.FEM_Solution.step == 46:
            if not IsOk(damage, 0.5600214207342531, tol):
                raise ValueError('The computed damage at step = 46 is not correct')
        elif self.FEM_Solution.step == 61:
            if not IsOk(damage, 0.5600214207342531, tol):
                raise ValueError('The computed damage at step = 61 is not correct')

        # Here we check the vertical displacement of a node
        node = self.FEM_Solution.main_model_part.GetNode(1)
        dy = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y)
        if self.FEM_Solution.step == 26:
            if not IsOk(dy, 1.971665114439254e-05, tol):
                raise ValueError('The computed displacement at step = 26 is not correct')
        elif self.FEM_Solution.step == 36:
            if not IsOk(dy, 2.7299978508712653e-05, tol):
                raise ValueError('The computed displacement at step = 36 is not correct')
        elif self.FEM_Solution.step == 46:
            if not IsOk(dy, 2.578331303724928e-05, tol):
                raise ValueError('The computed displacement at step = 46 is not correct')
        elif self.FEM_Solution.step == 61:
            if not IsOk(dy, 1.4408321991404051e-05, tol):
                raise ValueError('The computed displacement at step = 61 is not correct')

        # Here we check the stresses and strains at one FE
        element = self.FEM_Solution.main_model_part.GetElement(1)
        Sx = element.CalculateOnIntegrationPoints(KratosFemDem.STRESS_VECTOR_INTEGRATED,           self.FEM_Solution.main_model_part.ProcessInfo)[0][0]
        Ex = element.CalculateOnIntegrationPoints(KratosMultiphysics.GREEN_LAGRANGE_STRAIN_VECTOR, self.FEM_Solution.main_model_part.ProcessInfo)[0][0]

        if self.FEM_Solution.step == 26:
            if not IsOk(Sx, 1441812.5386046136, tol):
                raise ValueError('The computed stress at step = 26 is not correct')
            if not IsOk(Ex, 4.6561604584527234e-05, tol):
                raise ValueError('The computed strain at step = 26 is not correct')
        elif self.FEM_Solution.step == 36:
            if not IsOk(Sx, 1190806.4343404802, tol):
                raise ValueError('The computed stress at step = 36 is not correct')
            if not IsOk(Ex, 6.446991404011464e-05, tol):
                raise ValueError('The computed strain at step = 36 is not correct')
        elif self.FEM_Solution.step == 46:
            if not IsOk(Sx, 937633.4336071612, tol):
                raise ValueError('The computed stress at step = 46 is not correct')
            if not IsOk(Ex, 6.0888252148997134e-05, tol):
                raise ValueError('The computed strain at step = 46 is not correct')
        elif self.FEM_Solution.step == 61:
            if not IsOk(Sx, 523971.6246628269, tol):
                raise ValueError('The computed stress at step = 61 is not correct')
            if not IsOk(Ex, 3.4025787965616143e-05, tol):
                raise ValueError('The computed strain at step = 61 is not correct')
        


class TestAnalytics(KratosUnittest.TestCase):
    
    def setUp(self):
        pass

    @classmethod
    def tables(self):
        model = KratosMultiphysics.Model()
        MainCouplingFemDemForTestingSolution(model, "small_tests/table_bc_imposition/").Run()



if __name__ == "__main__":
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(Logger.Severity.WARNING)
    KratosUnittest.main()
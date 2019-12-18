import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.CoSimulationApplication.co_simulation_tools import ImportDataStructure
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools

class TestCoupledSolverIBQN(KratosUnittest.TestCase):

    def assertArrayAlmostEqual(self, a1, a2):
        ls1 = list(a1)
        ls2 = list(a2)
        try:
            self.assertEqual(ls1, ls2)
        except AssertionError:
            for i in range(len(ls1)):
                self.assertAlmostEqual(ls1[i], ls2[i])

    def assertArrayEqual(self, a1, a2):
        self.assertEqual(list(a1), list(a2))

    def test_coupled_solver_ibqn(self):
        parameter_file_name = "test_parameters.json"
        cs_data_structure = ImportDataStructure(parameter_file_name)

        parameter_file_name = "coupled_solvers/test_ibqn.json"
        with open(parameter_file_name, 'r') as parameter_file:
            settings = cs_data_structure.Parameters(parameter_file.read())

        coupled_solver = cs_tools.CreateInstance(settings)
        coupled_solver.Initialize()
        coupled_solver.Check()

        coupled_solver.InitializeSolutionStep()
        coupled_solver.SolveSolutionStep()
        sol_x = [0.78539914, 0.78539911, 0.78539908, 0.78539905, 0.78539902, 0.78539899, 0.78539896, 0.78539894,
                 0.78539891, 0.78539888]
        self.assertArrayAlmostEqual(coupled_solver.x.GetNumpyArray(), sol_x)

        coupled_solver.FinalizeSolutionStep()
        coupled_solver.OutputSolutionStep()

        coupled_solver.Finalize()

if __name__ == '__main__':
    KratosUnittest.main()

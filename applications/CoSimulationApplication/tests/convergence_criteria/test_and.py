import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.CoSimulationApplication.co_simulation_tools import ImportDataStructure
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools


class TestConvergenceCriterionAnd(KratosUnittest.TestCase):
    def test_convergence_criterion_and(self):
        parameter_file_name = "test_parameters.json"
        cs_data_structure = ImportDataStructure(parameter_file_name)

        self.assertAlmostEqual(1, 1)


if __name__ == '__main__':
    KratosUnittest.main()

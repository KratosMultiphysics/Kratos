import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as kratos_unittest

from KratosMultiphysics.kratos_utilities import DeleteFileIfExisting
from KratosMultiphysics.OptimizationApplication.optimization_analysis import OptimizationAnalysis
import csv, os

try:
    import nlopt
    nlopt_available = True
except ImportError:
    nlopt_available = False

@kratos_unittest.skipIf(not nlopt_available, "Missing nlopt python libraries ")
class TestNLOPTOptimizers(kratos_unittest.TestCase):
    # @classmethod
    # def setUpClass(cls):
    #     with kratos_unittest.WorkFolderScope(".", __file__):
    #         with open("optimization_parameters.json", "r") as file_input:
    #             cls.parameters = Kratos.Parameters(file_input.read())
    #         cls.model = Kratos.Model()
    #         cls.analysis = OptimizationAnalysis(cls.model, cls.parameters)

    def test_MMA_optimizer(self):
        pass
        with kratos_unittest.WorkFolderScope(".", __file__):
            with open("optimization_parameters.json", "r") as file_input:
                parameters = Kratos.Parameters(file_input.read())
            model = Kratos.Model()
            analysis = OptimizationAnalysis(model, parameters)
            analysis.Run()
            optimization_log_filename = "summary.csv"
            # Testing
            with open(optimization_log_filename, 'r') as csvfile:
                reader = csv.reader(csvfile, delimiter=',')
                last_line = None
                for i, row in enumerate(reader):
                    if i == 15:
                        last_line = row

                resulting_strain_energy = float(last_line[1].strip())
                resulting_mass = float(last_line[2].strip())

                self.assertAlmostEqual(resulting_mass, 7642.774084393494,2)
                self.assertAlmostEqual(resulting_strain_energy, 5.209266637299775e-05,2)

            os.remove(optimization_log_filename)

if __name__ == "__main__":
    kratos_unittest.main()

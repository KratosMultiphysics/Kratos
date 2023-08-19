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
class TestNLOPTOptimizer(kratos_unittest.TestCase):
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

                resulting_mass = float(last_line[1].strip())
                resulting_strain_energy = float(last_line[4].strip())

                self.assertAlmostEqual(resulting_mass, 7278.149435,2)
                self.assertAlmostEqual(resulting_strain_energy, 5.572267523e-5,2)

            os.remove(optimization_log_filename)

if __name__ == "__main__":
    kratos_unittest.main()

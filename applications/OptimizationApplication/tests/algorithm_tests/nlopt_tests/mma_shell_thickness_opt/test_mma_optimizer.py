import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as kratos_unittest

from KratosMultiphysics.OptimizationApplication.optimization_analysis import OptimizationAnalysis
import csv, os

try:
    import nlopt
except ImportError:
    nlopt = None

@kratos_unittest.skipIf(nlopt==None, "Missing nlopt python libraries ")
class TestNLOPTOptimizer(kratos_unittest.TestCase):
    def test_MMA_optimizer(self):
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
                    if i == 13:
                        last_line = row

                resulting_mass = float(last_line[1].strip())
                self.assertAlmostEqual(resulting_mass, 542.0733237,2)

            os.remove(optimization_log_filename)

if __name__ == "__main__":
    kratos_unittest.main()

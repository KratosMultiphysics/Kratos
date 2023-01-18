import KratosMultiphysics as Kratos

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as kratos_unittest

from KratosMultiphysics.OptimizationApplication.optimization_analysis import (
    OptimizationAnalysis,
)


class TestMassOptimization(kratos_unittest.TestCase):
    def test_SystemIdentification(self):
        model = Kratos.Model()
        with kratos_unittest.WorkFolderScope("system_identification", __file__):
            with open("optimization_parameters.json", "r") as file_input:
                parameters = Kratos.Parameters(file_input.read())

            analysis = OptimizationAnalysis(model, parameters)
            analysis.Run()


if __name__ == "__main__":
    Kratos.Tester.SetVerbosity(Kratos.Tester.Verbosity.PROGRESS)  # TESTS_OUTPUTS
    kratos_unittest.main()
